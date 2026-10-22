#!/usr/bin/env python3
"""
find_tdna_insertion.py
======================
Locate T-DNA / transgene insertion sites in a plant genome from 3'-mRNA-seq
(or any RNA-seq / DNA-seq) short-read data, given the transformation plasmid.

Principle
---------
A T-DNA inserted in the genome produces chimeric transcripts / reads that
contain BOTH T-DNA sequence and plant genomic sequence across a junction.
Such a read is (i) chimeric -- it does not map end-to-end to either the plant
genome or the plasmid -- and (ii) its T-DNA moiety abuts one of the T-DNA
borders (LB / RB repeat).  The plant moiety is then mapped to the genome, and
the mapping position is the insertion site.

Pipeline
--------
 1. parse plasmid (GenBank / FASTA); find LB & RB border repeats; make a
    circularly doubled reference so reads across the plasmid origin align.
 2. align all R1 reads to the doubled plasmid (minimap2 via mappy);
    keep reads that contain >= --min-plasmid-part bp of plasmid.
 3. split every hit into its plasmid moiety and its foreign ("tail") moiety,
    in a plasmid-forward frame; clean sequencing artefacts
    (template-switch GGG, polyA/polyT, Illumina adapter read-through,
    low complexity).
 4. map tails to the genome; map the FULL reads to the genome as well.
    Keep only junctions that are genuinely chimeric (the full read does
    NOT map end-to-end to the genome) and whose T-DNA coordinate lies
    within --border-window bp of a T-DNA border.
 5. cluster the surviving junctions by genomic position -> insertion sites.
 6. refine each breakpoint to base resolution by an exact-substring grid
    search over (plant position x T-DNA position), which also detects
    microhomology at the junction.
 7. build a locus contig (plant flank + integrated T-DNA + plant flank),
    realign all reads to it and report the number of reads spanning each
    junction, plus T-DNA / flank coverage, as validation.
 8. annotate (optionally) with the nearest genes from a GFF3.
 9. write reports + locus FASTA (ready for primer design / Sanger validation).

Requires: mappy (pip install mappy), pysam (pip install pysam)

Example
-------
  python find_tdna_insertion.py \
      --r1 C110M_UMI_R1.fastq.gz \
      --plasmid pTF101-dsRED_clv3.gb \
      --genome Gmax_880_v6.0.fa \
      --gff Gmax_880_Wm82.a6.v1.gene_exons.gff3 \
      --outdir out_C110M --threads 8

Author: generated for the soybean clv3 / miR-clv3 T-DNA line.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
import time

try:
    import mappy
except ImportError:  # pragma: no cover
    sys.exit("ERROR: mappy is required  (pip install mappy)")
try:
    import pysam
except ImportError:  # pragma: no cover
    sys.exit("ERROR: pysam is required  (pip install pysam)")

# ----------------------------------------------------------------------------
# small sequence helpers
# ----------------------------------------------------------------------------
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def wrap(s: str, n: int = 60) -> str:
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def kmer_counter(seq: str, k: int) -> dict:
    c = collections.Counter()
    for i in range(len(seq) - k + 1):
        c[seq[i:i + k]] += 1
    return c


ADAPTER_K = [
    "AGATCGGAAGAGC", "CTGTCTCTTATACACATCT", "TGGAATTCTCGGGTGCCAAGG",
    "GTTCAGAGTTCTACAGTCCGACGATC", "AGGTGGTTCTAC", "GTAGAACCACCT",
    "AAGCAGTGGTATCAACGCAGAGTAC", "AAAAAAAAAAAAAAAAAAAA",
]


def looks_like_adapter(s: str) -> bool:
    return any(a in s for a in ADAPTER_K)


def clean_tail(x: str) -> str:
    """remove template-switch G, poly(A/T), adapter remnants, from a read tail"""
    y = re.sub(r"^G{3,}", "", x)
    y = re.sub(r"[AT]{6,}$", "", y)
    y = re.sub(r"^[AT]{6,}", "", y)
    y = re.sub(r"G{4,}$", "", y)          # template-switch G on the other side
    return y


def low_complexity(s: str, max_frac: float = 0.70, k: int = 3) -> bool:
    if len(s) < 12:
        return True
    if len(set(s)) < 3:
        return True
    c = collections.Counter(s)
    return c.most_common(1)[0][1] / len(s) > max_frac


# ----------------------------------------------------------------------------
# plasmid handling
# ----------------------------------------------------------------------------
class Plasmid:
    """plasmid sequence + feature table + T-DNA borders"""

    def __init__(self, path: str):
        self.path = path
        self.seq = ""
        self.features = []           # (type, start1, end1, strand, label, raw_quals)
        if path.lower().endswith((".gb", ".gbk", ".genbank", ".gbff")):
            self._parse_genbank()
        else:
            self._parse_fasta()
        self.L = len(self.seq)
        if self.L == 0:
            raise ValueError("empty plasmid sequence")
        self.borders = self._find_borders()
        self.tdna = self._tdna_bounds()

    # -- parsers ------------------------------------------------------------
    def _parse_fasta(self):
        self.seq = "".join(l.strip() for l in open(self.path) if not l.startswith(">")).upper()

    def _parse_genbank(self):
        buf, in_origin, cur = [], False, None
        for line in open(self.path):
            if line.startswith("LOCUS"):
                m = re.search(r"(\d+)\s*bp", line)
                if m:
                    pass
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if line.startswith("//"):
                in_origin = False
                break
            if in_origin:
                buf.append(re.sub(r"[^acgtnACGTN]", "", line))
                continue
            if line.startswith("FEATURES"):
                continue
            m = re.match(r"^\s{5}(\S+)\s+(\S+)", line)
            if m and not line.startswith("                     /"):
                if cur:
                    self.features.append(cur)
                ftype, loc = m.group(1), m.group(2)
                cur = [ftype, loc, None, None, None, {}]
                continue
            if cur is not None:
                q = re.match(r'^\s+/(\w+)=(.*)$', line.rstrip())
                if q:
                    cur[5][q.group(1)] = q.group(2).strip('"')
                elif line.startswith("                     /"):
                    pass
        if cur:
            self.features.append(cur)
        self.seq = "".join(buf).upper()

        parsed = []
        for ftype, loc, _, _, _, quals in self.features:
            strand = "+"
            if loc.startswith("complement"):
                strand = "-"
            nums = re.findall(r"(\d+)\.\.(\d+)", loc)
            if not nums:
                continue
            s, e = int(nums[0][0]), int(nums[0][1])
            label = quals.get("label") or quals.get("gene") or quals.get("product") or ""
            parsed.append((ftype, s, e, strand, label, quals))
        self.features = parsed

    # -- borders ------------------------------------------------------------
    def _find_borders(self):
        """locate LB / RB border repeats, either from feature labels or by motif"""
        lb = rb = None
        for ftype, s, e, strand, label, quals in self.features:
            lab = label.lower()
            if "left border" in lab or re.search(r"\bLB\b", label):
                lb = (s, e, label)
            if "right border" in lab or re.search(r"\bRB\b", label):
                rb = (s, e, label)
        if lb is None or rb is None:
            # fall back to the conserved border repeat motif
            core = "TGGCAGGATATATTGTGGTGTAAAC"
            core2 = revcomp(core)
            hits = [m.start() + 1 for m in re.finditer(core, self.seq)]
            hits += [m.start() + 1 for m in re.finditer(core2, self.seq)]
            hits.sort()
            if len(hits) >= 2:
                if lb is None:
                    lb = (hits[0], hits[0] + len(core) - 1, "LB(repeat-motif)")
                if rb is None:
                    rb = (hits[-1], hits[-1] + len(core) - 1, "RB(repeat-motif)")
        return {"LB": lb, "RB": rb}

    def _tdna_bounds(self):
        """integrated region of interest = from LB start to RB end (inclusive)"""
        lb, rb = self.borders["LB"], self.borders["RB"]
        if lb and rb:
            s = min(lb[0], rb[0])
            e = max(lb[1], rb[1])
            return (s, e)
        return (1, self.L)

    # -- convenience --------------------------------------------------------
    def sub(self, a: int, b: int) -> str:
        """1-based inclusive, supports wrap-around"""
        a -= 1
        if b <= self.L and a >= 0 and b > a:
            return self.seq[a:b]
        out = []
        for p in range(a, b):
            out.append(self.seq[p % self.L])
        return "".join(out)

    def write_fasta(self, path: str):
        with open(path, "w") as f:
            f.write(">plasmid_%s circular_%dbp\n" % (os.path.basename(self.path), self.L))
            f.write(wrap(self.seq) + "\n")
            f.write(">plasmid_2x_doubled\n")
            f.write(wrap(self.seq + self.seq) + "\n")

    def feature_table(self) -> list:
        return [{"type": t, "start": s, "end": e, "strand": st, "label": lab}
                for (t, s, e, st, lab, q) in self.features]


# ----------------------------------------------------------------------------
# read loading & plasmid screening
# ----------------------------------------------------------------------------
def open_fastq(path: str):
    import gzip
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def screen(reads_path: str, plasmid: Plasmid, threads: int,
           min_plasmid_part: int = 30, progress: int = 20_000_000):
    """align every read to the doubled plasmid; return hit records"""
    idx_path = None
    aligner = mappy.Aligner(seq=plasmid.seq + plasmid.seq, preset="sr", n_threads=threads)
    if aligner is None:
        raise RuntimeError("could not build plasmid index")
    hits = []
    total = 0
    t0 = time.time()
    with open_fastq(reads_path) as f:
        name = None
        for i, line in enumerate(f):
            m = i & 3
            if m == 0:
                name = line.split()[0][1:]
            elif m == 1:
                s = line.strip().rstrip("N")
                total += 1
                if len(s) < 40:
                    continue
                best = None
                for h in aligner.map(s):
                    if h.is_primary and (best is None or h.mlen > best.mlen):
                        best = h
                if best is not None and best.mlen >= min_plasmid_part:
                    hits.append({"name": name, "seq": s,
                                 "r_st": best.r_st, "r_en": best.r_en,
                                 "q_st": best.q_st, "q_en": best.q_en,
                                 "mlen": best.mlen, "blen": best.blen,
                                 "NM": best.NM, "strand": best.strand,
                                 "cigar": best.cigar_str, "mapq": best.mapq})
                if progress and total % progress == 0:
                    print(f"    ... {total:,} reads scanned, {len(hits):,} plasmid hits "
                          f"({time.time()-t0:.0f}s)", flush=True)
    print(f"  scanned {total:,} reads; {len(hits):,} contain plasmid sequence "
          f"({100.0*len(hits)/max(1,total):.4f}%)  [{time.time()-t0:.0f}s]")
    return total, hits


# ----------------------------------------------------------------------------
# split hits into plasmid / foreign moieties
# ----------------------------------------------------------------------------
def split_hits(hits, plasmid: Plasmid, border_window: int,
               min_plant_tail: int = 25):
    """Return candidate junction records (one per read end that has a tail)."""
    L = plasmid.L
    cand = []
    for h in hits:
        s = len(h["seq"])
        if h["strand"] == 1:
            T = h["seq"]
            qs, qe = h["q_st"], h["q_en"]
        else:
            T = revcomp(h["seq"])
            qs, qe = s - h["q_en"], s - h["q_st"]
        pL = h["r_st"] % L + 1        # 1-based plasmid coord of first aligned base
        pR = (h["r_en"] - 1) % L + 1  # 1-based plasmid coord of last aligned base
        tails = []
        if qs >= 20:
            tails.append(("up", T[:qs], pL, "minus"))    # extends to lower plasmid coords
        if s - qe >= 20:
            tails.append(("down", T[qe:], pR, "plus"))   # extends to higher plasmid coords
        for direction, raw, pcoord, plas_dir in tails:
            tail = clean_tail(raw)
            if len(tail) < min_plant_tail:
                continue
            if low_complexity(tail) or looks_like_adapter(raw):
                continue
            # T-DNA border proximity (allow either side of the border)
            near = None
            for bname in ("LB", "RB"):
                b = plasmid.borders[bname]
                if b is None:
                    continue
                if min(b[0], pcoord) - border_window <= max(b[1], pcoord) \
                   and min(b[0], pcoord) <= max(b[1], pcoord) + border_window:
                    d = 0 if b[0] <= pcoord <= b[1] else (
                        b[0] - pcoord if pcoord < b[0] else pcoord - b[1])
                    if d <= border_window and (near is None or d < near[1]):
                        near = (bname, d)
            if near is None:
                continue
            cand.append({"name": h["name"], "read": h["seq"], "T": T,
                         "strand_plasmid": h["strand"], "tail": tail,
                         "raw_tail": raw, "plasmid_coord": pcoord,
                         "plasmid_dir": plas_dir, "tail_dir": direction,
                         "border": near[0], "border_dist": near[1],
                         "plasmid_part": h["mlen"], "read_len": s})
    print(f"  {len(cand):,} border-proximal candidate junctions "
          f"(>= {min_plant_tail} bp plant tail) out of {len(hits):,} plasmid reads")
    return cand


# ----------------------------------------------------------------------------
# genome mapping of tails & chimerism filter
# ----------------------------------------------------------------------------
def map_to_genome(cand, genome_fa: str, threads: int, min_tail_map: int = 25):
    aligner = mappy.Aligner(genome_fa, preset="sr", n_threads=threads)
    if aligner is None:
        raise RuntimeError("could not build genome index")
    kept = []
    for c in cand:
        best = None
        for h in aligner.map(c["tail"]):
            if h.is_primary and (best is None or h.mlen > best.mlen):
                best = h
        if best is None or best.mlen < min_tail_map:
            continue
        # genomic coordinate of the base that abuts the T-DNA
        if c["tail_dir"] == "down":
            t = best.q_st                      # tail 5' end abuts T-DNA
        else:
            t = best.q_en - 1                  # tail 3' end abuts T-DNA
        gpos = best.r_st + t if best.strand == 1 else best.r_en - 1 - t
        c.update({"chrom": best.ctg, "gpos0": gpos, "tail_strand": best.strand,
                  "tail_mlen": best.mlen, "tail_NM": best.NM,
                  "tail_r_st": best.r_st, "tail_r_en": best.r_en})
        # chimerism: the full read must NOT map end-to-end to the genome
        full = None
        for h in aligner.map(c["read"]):
            if h.is_primary and (full is None or h.mlen > full.mlen):
                full = h
        c["full_genome_mlen"] = full.mlen if full else 0
        if full is not None and full.mlen >= 0.90 * c["read_len"]:
            continue                            # whole read is genomic -> artefact
        kept.append(c)
    print(f"  {len(kept):,} junctions map to the genome and are genuinely chimeric")
    return kept


# ----------------------------------------------------------------------------
# clustering
# ----------------------------------------------------------------------------
def cluster_candidates(kept, slack: int = 250):
    groups = collections.defaultdict(list)
    for c in kept:
        groups[c["chrom"]].append(c)
    sites = []
    for chrom, cs in groups.items():
        cs.sort(key=lambda x: x["gpos0"])
        cur = []
        last = None
        for c in cs:
            if last is not None and c["gpos0"] - last > slack:
                sites.append(cur)
                cur = []
            cur.append(c)
            last = c["gpos0"]
        if cur:
            sites.append(cur)
    out = []
    for site in sites:
        names = {c["name"] for c in site}
        out.append({"chrom": site[0]["chrom"],
                    "start0": min(c["gpos0"] for c in site),
                    "end0": max(c["gpos0"] for c in site),
                    "n_reads": len(names), "junction_records": site})
    out.sort(key=lambda x: -x["n_reads"])
    return out


# ----------------------------------------------------------------------------
# base-resolution breakpoint refinement
# ----------------------------------------------------------------------------
def refine(site, plasmid: Plasmid, genome: pysam.FastaFile, pad: int = 12,
           max_reads: int = 400, win: int = 200):
    """Base-resolution breakpoint search.

    For every candidate junction orientation we build the junction reference
    (200 bp plant | 200 bp T-DNA) and grid-search the exact split point that
    maximises the number of reads which are an exact substring of it.  All
    (plant, T-DNA) pairs giving the same maximum support are equivalent
    solutions (i.e. they sit inside a stretch of microhomology at the
    junction); they are reported together.
    """
    L = plasmid.L
    d = plasmid.seq * 3                      # triple copy: safe indices at the origin
    chrom = site["chrom"]
    clen = genome.get_reference_length(chrom)
    lo = max(1, site["start0"] - win - pad)  # 1-based inclusive
    hi = min(clen, site["end0"] + win + pad)
    G = genome.fetch(chrom, lo - 1, hi)

    def gseq(p, n):
        if p < lo or p + n - 1 > hi:
            return None
        return G[p - lo:p - lo + n]

    def variants(T):
        out = []
        for a in range(0, 5):
            for b in range(0, 21):
                x = T[a:len(T) - b] if b else T[a:]
                if len(x) >= 100:
                    out.append(x)
        return out

    refined = []
    for order in ("up", "down"):
        recs = [c for c in site["junction_records"]
                if (c["tail_dir"] == "down") == (order == "down")]
        if not recs:
            continue
        # one read per (name, direction), prefer the longest plant tail
        uniq = {}
        for c in recs:
            uniq.setdefault(c["name"], c)
        recs = sorted(uniq.values(), key=lambda c: -len(c["tail"]))[:max_reads]
        var_cache = {c["name"]: variants(c["T"]) for c in recs}
        p0 = sorted(c["plasmid_coord"] for c in recs)[len(recs) // 2]
        U0 = sorted(c["gpos0"] + 1 for c in recs)[len(recs) // 2]

        def count(p, U):
            if order == "down":
                pl_part = d[L + p - win:L + p]
                g_part = gseq(U, win)
                J = None if g_part is None else pl_part + g_part
            else:
                pl_part = d[L + p - 1:L + p - 1 + win]
                g_part = gseq(U - win + 1, win)
                J = None if g_part is None else g_part + pl_part
            if not J:
                return 0, 0
            n = tot = 0
            for vs in var_cache.values():
                best = 0
                for v in vs:
                    if len(v) > best and v in J:
                        best = len(v)
                if best >= 100:
                    n += 1
                    tot += best
            return n, tot

        scored = []
        for p in range(p0 - pad, p0 + pad + 1):
            pp = ((p - 1) % L) + 1
            for U in range(U0 - pad, U0 + pad + 1):
                n, tot = count(pp, U)
                if n:
                    scored.append((n, tot, pp, U))
        if not scored:
            continue
        scored.sort(reverse=True)
        best_n = scored[0][0]
        eq = [x for x in scored if x[0] == best_n and x[1] == scored[0][1]]
        # keep only the contiguous diagonal of equivalent solutions
        eq.sort(key=lambda x: x[3])
        refined.append({"order": order, "support": best_n, "bp": scored[0][1],
                        "plasmid_pos": scored[0][2], "plant_pos": scored[0][3],
                        "plant_positions": [x[3] for x in eq],
                        "plasmid_positions": [x[2] for x in eq],
                        "n_ties": len(eq)})
    refined.sort(key=lambda r: -r["support"])
    site["refined"] = refined
    return site


# ----------------------------------------------------------------------------
# locus construction & validation
# ----------------------------------------------------------------------------
def build_locus(site, plasmid: Plasmid, genome: pysam.FastaFile,
                flank: int = 1500):
    """build plant|T-DNA|plant contig from the best refined junction pair"""
    if not site.get("refined"):
        return None
    chrom = site["chrom"]
    clen = genome.get_reference_length(chrom)
    ups = [r for r in site["refined"] if r["order"] == "up"]
    downs = [r for r in site["refined"] if r["order"] == "down"]
    if not ups and not downs:
        return None

    # ---- one-sided evidence: RNA-seq very often shows only one border ----
    one_sided = None
    if downs and not ups:
        one_sided = "LB"                     # LB junction not observed
        d0 = downs[0]
        u_right = min(d0["plant_positions"])
        t_end = d0["plasmid_positions"][d0["plant_positions"].index(u_right)]
        u_left = u_right - 1                 # other flank assumed adjacent
        t_start = plasmid.tdna[0]            # inferred from the annotation
        if t_end <= t_start:
            t_end += plasmid.L
        tdna = plasmid.sub(t_start, t_end)
        deleted = ""
        del_options = None
        lb_pos_span = (u_left, u_left)
        rb_pos_span = (min(d0["plant_positions"]), max(d0["plant_positions"]))
        lb_sup, rb_sup = 0, d0["support"]
        lflank = genome.fetch(chrom, max(0, u_left - flank), u_left)
        rflank = genome.fetch(chrom, u_right - 1, min(clen, u_right - 1 + flank))
        return {"chrom": chrom, "u_left": u_left, "u_right": u_right,
                "one_sided": one_sided,
                "lb_pos_span": lb_pos_span, "rb_pos_span": rb_pos_span,
                "lb_microhomology_bp": 1, "rb_microhomology_bp":
                    rb_pos_span[1] - rb_pos_span[0] + 1,
                "alternative_lb": [], "del_options": del_options,
                "tdna_start": t_start, "tdna_end": (t_end - 1) % plasmid.L + 1,
                "tdna_len": len(tdna), "deleted_plant": deleted,
                "contig": lflank + tdna + rflank,
                "lflank_len": len(lflank), "rflank_len": len(rflank),
                "lb_support": lb_sup, "rb_support": rb_sup,
                "junction_LB": (lflank[-30:] + "|" + tdna[:30]) if one_sided == "LB"
                               else (lflank[-30:] + "|" + tdna[:30]),
                "junction_RB": tdna[-30:] + "|" + rflank[:30]}
    if ups and not downs:
        one_sided = "RB"                     # RB junction not observed
        u0 = ups[0]
        u_left = max(u0["plant_positions"])
        t_start = u0["plasmid_positions"][u0["plant_positions"].index(u_left)]
        u_right = u_left + 1
        t_end = plasmid.tdna[1]              # inferred from the annotation
        if t_end < t_start:
            t_end += plasmid.L
        tdna = plasmid.sub(t_start, t_end)
        deleted = ""
        del_options = None
        lb_pos_span = (min(u0["plant_positions"]), max(u0["plant_positions"]))
        rb_pos_span = (u_right, u_right)
        lflank = genome.fetch(chrom, max(0, u_left - flank), u_left)
        rflank = genome.fetch(chrom, u_right - 1, min(clen, u_right - 1 + flank))
        return {"chrom": chrom, "u_left": u_left, "u_right": u_right,
                "one_sided": one_sided,
                "lb_pos_span": lb_pos_span, "rb_pos_span": rb_pos_span,
                "lb_microhomology_bp": lb_pos_span[1] - lb_pos_span[0] + 1,
                "rb_microhomology_bp": 1,
                "alternative_lb": u0["plant_positions"], "del_options": del_options,
                "tdna_start": t_start, "tdna_end": (t_end - 1) % plasmid.L + 1,
                "tdna_len": len(tdna), "deleted_plant": deleted,
                "contig": lflank + tdna + rflank,
                "lflank_len": len(lflank), "rflank_len": len(rflank),
                "lb_support": u0["support"], "rb_support": 0,
                "junction_LB": lflank[-30:] + "|" + tdna[:30],
                "junction_RB": tdna[-30:] + "|" + rflank[:30]}
    # Each junction has an equivalence class of (plant, T-DNA) positions when the
    # breakpoint sits inside a stretch of microhomology.  Choose the assignment
    # that (a) is internally consistent, u_left < u_right, and (b) loses the
    # fewest plant bases.
    max_del = 50
    up, down = ups[0], downs[0]
    u_right = min(down["plant_positions"])          # earliest plant base after the T-DNA
    cand = [p for p in up["plant_positions"] if u_right - max_del - 1 <= p < u_right]
    if not cand:
        return None
    u_left = max(cand)                              # latest plant base kept before the T-DNA
    t_start = up["plasmid_positions"][up["plant_positions"].index(u_left)]
    t_end = down["plasmid_positions"][down["plant_positions"].index(u_right)]
    if t_end < t_start:
        t_end += plasmid.L
    tdna = plasmid.sub(t_start, t_end)
    lflank = genome.fetch(chrom, max(0, u_left - flank), u_left)
    rflank = genome.fetch(chrom, u_right - 1, min(clen, u_right - 1 + flank))
    contig = lflank + tdna + rflank
    deleted = genome.fetch(chrom, u_left, u_right - 1) if u_right > u_left + 1 else ""
    del_options = sorted({u_right - p - 1 for p in up["plant_positions"] if p < u_right})
    # microhomology span at each junction
    lb_pos_span = (min(up["plant_positions"]), max(up["plant_positions"]))
    rb_pos_span = (min(down["plant_positions"]), max(down["plant_positions"]))
    return {"chrom": chrom, "u_left": u_left, "u_right": u_right,
            "one_sided": one_sided,
            "lb_pos_span": lb_pos_span, "rb_pos_span": rb_pos_span,
            "lb_microhomology_bp": lb_pos_span[1] - lb_pos_span[0] + 1,
            "rb_microhomology_bp": rb_pos_span[1] - rb_pos_span[0] + 1,
            "alternative_lb": up["plant_positions"],
            "del_options": del_options,
            "tdna_start": t_start, "tdna_end": (t_end - 1) % plasmid.L + 1,
            "tdna_len": len(tdna), "deleted_plant": deleted,
            "contig": contig, "lflank_len": len(lflank), "rflank_len": len(rflank),
            "lb_support": up["support"], "rb_support": down["support"],
            "junction_LB": lflank[-30:] + "|" + tdna[:30],
            "junction_RB": tdna[-30:] + "|" + rflank[:30]}


def _al_from_records(records):
    """build a mappy index on a small in-memory multi-record fasta"""
    import tempfile
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    for name, seq in records:
        tmp.write(f">{name}\n{seq}\n")
    tmp.close()
    al = mappy.Aligner(tmp.name, preset="sr")
    os.unlink(tmp.name)
    return al


def validate_locus(locus, reads_path: str, threads: int, genome=None, flank: int = 300):
    """Realign every read to two small references in one pass:

      (1) the reconstructed transgene locus  -> read support per junction, coverage;
      (2) the intact WILD-TYPE sequence across the same insertion point, plus two
          length-matched control windows ~50 kb away -> is the site transcribed at all,
          and are there transcripts running *through* the undisrupted site?

    (2) is the only allele-dosage signal an RNA-seq library can offer, and it is weak:
    the absence of wild-type transcripts can mean homozygosity OR simply that the site
    is not transcribed.  Two separate aligners are used on purpose -- with a shared
    reference, reads from the flanks are reported as primary for one record only and
    the depth of the other record is starved.
    """
    al = _al_from_records([("locus", locus["contig"])])

    wt_recs, wt_len, wt_junction = [], None, None
    if genome is not None:
        chrom, ul, ur = locus["chrom"], locus["u_left"], locus["u_right"]
        if ur > ul:
            wt = genome.fetch(chrom, max(0, ul - 1 - flank), ur + flank)
            wt_recs.append(("wt", wt))
            wt_len, wt_junction = len(wt), flank
            clen = genome.get_reference_length(chrom)
            for tag, off in (("ctrl_up", -50000), ("ctrl_down", 50000)):
                a = ul - 1 + off - flank // 2
                if 0 <= a and a + wt_len <= clen:
                    wt_recs.append((tag, genome.fetch(chrom, a, a + wt_len)))
    al_wt = _al_from_records(wt_recs) if wt_recs else None

    lb_pos = locus["lflank_len"]
    rb_pos = locus["lflank_len"] + locus["tdna_len"]
    st = collections.Counter()
    depth = [0] * (len(locus["contig"]) + 1)
    wt_depth = [0] * (wt_len + 1) if wt_len else None
    ctl_depth = {}

    with open_fastq(reads_path) as f:
        for i, line in enumerate(f):
            if (i & 3) != 1:
                continue
            s = line.strip()
            best = None
            for h in al.map(s):
                if h.is_primary and (best is None or h.mlen > best.mlen):
                    best = h
            if best is not None and best.mlen >= 0.85 * len(s) and best.NM <= 8:
                st["locus_hits"] += 1
                slb = best.r_st < lb_pos < best.r_en
                srb = best.r_st < rb_pos < best.r_en
                if slb:
                    st["span_LB"] += 1
                if srb:
                    st["span_RB"] += 1
                if not slb and not srb:
                    for p in range(best.r_st, best.r_en):
                        depth[p] += 1

            if al_wt is None:
                continue
            bw = None
            for h in al_wt.map(s):
                if h.is_primary and (bw is None or h.mlen > bw.mlen):
                    bw = h
            if bw is None or bw.NM > 3 or bw.mlen < 0.95 * len(s):
                continue
            if bw.ctg == "wt":
                if bw.r_st <= wt_junction - 1 and bw.r_en >= wt_junction + 1:
                    st["wt_spanning"] += 1
                    continue
                for p in range(bw.r_st, bw.r_en):
                    wt_depth[p] += 1
            else:
                ctl_depth.setdefault(bw.ctg, [0] * (wt_len + 1))
                for p in range(bw.r_st, bw.r_en):
                    ctl_depth[bw.ctg][p] += 1

    # a one-sided locus has an *assumed* junction on the missing side; reads crossing
    # that assumed point carry no information, so do not report a support number
    one_sided = locus.get("one_sided")
    locus["support_span_LB"] = None if one_sided == "LB" else st["span_LB"]
    locus["support_span_RB"] = None if one_sided == "RB" else st["span_RB"]
    locus["n_good_alignments"] = st["locus_hits"]

    def mean(arr, a, b):
        a = max(0, a); b = min(len(arr) - 1, b)
        v = arr[a:b]
        return round(sum(v) / max(1, len(v)), 2)

    locus["depth_TDNA"] = mean(depth, lb_pos + 200, rb_pos - 200)
    locus["depth_left_flank"] = mean(depth, 0, lb_pos - 200)
    locus["depth_right_flank"] = mean(depth, rb_pos + 200, len(depth) - 1)

    # --------------- allele-dosage / zygosity diagnostic (weak by design) ---------------
    if wt_depth is not None:
        left_d = mean(wt_depth, 0, wt_junction - 20)
        right_d = mean(wt_depth, wt_junction + 21, wt_len - 1)
        local = max(left_d, right_d)
        bg = max([mean(v, 0, wt_len - 1) for v in ctl_depth.values()] or [0.0])
        transcribed = local > max(1.5, 3.0 * bg)
        locus["wt_junction_spanning_reads"] = st["wt_spanning"]
        locus["wt_flank_depth"] = [left_d, right_d]
        locus["wt_background_depth"] = bg
        if not transcribed:
            verdict = "undetermined_low_coverage"
        elif st["wt_spanning"] >= max(3, 0.5 * local):
            verdict = "wt_allele_transcribed"
        else:
            verdict = "no_wt_transcript_detected"
        expected = int(round(local * 150))
        locus["zygosity_call"] = {
            "verdict": verdict,
            "wt_spanning_reads": st["wt_spanning"],
            "expected_spanning_if_wt_present": expected,
            "wt_flank_depth": [left_d, right_d],
            "matched_control_depth": bg,
            "note": "suggestive only -- RNA-seq cannot prove zygosity, and absence of "
                    "wild-type transcripts may mean homozygosity or simply that the site "
                    "is not transcribed. Confirm on genomic DNA."}
    else:
        locus["zygosity_call"] = {"verdict": "not_assessed",
                                  "note": "genome not supplied to validate_locus()"}
    return locus


def annotate(locus, gff_path: str, max_genes: int = 5):
    if not gff_path or not os.path.exists(gff_path):
        return []
    chrom, pos = locus["chrom"], locus["u_left"]
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom or p[2] != "gene":
                continue
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            if e < pos - 200_000 or s > pos + 200_000:
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            genes.append({"id": attrs.get("Name") or attrs.get("ID", "?"),
                          "start": s, "end": e, "strand": p[6],
                          "distance": 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e)),
                          "inside": s <= pos <= e})
    genes.sort(key=lambda g: g["distance"])
    return genes[:max_genes]


# ----------------------------------------------------------------------------
# optional: junction-spanning primer pairs for Sanger validation
# ----------------------------------------------------------------------------
def _tm(seq: str) -> float:
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = len(s) - gc
    return 64.9 + 41 * (gc - 16.4) / len(s) if len(s) else 0.0


def _ok_primer(p: str) -> bool:
    if not (18 <= len(p) <= 26):
        return False
    gc = (p.count("G") + p.count("C")) / len(p)
    if not (0.35 <= gc <= 0.65):
        return False
    if p[-1] not in "GC":
        return False
    if re.search(r"(.)\1{3,}", p):
        return False
    if p.count("G") + p.count("C") == 0:
        return False
    return True


def _pick(seq: str, forward: bool, tmin: float = 55.0, tmax: float = 65.0):
    """pick the first acceptable primer in seq (5'->3' of the given strand)"""
    best = None
    for L in (20, 22, 24):
        for i in range(0, max(1, len(seq) - L)):
            p = seq[i:i + L] if forward else revcomp(seq[len(seq) - i - L:len(seq) - i])
            if not _ok_primer(p):
                continue
            t = _tm(p)
            if tmin <= t <= tmax:
                return {"seq": p, "len": L, "tm": round(t, 1),
                        "pos": i if forward else len(seq) - i - L}
            if best is None:
                best = {"seq": p, "len": L, "tm": round(t, 1),
                        "pos": i if forward else len(seq) - i - L}
    return best


def design_primers(locus):
    """4 primers: plant-F/T-DNA-R (LB junction) and T-DNA-F/plant-R (RB junction)"""
    contig = locus["contig"]
    lb = locus["lflank_len"]
    rb = locus["lflank_len"] + locus["tdna_len"]
    out = {}
    p = _pick(contig[max(0, lb - 500):lb - 150], True)
    if p:
        out["LB_plant_F"] = p
    p = _pick(contig[lb + 120:lb + 420], False)
    if p:
        out["LB_TDNA_R"] = p
    p = _pick(contig[rb - 420:rb - 120], True)
    if p:
        out["RB_TDNA_F"] = p
    p = _pick(contig[rb + 150:rb + 500], False)
    if p:
        out["RB_plant_R"] = p
    # expected product sizes
    if "LB_plant_F" in out and "RB_plant_R" in out:
        # both plant primers -> wild-type allele (T-DNA absent)
        f = locus["u_left"] - 500 + out["LB_plant_F"]["pos"]
        r = locus["u_right"] + 150 + out["RB_plant_R"]["pos"] + out["RB_plant_R"]["len"]
        out["WT_allele_product_bp"] = r - f
    if "LB_plant_F" in out and "LB_TDNA_R" in out:
        a = max(0, lb - 500) + out["LB_plant_F"]["pos"]
        b = lb + 120 + out["LB_TDNA_R"]["pos"] + out["LB_TDNA_R"]["len"]
        out["LB_expected_product_bp"] = b - a
    if "RB_TDNA_F" in out and "RB_plant_R" in out:
        a = rb - 420 + out["RB_TDNA_F"]["pos"]
        b = rb + 150 + out["RB_plant_R"]["pos"] + out["RB_plant_R"]["len"]
        out["RB_expected_product_bp"] = b - a
    return out


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True, help="R1 fastq(.gz) of the 3' mRNA-seq library")
    ap.add_argument("--plasmid", required=True, help="plasmid GenBank (.gb) or FASTA")
    ap.add_argument("--genome", required=True, help="genome FASTA (needs .fai -> samtools faidx)")
    ap.add_argument("--gff", default=None, help="optional GFF3 for annotation")
    ap.add_argument("--outdir", default="tdna_out")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-plasmid-part", type=int, default=30,
                    help="min bp of plasmid required in a read [30]")
    ap.add_argument("--min-plant-tail", type=int, default=25,
                    help="min bp of plant sequence required at the junction [25]")
    ap.add_argument("--border-window", type=int, default=100,
                    help="max distance (bp) of the T-DNA breakpoint from a border repeat [100]")
    ap.add_argument("--cluster-slack", type=int, default=250)
    ap.add_argument("--validate", action="store_true", default=True,
                    help="realign all reads to the locus to count junction support")
    ap.add_argument("--no-validate", dest="validate", action="store_false")
    ap.add_argument("--reuse-hits", action="store_true",
                    help="reuse <outdir>/plasmid_hits.pkl from a previous run")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()
    log = open(os.path.join(args.outdir, "run.log"), "w")

    def say(*a):
        msg = " ".join(str(x) for x in a)
        print(msg)
        log.write(msg + "\n")
        log.flush()

    say("=== 1. plasmid ===")
    plasmid = Plasmid(args.plasmid)
    say(f"  length {plasmid.L} bp; borders: "
        f"LB={plasmid.borders['LB']}, RB={plasmid.borders['RB']}")
    say(f"  T-DNA span (LB..RB): {plasmid.tdna[0]}..{plasmid.tdna[1]} "
        f"({plasmid.tdna[1]-plasmid.tdna[0]+1} bp)")
    plasmid.write_fasta(os.path.join(args.outdir, "plasmid.fasta"))
    json.dump(plasmid.feature_table(),
              open(os.path.join(args.outdir, "plasmid_features.json"), "w"), indent=1)

    say("=== 2. screen reads against the plasmid ===")
    cache = os.path.join(args.outdir, "plasmid_hits.pkl")
    if args.reuse_hits and os.path.exists(cache):
        import pickle
        total, hits = pickle.load(open(cache, "rb"))
        say(f"  reused cached screen: {total:,} reads, {len(hits):,} plasmid hits")
    else:
        total, hits = screen(args.r1, plasmid, args.threads, args.min_plasmid_part)
        import pickle
        pickle.dump((total, hits), open(cache, "wb"))

    say("=== 3. split reads into plasmid / plant moieties ===")
    cand = split_hits(hits, plasmid, args.border_window, args.min_plant_tail)
    if not cand:
        say("NO border-proximal chimeric reads found -> no insertion detected.")
        return

    say("=== 4. map plant moieties to the genome ===  (whole-genome index)")
    kept = map_to_genome(cand, args.genome, args.threads, args.min_plant_tail)
    if not kept:
        say("no genuinely chimeric junction survived -> no insertion detected.")
        return

    say("=== 5. cluster ===")
    sites = cluster_candidates(kept, args.cluster_slack)
    for i, s in enumerate(sites[:10]):
        say(f"  site {i+1}: {s['chrom']}:{s['start0']+1}-{s['end0']+1} "
            f"reads={s['n_reads']}")

    # drop singleton / background clusters: keep sites with >= max(3, 2% of the top) reads
    top = sites[0]["n_reads"] if sites else 0
    min_support = max(3, int(0.02 * top))
    sites = [x for x in sites if x["n_reads"] >= min_support]
    say(f"  keeping {len(sites)} site(s) with >= {min_support} supporting reads")

    say("=== 6. refine breakpoints ===")
    genome = pysam.FastaFile(args.genome)
    for s in sites[:5]:
        refine(s, plasmid, genome)
        r = s.get("refined", [])
        if r:
            b = r[0]
            say(f"  site {s['chrom']}:{s['start0']+1}: best {b['order']} "
                f"plant={b['plant_pos']} plasmid={b['plasmid_pos']} support={b['support']}")

    say("=== 7. build + validate loci ===")
    loci = []
    built = 0
    for s in sites[:5]:
        lc = build_locus(s, plasmid, genome)
        if lc is None:
            say(f"  WARNING: junction(s) found at {s['chrom']}:{s['start0']+1} "
                f"({s['n_reads']} reads) but no locus could be assembled")
            continue
        built += 1
        if args.validate:
            validate_locus(lc, args.r1, args.threads, genome)
        lc["nearest_genes"] = annotate(lc, args.gff)
        lc["primers"] = design_primers(lc)
        lc["site"] = {"chrom": s["chrom"], "n_reads": s["n_reads"]}
        loci.append(lc)
        say(f"  {lc['chrom']}:{lc['u_left']}|{lc['u_right']}  "
            f"T-DNA plasmid {lc['tdna_start']}..{lc['tdna_end']} ({lc['tdna_len']} bp)  "
            f"spanLB={lc.get('support_span_LB')} spanRB={lc.get('support_span_RB')}")

    say("=== 8. write reports ===")
    for i, lc in enumerate(loci):
        fn = os.path.join(args.outdir, f"locus_{i+1}_{lc['chrom']}_{lc['u_left']}.fasta")
        with open(fn, "w") as f:
            f.write(f">locus{i+1} {lc['chrom']}:{lc['u_left']}|{lc['u_right']} "
                    f"T-DNA={lc['tdna_start']}..{lc['tdna_end']} "
                    f"deleted_plant={len(lc['deleted_plant'])}bp\n")
            f.write(wrap(lc["contig"]) + "\n")
        say(f"  wrote {fn}")

    with open(os.path.join(args.outdir, "junctions.tsv"), "w") as f:
        f.write("read\tchrom\tplant_pos\tplant_strand\tplasmid_pos\tplasmid_strand\t"
                "border\tborder_dist\tplasmid_part_bp\ttail_bp\tfull_genome_mlen\n")
        for c in kept:
            f.write(f"{c['name']}\t{c['chrom']}\t{c['gpos0']+1}\t{c['tail_strand']}\t"
                    f"{c['plasmid_coord']}\t{c['plasmid_dir']}\t{c['border']}\t"
                    f"{c['border_dist']}\t{c['plasmid_part']}\t{len(c['tail'])}\t"
                    f"{c['full_genome_mlen']}\n")

    for i, lc in enumerate(loci):
        if not lc.get("primers"):
            continue
        fn = os.path.join(args.outdir, f"primers_locus{i+1}.tsv")
        with open(fn, "w") as f:
            f.write("name\tsequence\tlength\tTm\n")
            for k, v in lc["primers"].items():
                if isinstance(v, dict):
                    f.write(f"{k}\t{v['seq']}\t{v['len']}\t{v['tm']}\n")
            for key in ("LB_expected_product_bp", "RB_expected_product_bp",
                        "WT_allele_product_bp"):
                if key in lc["primers"]:
                    f.write(f"{key}\t{lc['primers'][key]}\t\t\n")
        say(f"  wrote {fn}")

    summary = {"reads_scanned": total, "plasmid_reads": len(hits),
               "candidate_junctions": len(cand), "chimeric_junctions": len(kept),
               "input": {"r1": os.path.abspath(args.r1),
                         "plasmid": os.path.abspath(args.plasmid),
                         "genome": os.path.abspath(args.genome)},
               # everything a downstream plotting / reporting tool needs to know
               # about the vector, so nothing has to be hard-coded per construct
               "plasmid": {"name": os.path.basename(args.plasmid),
                           "length": plasmid.L,
                           "borders": {k: (list(v) if v else None)
                                       for k, v in plasmid.borders.items()},
                           "tdna_span": list(plasmid.tdna),
                           "features": plasmid.feature_table()},
               "params": {"border_window": args.border_window,
                          "min_plasmid_part": args.min_plasmid_part,
                          "min_plant_tail": args.min_plant_tail},
               "loci": [{k: v for k, v in lc.items() if k != "contig"} for lc in loci]}
    json.dump(summary, open(os.path.join(args.outdir, "summary.json"), "w"),
              indent=1, default=str)

    with open(os.path.join(args.outdir, "report.md"), "w") as f:
        f.write("# T-DNA insertion site report\n\n")
        f.write(f"- input reads: `{args.r1}`\n")
        f.write(f"- reads scanned: {total:,}\n")
        f.write(f"- reads containing plasmid sequence: {len(hits):,} "
                f"({100*len(hits)/max(1,total):.4f}%)\n")
        f.write(f"- chimeric, border-proximal junctions: {len(kept)}\n\n")
        if not loci:
            if kept:
                f.write("Junction reads were found but no locus could be assembled — "
                        "see run.log.\n")
            else:
                f.write("**No insertion site detected.**\n")
        for i, lc in enumerate(loci):
            f.write(f"## Insertion {i+1}: {lc['chrom']}:{lc['u_left']:,}"
                    f"–{lc['u_right']:,}\n\n")
            f.write(f"- T-DNA: plasmid {lc['tdna_start']}..{lc['tdna_end']} "
                    f"({lc['tdna_len']:,} bp)"
                    + ("   *(other border inferred from the annotation)*"
                       if lc.get("one_sided") else "") + "\n")
            if lc.get("one_sided"):
                f.write(f"- **one-sided evidence**: the {lc['one_sided']} junction was not "
                        f"observed in this library, so the {lc['one_sided']} T-DNA end and any "
                        f"deletion on that side are inferred from the plasmid annotation. "
                        f"The insertion **position** is still known to the base from the "
                        f"observed border.\n")
            if lc.get("lb_microhomology_bp", 1) > 1:
                f.write(f"- LB microhomology: {lc['lb_microhomology_bp']} bp -- the junction "
                        f"sequence is unambiguous, but the breakpoint can be drawn anywhere in "
                        f"{lc['lb_pos_span'][0]:,}..{lc['lb_pos_span'][1]:,}; "
                        f"0..{lc['lb_microhomology_bp']-1} bp of plant DNA may have been lost "
                        f"(constrained by the RB junction: {lc['del_options']} bp); "
                        f"breakpoint reported at {lc['u_left']:,}\n")
            dellen = len(lc["deleted_plant"])
            extra = " (assumed 0 — not measured on the unobserved side)" \
                if lc.get("one_sided") else \
                (f" ({lc['deleted_plant']})" if lc["deleted_plant"] else "")
            f.write(f"- deleted plant bases: {dellen}{extra}\n")
            for side in ("LB", "RB"):
                v = lc.get(f"support_span_{side}")
                if v is None:
                    f.write(f"- reads spanning {side} junction: not applicable "
                            f"(junction not observed)\n")
                else:
                    f.write(f"- reads spanning {side} junction: {v}\n")

            if lc.get("one_sided") == "LB":
                f.write(f"- RB junction (observed): `{lc['junction_RB']}`\n")
                f.write(f"- LB junction (inferred from annotation, not observed): "
                        f"`{lc['junction_LB']}`\n")
            elif lc.get("one_sided") == "RB":
                f.write(f"- LB junction (observed): `{lc['junction_LB']}`\n")
                f.write(f"- RB junction (inferred from annotation, not observed): "
                        f"`{lc['junction_RB']}`\n")
            else:
                f.write(f"- LB junction: `{lc['junction_LB']}`\n")
                f.write(f"- RB junction: `{lc['junction_RB']}`\n")

            # ---- zygosity: explicitly state what this library cannot answer ----
            z = lc.get("zygosity_call") or {}
            f.write("\n### Zygosity (homozygous vs heterozygous)\n\n")
            f.write("**Not determinable from a 3'mRNA-seq library.** RNA-seq reads come from "
                    "transcripts, and both chromosomes are transcribed in the same cells, so "
                    "allele dosage cannot be read out directly.\n\n")
            if "wt_junction_spanning_reads" in lc:
                f.write(f"- reads spanning the intact wild-type insertion site: "
                        f"{lc['wt_junction_spanning_reads']}\n")
                f.write(f"- transcript depth flanking that site: "
                        f"{lc['wt_flank_depth'][0]} / {lc['wt_flank_depth'][1]} reads per bp\n")
                if "wt_background_depth" in lc:
                    f.write(f"- matched control windows ~50 kb away: "
                            f"{lc['wt_background_depth']} reads per bp (background)\n")
                exp = z.get("expected_spanning_if_wt_present")
                if exp:
                    f.write(f"- if a wild-type transcript crossed the site at that depth, "
                            f"~{exp} spanning reads would be expected; "
                            f"{lc['wt_junction_spanning_reads']} were observed\n")
            verdict = z.get("verdict")
            if verdict == "undetermined_low_coverage":
                f.write("- verdict: **undetermined** — the site is not transcribed above the "
                        "matched background, so the absence of wild-type junction reads is "
                        "uninformative (a homozygote and a heterozygote would look identical).\n")
            elif verdict == "no_wt_transcript_detected":
                f.write("- verdict: no transcript crossing the site was detected — "
                        "*suggestive* of homozygosity (no wild-type allele transcribed), but "
                        "not proof: it can also mean the wild-type allele is silent or "
                        "untranscribed there.\n")
            elif verdict == "wt_allele_transcribed":
                f.write("- verdict: wild-type transcripts across the site were detected, so at "
                        "least one chromosome carries the wild-type allele (heterozygote, or "
                        "chimeric/mixed tissue).\n")
            f.write("\nConfirm on genomic DNA with the three-primer assay below "
                    "(a homozygote gives no wild-type band).\n")

            if lc.get("primers"):
                pr = lc["primers"]
                f.write("\n### Genotyping primers\n\n")
                f.write("| allele | primer pair | expected product |\n|---|---|---|\n")
                pairs = [("wild type", "LB_plant_F", "RB_plant_R", "WT_allele_product_bp"),
                         ("T-DNA (LB junction)", "LB_plant_F", "LB_TDNA_R", "LB_expected_product_bp"),
                         ("T-DNA (RB junction)", "RB_TDNA_F", "RB_plant_R", "RB_expected_product_bp")]
                for lab, a, b, key in pairs:
                    if a in pr and b in pr:
                        size = f"{pr[key]} bp" if key in pr else "n/a"
                        f.write(f"| {lab} | {a} + {b} | {size} |\n")
                f.write("\nA homozygote gives only the T-DNA bands; a heterozygote gives a "
                        "T-DNA band and the wild-type band; a wild type gives only the "
                        "wild-type band.\n\n")
                f.write("- validation primers (junction-spanning):\n")
                for k, v in lc["primers"].items():
                    if isinstance(v, dict):
                        f.write(f"  - {k}: `{v['seq']}` (Tm {v['tm']} C)\n")
                for key, lab in (("LB_expected_product_bp", "LB-junction product (plant-F + T-DNA-R)"),
                                 ("RB_expected_product_bp", "RB-junction product (T-DNA-F + plant-R)"),
                                 ("WT_allele_product_bp", "wild-type allele (both plant primers)")):
                    if key in lc["primers"]:
                        f.write(f"  - expected {lab}: {lc['primers'][key]} bp\n")
            if lc["nearest_genes"]:
                f.write("- nearest genes:\n")
                for g in lc["nearest_genes"]:
                    f.write(f"  - {g['id']} ({g['start']:,}–{g['end']:,}, {g['strand']})"
                            f" distance {g['distance']:,} bp"
                            f"{' (contains insertion)' if g['inside'] else ''}\n")
            f.write("\n")
    say(f"done in {time.time()-t0:.0f}s -> {args.outdir}/report.md")
    log.close()


if __name__ == "__main__":
    main()
