#!/usr/bin/env python3
"""
plot_insertion.py
=================
Draw the summary figure for a T-DNA insertion site, for ANY construct / sample,
using only the files that ``find_tdna_insertion.py`` writes:

    <outdir>/summary.json          plasmid length, borders, features, loci, counts
    <outdir>/junctions.tsv         every chimeric junction read
    <outdir>/plasmid_hits.pkl      screening cache -> plasmid read coverage

Nothing about a particular plasmid, gene, coordinate or read number is built in.

Panels
  A  plasmid read depth of the library, integrated region shaded, LB/RB marked
  B  vector map from the GenBank feature table (omitted for a FASTA-only plasmid)
  C  the chimeric junction reads at each border -- plant moiety vs T-DNA moiety

Usage
  python3 plot_insertion.py --outdir out_C1R
  python3 plot_insertion.py --outdir out_C1R --locus 1 --dpi 200 --title "line C1R"
  python3 plot_insertion.py --outdir out_miR156a --sample-reads 0     # skip panel A

Requires only matplotlib (+ pickle/json from the stdlib).
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys

os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), ".mplcache"))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except ImportError:  # pragma: no cover
    sys.exit("matplotlib is required:  pip install matplotlib")

GOLD, BLUE, RED, GREY = "#E8B33A", "#4C72B0", "#C0392B", "#9a9a9a"


# ----------------------------------------------------------------------------
# input
# ----------------------------------------------------------------------------
def load(outdir: str):
    fj = os.path.join(outdir, "summary.json")
    if not os.path.exists(fj):
        sys.exit(f"error: {fj} not found -- run find_tdna_insertion.py first")
    summary = json.load(open(fj))
    if not summary.get("loci"):
        sys.exit("error: summary.json contains no insertion locus")

    junc = []
    ft = os.path.join(outdir, "junctions.tsv")
    if os.path.exists(ft):
        with open(ft) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                junc.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    cov = None
    fc = os.path.join(outdir, "plasmid_hits.pkl")
    if os.path.exists(fc):
        try:
            _, hits = pickle.load(open(fc, "rb"))
            L = summary["plasmid"]["length"]
            arr = [0] * (L + 1)
            for h in hits:
                for p in range(h["r_st"], h["r_en"]):
                    arr[(p % L) + 1] += 1
            cov = arr[1:]
        except Exception as exc:                      # cache from an older version
            print(f"  (plasmid coverage skipped: {exc})")
    return summary, junc, cov


def pick_features(features, max_labels=14):
    """largest annotated features, borders always kept"""
    def score(f):
        lab = (f.get("label") or "").lower()
        border = "left border" in lab or "right border" in lab
        return (border, f["end"] - f["start"])
    feats = [f for f in features if (f.get("label") or "").strip()]
    return sorted(feats, key=score, reverse=True)[:max_labels]


# ----------------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--locus", type=int, default=1, help="1-based insertion index [1]")
    ap.add_argument("--title", default=None, help="main title [sample = outdir name]")
    ap.add_argument("--out", default=None, help="output png [<outdir>/insertion_summary.png]")
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--sample-reads", type=int, default=60,
                    help="junction reads drawn per border [60]; 0 disables panel A")
    ap.add_argument("--max-labels", type=int, default=14,
                    help="max feature labels on the vector map [14]")
    args = ap.parse_args()

    summary, junc, cov = load(args.outdir)
    pl_info = summary["plasmid"]
    L = pl_info["length"]
    borders = pl_info.get("borders") or {}
    feats_all = pl_info.get("features") or []

    loci = summary["loci"]
    if not (1 <= args.locus <= len(loci)):
        sys.exit(f"error: --locus {args.locus} out of range (1..{len(loci)})")
    lc = loci[args.locus - 1]

    chrom = lc["chrom"]
    u_left, u_right = lc["u_left"], lc["u_right"]
    t_start, t_end = lc["tdna_start"], lc["tdna_end"]
    t_len = lc["tdna_len"]
    n_plasmid_reads = summary.get("plasmid_reads", 0)
    sample = args.title or os.path.basename(os.path.normpath(args.outdir))

    # junction reads belonging to this insertion
    def mine(r):
        try:
            pos = int(r["plant_pos"])
        except (KeyError, ValueError):
            return False
        return r.get("chrom") == chrom and (u_left - 250) <= pos <= (u_right + 250)
    junc = [r for r in junc if mine(r)]
    lb = [r for r in junc if r.get("border") == "LB"] or \
         [r for r in junc if int(r["plant_pos"]) <= u_left]
    rb = [r for r in junc if r not in lb]

    has_cov = cov is not None and args.sample_reads != 0
    ratios = [1.35, 0.42, 1.5]
    fig = plt.figure(figsize=(15.5, 10.5 if has_cov else 7.6))
    gs = fig.add_gridspec(3 if has_cov else 2, 2,
                          height_ratios=ratios if has_cov else ratios[1:],
                          hspace=0.5 if has_cov else 0.32, wspace=0.14)
    r0 = 0

    # ---------------------------------------------------------- panel A: plasmid
    if has_cov:
        ax = fig.add_subplot(gs[r0, :]); r0 += 1
        step = max(10, L // 600)
        xs, ys = [], []
        for i in range(1, L, step):
            w = cov[i:i + step]
            xs.append(i + step // 2)
            ys.append(sum(w) / len(w))
        ax.fill_between(xs, ys, color=BLUE, lw=0)
        ax.axvspan(t_start, t_end, color=GOLD, alpha=0.22, zorder=0)
        for name, col in (("LB", RED), ("RB", RED)):
            b = borders.get(name)
            if b:
                ax.axvline(b[0] if name == "LB" else b[1], color=col, ls="--", lw=1.1)
        ax.set_xlim(0, L)
        ax.set_ylim(0, max(ys) * 1.30 if max(ys) else 1)
        ax.set_ylabel("read depth\n(reads / bp)", fontsize=9)
        ax.set_xlabel(f"{pl_info.get('name','plasmid')} position (bp)", fontsize=9)
        ax.set_title(f"A   T-DNA expression — {n_plasmid_reads:,} library reads map to the "
                     f"plasmid; integrated region (plasmid {t_start:,}–{t_end:,}) shaded",
                     fontsize=10.5, loc="left")
        ax.tick_params(labelsize=8)
        for name in ("LB", "RB"):
            b = borders.get(name)
            if b:
                ax.annotate(name, xy=((b[0] + b[1]) / 2, ax.get_ylim()[1] * 0.97),
                            ha="center", fontsize=8, color=RED, fontweight="bold")
        # a few informative labels inside the shaded region
        inner = [f for f in pick_features(feats_all, args.max_labels)
                 if t_start <= f["start"] <= t_end]
        inner.sort(key=lambda f: f["start"])
        for k, f in enumerate(inner):
            h = ax.get_ylim()[1] * (0.42 + 0.14 * (k % 3))
            ax.annotate(f["label"][:26], xy=((f["start"] + f["end"]) / 2, 0),
                        xytext=((f["start"] + f["end"]) / 2, h), ha="center",
                        fontsize=7, color="#444444", rotation=90, va="center", alpha=0.85)

    # ---------------------------------------------------------- panel B: vector map
    axb = fig.add_subplot(gs[r0, :]); r0 += 1
    axb.set_xlim(0, L); axb.set_ylim(0, 1); axb.axis("off")
    axb.add_patch(Rectangle((1, .30), L, .28, fc="#e6e6e6", ec="#999999", lw=.6))
    axb.add_patch(Rectangle((t_start, .30), max(t_end - t_start, 1), .28,
                            fc=GOLD, ec="none", alpha=.75))
    axb.text((t_start + t_end) / 2, .10,
             f"T-DNA transferred: plasmid {t_start:,}–{t_end:,}   ({t_len:,} bp)",
             ha="center", fontsize=8.5, color="#8a6d1f", fontweight="bold")
    feats = pick_features(feats_all, args.max_labels)
    if feats:
        for f in feats:
            col = RED if (f.get("label") or "").lower().find("border") >= 0 else "#3b6ea5"
            axb.add_patch(Rectangle((f["start"], .34), max(f["end"] - f["start"], L * 0.001),
                                    .20, fc=col, ec="none", alpha=.6))
        feats.sort(key=lambda f: f["start"])
        for k, f in enumerate(feats):
            up = (k % 2 == 0)
            y = .68 + (0.16 if up else 0.0)
            axb.plot([(f["start"] + f["end"]) / 2] * 2, [.54, y], color="#bbbbbb", lw=.5)
            axb.annotate(f["label"][:30], xy=((f["start"] + f["end"]) / 2, y + .05),
                         ha="center", fontsize=7.5, color="#222222")
        axb.set_title("B   vector map (annotation of the supplied plasmid)", fontsize=10.5,
                      loc="left")
    else:
        axb.set_title("B   vector map (no feature table supplied — FASTA plasmid)",
                      fontsize=10.5, loc="left")

    # ---------------------------------------------------------- panel C: junctions
    span = 220
    for col, (name, recs) in enumerate([("LB", lb), ("RB", rb)]):
        axc = fig.add_subplot(gs[r0, col])
        if not recs:
            axc.axis("off")
            axc.set_title(f"C   {name} junction — no reads", fontsize=10, loc="left")
            continue
        recs = sorted(recs, key=lambda r: int(r.get("tail_bp", 0)))
        keep = max(1, len(recs) // args.sample_reads) if args.sample_reads else 1
        show = recs[::keep][:args.sample_reads] if args.sample_reads else []
        for i, r in enumerate(show):
            try:
                tail = int(r["tail_bp"]); plpart = int(r["plasmid_part_bp"])
            except (KeyError, ValueError):
                continue
            if name == "LB":
                axc.plot([-tail, 0], [i, i], color=BLUE, lw=1.1, solid_capstyle="butt")
                axc.plot([0, plpart], [i, i], color=GOLD, lw=1.1, solid_capstyle="butt")
            else:
                axc.plot([-plpart, 0], [i, i], color=GOLD, lw=1.1, solid_capstyle="butt")
                axc.plot([0, tail], [i, i], color=BLUE, lw=1.1, solid_capstyle="butt")
        axc.axvline(0, color=RED, ls="--", lw=1.2)
        axc.set_xlim(-span, span)
        axc.set_ylim(-2, len(show) + 3)
        axc.set_yticks([])
        axc.set_xlabel("distance from the junction (bp)", fontsize=9)
        axc.tick_params(labelsize=8)
        for sp in ("top", "right", "left"):
            axc.spines[sp].set_visible(False)
        if name == "LB":
            head = (f"C   LB junction — {chrom}:{u_left:,} | plasmid {t_start:,}")
            axc.text(-span + 10, len(show) + 1.2, "plant", fontsize=8, color=BLUE)
            axc.text(10, len(show) + 1.2, "T-DNA", fontsize=8, color="#8a6d1f")
        else:
            head = (f"RB junction — plasmid {t_end:,} | {chrom}:{u_right:,}")
            axc.text(-span + 10, len(show) + 1.2, "T-DNA", fontsize=8, color="#8a6d1f")
            axc.text(10, len(show) + 1.2, "plant", fontsize=8, color=BLUE)
        if args.sample_reads == 0:
            head += "   (panel C disabled)"
        axc.set_title(head + f"\n({len(show)} junction reads sampled of {len(recs)})",
                      fontsize=10, loc="left")

    n_del = len(lc.get("deleted_plant") or "")
    n_del = n_del if n_del else (u_right - u_left - 1 if u_right > u_left else 0)
    fig.suptitle(f"T-DNA insertion site — {sample}:  {chrom}:{u_left:,} / {u_right:,}"
                 f"   ({n_del} bp plant deleted)",
                 fontsize=12.5 if len(sample) < 40 else 11, y=0.975)

    out = args.out or os.path.join(args.outdir, "insertion_summary.png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  locus {args.locus}/{len(loci)}: {chrom}:{u_left:,}/{u_right:,}  "
          f"T-DNA {t_start:,}-{t_end:,} ({t_len:,} bp)  "
          f"junction reads: {len(lb)} LB / {len(rb)} RB")


if __name__ == "__main__":
    main()
