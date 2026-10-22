# T-DNA insertion-site finder from 3′-mRNA-seq

`find_tdna_insertion.py` locates T‑DNA / transgene insertion sites in a plant genome
starting from a 3′‑mRNA‑seq (or any RNA‑seq / DNA‑seq) short‑read library plus the
transformation plasmid. It is written to be reused unchanged on other transformation
lines — only the reads, the plasmid and the genome change.

---

## 1. Principle

A T‑DNA inserted in the genome produces transcripts that run from plant DNA into
T‑DNA DNA (and vice‑versa). Such a transcript yields **chimeric reads** that are

* partially T‑DNA, partially plant genome,
* not mappable end‑to‑end to either genome or plasmid, and
* positioned at a **T‑DNA border** (LB or RB repeat) on the plasmid side.

The plant moiety of each chimeric read is mapped back to the genome; the mapping
position *is* the insertion site. Both borders of the same insertion land on the same
plant coordinate, which is the internal consistency check that separates real
insertions from artefacts.

Two artefact classes are explicitly filtered out, because both occur in practice:

| artefact | why it happens | filter |
|---|---|---|
| reads that are 100 % plant DNA but share a few dozen bp with the plasmid | the plasmid carries plant‑derived sequence (a native promoter, a miRNA precursor) | the **full read** is aligned to the genome; if it maps end‑to‑end it is discarded |
| poly(A) / adapter read‑through, template‑switch GGG | library chemistry | poly(A/T) trimming, adapter motifs, low‑complexity filter |
| T‑DNA coordinates far from any border | rearranged vector backbone, homology | only junctions within `--border-window` bp of LB/RB are kept |

---

## 2. Requirements

* `mappy` — minimap2 bindings (alignment, builds the genome index in memory)
* `pysam` — indexed genome access (`samtools faidx genome.fa` must have been run)
* `matplotlib` — only for `plot_insertion.py`
* ~4 GB RAM per Gb of genome; add `--threads` for speed

Tested with Python 3.15 · mappy 2.31 · pysam 0.24.1 · matplotlib 3.11.2.

Two ways to install, pick one:

```bash
# A) into the Homebrew Python already on your PATH  (what this machine uses)
/opt/homebrew/bin/python3 -m pip install --break-system-packages mappy pysam

# B) into an isolated virtual environment  (does not touch the system python)
python3 -m venv .venv && . .venv/bin/activate
pip install mappy pysam matplotlib
```

With option A everything runs as plain `python3 find_tdna_insertion.py ...`.
Note that a Homebrew Python upgrade (`brew upgrade python@3.15`) moves to a new
`site-packages` and the pip-installed modules have to be reinstalled; option B is
immune to that.

If matplotlib complains `~/.matplotlib is not a writable directory`, either
`mkdir -p ~/.matplotlib` or export `MPLCONFIGDIR` to a writable path
(`plot_insertion.py` already sets it to `.mplcache/` next to itself).

---

## 3. Usage

```bash
python3 find_tdna_insertion.py \
    --r1        C110M_UMI_R1.fastq.gz \
    --plasmid   pTF101-dsRED_clv3.gb \
    --genome    Gmax_880_v6.0.fa \
    --gff       Gmax_880_Wm82.a6.v1.gene_exons.gff3 \
    --outdir    out_C110M \
    --threads   8
```

Useful options

| option | default | meaning |
|---|---|---|
| `--min-plasmid-part` | 30 | minimum T‑DNA bases a read must contain to be considered |
| `--min-plant-tail` | 25 | minimum plant bases at the junction |
| `--border-window` | 100 | max distance of the breakpoint from the LB/RB repeat |
| `--cluster-slack` | 250 | bp within which junctions are merged into one locus |
| `--reuse-hits` | off | reuse `<outdir>/plasmid_hits.pkl` (skips the slow screen) |
| `--no-validate` | off | skip the second pass that counts junction‑spanning reads |

The plasmid may be given as **GenBank** (borders are taken from the `LB/RB T-DNA repeat`
features) or as **FASTA** (borders are then located by the conserved
`TGGCAGGATATATTGTGGTGTAAAC` border motif).

---

## 4. Output

```
outdir/
├── report.md                  human-readable summary (the answer)
├── summary.json               machine-readable
├── junctions.tsv              every chimeric junction read
├── locus_1_<chr>_<pos>.fasta  reconstructed transgene locus  (plant|T-DNA|plant)
├── primers_locus1.tsv         junction-spanning primers + expected product sizes
├── plasmid.fasta              plasmid + circularly doubled reference
├── plasmid_features.json      parsed feature table
├── plasmid_hits.pkl           cache of the screening step
└── run.log
```

`junctions.tsv` columns:
`read, chrom, plant_pos, plant_strand, plasmid_pos, plasmid_strand, border, border_dist,
plasmid_part_bp, tail_bp, full_genome_mlen`

`report.md` contains, per insertion: the plant coordinate, the integrated plasmid
coordinates, the number of reads spanning each junction, the junction sequences, any
microhomology at the breakpoint, the nearest genes, and validation primers.

---

## 4b. Plotting

`plot_insertion.py` draws the summary figure for **any** construct or sample. It is
driven entirely by the files above — plasmid length, borders, the GenBank feature
table, the T-DNA span, the plant coordinates and all read counts are read from
`summary.json` / `junctions.tsv` / `plasmid_hits.pkl`. Nothing about a particular
vector, gene or coordinate is hard-coded.

```bash
python3 plot_insertion.py --outdir out_C1R
python3 plot_insertion.py --outdir out_myLine --title "line myLine T2"
python3 plot_insertion.py --outdir out_myLine --locus 2        # second insertion
python3 plot_insertion.py --outdir out_myLine --sample-reads 0 # skip panel A
```

| option | default | meaning |
|---|---|---|
| `--locus` | 1 | which insertion to draw |
| `--title` | outdir name | figure title / sample label |
| `--sample-reads` | 60 | junction reads drawn per border; `0` drops panel A |
| `--max-labels` | 14 | max features labelled on the vector map |
| `--out`, `--dpi` | `<outdir>/insertion_summary.png`, 160 | output file |

Panel A is skipped automatically when `plasmid_hits.pkl` is absent (e.g. when the
screen was run elsewhere), and panel B degrades to a bare T-DNA bar when the plasmid
was supplied as FASTA and therefore has no feature table.

---

## 5. Interpretation notes

* **Microhomology.** T‑DNA integration frequently deletes a few plant bases and leaves
  a short stretch where plant and T‑DNA sequence coincide. The junction *sequence* is
  then unambiguous, but the coordinate can be drawn anywhere inside that stretch. The
  script reports the whole equivalence interval (`lb_pos_span`) and the most
  parsimonious breakpoint.
* **Border trimming.** The integrated T‑DNA usually does not start/end exactly at the
  annotated repeats; the report gives the actually integrated plasmid interval, so the
  number of bases lost from each border can be read off directly.
* **One‑sided evidence.** An RNA‑seq library often shows only one border (the read‑through
  transcript may only exist on one side). The script then reports `one_sided`, states that
  the missing border is inferred from the plasmid annotation, and marks that side's
  "reads spanning … junction" as *not applicable* rather than inventing a number. The
  insertion **position** is still base‑resolved from the observed border.
* **Multiple insertions.** Every independent locus appears as its own `## Insertion N`
  section. Tandem or inverted T‑DNA arrays show up as several junctions at the same
  plant position; unintegrated plasmid would show a plasmid‑only circle junction and no
  plant junction.
* **Transcript orientation.** On a 3′‑mRNA‑seq library, junction reads are only produced
  where the transgene or a neighbouring plant gene is transcribed across the junction.
  One border therefore often has far fewer reads than the other — that is expected and
  does not mean the border is absent. If a line is completely silent, use DNA‑seq
  (whole‑genome or capture) with the same script; the algorithm is library‑agnostic.

---

## 5b. Zygosity is *not* determined by this tool

Zygosity (homozygous vs heterozygous) is a property of the **genome**; 3′‑mRNA‑seq
sequences **transcripts**. Both chromosomes sit in the same cells and both are
transcribed, so allele dosage cannot be read out of junction read counts. A further
practical problem: the junction reads here come from the transgene's own strong
promoters (2x35S, Ubi10), not from the plant locus, so their number reflects
expression, not copy number.

The report states this explicitly and, as a diagnostic, counts reads that span the
**intact wild‑type** insertion site (a reconstructed contig with the T‑DNA removed):

| finding | report verdict |
|---|---|
| flanking depth ≤ 3× a length‑matched control ~50 kb away | `undetermined_low_coverage` — the site is untranscribed; the absence of WT reads means nothing |
| flanking depth adequate, WT‑spanning reads ≥ threshold | `wt_allele_transcribed` — a wild‑type allele is present (heterozygote / chimera) |
| flanking depth adequate, WT‑spanning reads ≈ 0 | `no_wt_transcript_detected` — consistent with homozygosity, still not proof |

**To actually call zygosity, use genomic DNA.** The script already designs the
three‑primer assay and prints the expected sizes (`<outdir>/primers_locus1.tsv` and
the *Genotyping primers* table in `report.md`):

```
wild type     LB_plant_F + RB_plant_R   817 bp    present only if a WT chromosome exists
T-DNA (LB)    LB_plant_F + LB_TDNA_R    636 bp
T-DNA (RB)    RB_TDNA_F  + RB_plant_R   829 bp
```

homozygote → 636/829 only; heterozygote → 636/829 **and** 817; wild type → 817 only.

Measured on this line (both libraries agree):

| | C110M (9.8 M reads) | C1R (170 M reads) |
|---|---|---|
| reads spanning the intact WT site | 1 | 6 |
| expected if a WT transcript crossed it | ~315 | ~8,800 |
| flank depth / matched control | 2.1 : 0.49 | 58.6 : 4.9 |
| verdict | `no_wt_transcript_detected` | `no_wt_transcript_detected` |

So **no transcript runs through the undisrupted site**, which is *suggestive* of
homozygosity — but note the caveat that makes this weak evidence here: the local
transcription is a 3′ end that stops ~360 bp upstream of the insertion
(Gm04:2,707,150), so even a heterozygous plant might produce no read crossing the
site. Only a DNA assay settles it.
Alternatives: ddPCR/qPCR of T‑DNA against a single‑copy reference (copy number 1 vs 2),
WGS/target‑capture read depth at the locus, or T2/T3 segregation of the selectable
marker (3:1 vs 1:0).

A `copy_number_estimate` field was removed from earlier versions of `summary.json`:
it divided T‑DNA expression depth by flank depth and produced a meaningless number
(8.94 for C1R) because the transgene is driven by strong promoters.

---

## 6. Result obtained for C1R / C110M

| | |
|---|---|
| insertion site | **Gm04:2,707,513 / 2,707,514** (Wm82.a6, `Gmax_880_v6.0`) |
| integrated T‑DNA | plasmid **6,433–14,573** (8,141 bp) |
| borders lost | 36 bp at the LB, 22 bp at the RB |
| plant DNA lost | 0–2 bp (5‑bp microhomology `CCTTA` around the LB breakpoint; junction sequence unambiguous) |
| support | 638 chimeric junction reads at the locus (610 at 2,707,514, 27 at 2,707,511, 1 at 2,707,516) + 1 unlinked singleton; after re‑alignment 1,710 reads span the RB junction and 77 the LB junction |
| copies | single, no tandem/inverted array |
| context | intergenic — 3.2 kb downstream of `Glyma.04G033900`, 2.8 kb upstream of `Glyma.04G034000` |

---

## 7. Reproducibility

Both libraries give the identical answer:

| library | reads | plasmid reads | chimeric junctions | LB span | RB span |
|---|---|---|---|---|---|
| `C1R_R1_001.fastq.gz` (raw) | 170,091,432 | 349,926 | 638 (+1 singleton) | 77 | 1,710 |
| `C110M_UMI_R1.fastq.gz` (10 M subset) | 9,764,529 | 20,042 | 35 | 8 | 93 |

The BAM of the 10 M subset independently contains soft‑clipped reads whose clipped
part is T‑DNA at the same coordinate, e.g. `Gm04:2,707,514` with a 43 bp clip matching
plasmid 14,531–14,573.

---

## 8. Result obtained for DA1

`DA1_R1_001.fastq.gz` — 30,694,140 reads, 126,356 (0.41 %) mapping to the plasmid.

| | |
|---|---|
| insertion site | **Gm04:2,707,513 / 2,707,514** — identical to C1R / C110M |
| observed junction | RB: plasmid 14,573 \| Gm04:2,707,514 (12 junction reads, 194 reads span it on re-alignment) |
| LB junction | **not observed** in this library → T-DNA end inferred from the annotation (6,397–14,573), flagged `one_sided` in the report |
| copies | single |
| zygosity | `undetermined` — the site is not transcribed above the matched control |

The identical coordinate in two independent sequencing runs is either confirmation
(same transformation event, sibling or re-sequenced plant) or a sample-identity
problem if DA1 was meant to be an independent transformant. Library overlap with C1R:
54.5 % of DA1's distinct 40‑mers occur in C1R (68.9 % of DA1 reads), yet the T‑DNA read
fraction differs 2× (0.41 % vs 0.21 %) and the junction yield differs ~50× (12 vs 639),
so the two libraries are not simply the same pool sequenced to different depths.

Two tool bugs were found and fixed through DA1:

1. it printed *"No insertion site detected"* even though a strong RB junction existed —
   `build_locus` demanded both borders. One‑sided loci are now assembled and reported.
2. it displayed the annotation‑inferred border as if it had been observed, and printed a
   "reads spanning" count for that assumed position. Both are now marked as inferred /
   not applicable.
