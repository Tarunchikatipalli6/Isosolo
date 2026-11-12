# isosolo

Design PCR primer pairs that amplify **one transcript isoform** and not its siblings.
Most primer tools solve the opposite problem. Primer3 has no concept of isoforms.
Primer-BLAST screens genome-wide but will happily return a pair hitting every
transcript of your gene. 

isosolo has an exclusion objective: amplify the target, exclude the siblings.

Two scripts, no packaging:
 
| File | What it does |
|---|---|
| `isosolo.py` | designs primer pairs |
| `check_primers.py` | checks a pair you already have, by brute force |

---
 
## Install
 
```bash
pip install primer3-py requests
```
 
Sequences come from the Ensembl REST API, so nothing needs downloading. The first
run of a gene needs internet; after that its sequences are cached in
`~/.ensembl_cache` and reruns are instant.

### Optional, but recommended
 
Without these two files the sibling screen still runs, but the **repeat filter** and
the **off-gene check** do not. The repeat filter is what catches Alu elements —
unique within a gene, present a million times genome-wide.
 
```bash
mkdir -p ~/refs && cd ~/refs
 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz
gunzip gencode.v50.transcripts.fa.gz
 
# GENCODE headers are 74 characters; -parse_seqids caps IDs at 50.
sed 's/^\(>[^|]*\).*/\1/' gencode.v50.transcripts.fa \
    > gencode.v50.transcripts.clean.fa
rm gencode.v50.transcripts.fa
 
makeblastdb -in gencode.v50.transcripts.clean.fa \
    -dbtype nucl -parse_seqids -out gencode_v50_tx
```
 
About 3 GB. Then set two paths near the top of `isosolo.py`:
 
```python
TRANSCRIPTOME_FASTA = "/path/to/refs/gencode.v50.transcripts.clean.fa"
BLAST_DB            = "/path/to/refs/gencode_v50_tx"
```
 
---

## Run
 
Edit one block:
 
```python
TARGET_JOBS = [
    {"gene": "MET", "transcript": "ENST00000436117"},
]
 
OUTPUT_ROOT = "/path/to/isosolo_output"
```
 
```bash
python isosolo.py
```
 
Or pass one target on the command line, which gives it its own output folder:
 
```bash
python isosolo.py MET ENST00000436117
```
 
Runtime is 1–3 minutes without the reference files, about 8 minutes with them. The
repeat filter dominates — it streams the whole 1.5 GB transcriptome once per target.
 
> **`OUTPUT_ROOT` is wiped at the start of every run.** Running targets one at a time
> with the same path deletes the previous result. Use the command-line form, which
> writes each target to its own folder.
 
---

## Reading the log
 
```
MET | ENST00000436117
target MET-204 (ENST00000436117), 6722 nt
19 exon-exon junction(s)
11 sibling isoform(s)
3 region(s) -> 1 interval(s), 39 bp unique
discriminating anchors: 19 of 39 interval positions have a 3' window
                        that differs from every sibling
12801 thermodynamically valid candidate(s) in 15.7s
  repeat filter 248s; 0/40 primer pairs sit in a repeat
  BLAST 5s; 0/40 hit another gene, 0/40 form a product with one primer alone
10/10 reported pairs amplify nothing else
```
 
Two lines matter:
 
- **`N bp unique`** — zero means no assay is possible.
- **`discriminating anchors: N of M`** — how many of those unique positions can
  actually hold a primer's 3′ end. Usually far fewer than M. See below.
---


## Output
 
One folder per target, named after the isoform.
 
| File | Contents |
|---|---|
| `*_assay.csv` | ranked primer pairs, ready to order |
| `*_assay.html` | two sequence tracks — unique regions above, primers below |
| `*_regions.html` | full cDNA with unique regions colour-coded |
| `*_regions.tsv` / `.txt` | each region with 50 bp of flanking context |
| `*_positions.csv` | per-base audit: every position against every sibling |
| `summary.csv` | one row per target, when more than one ran |
 
A run that fails at the uniqueness stage writes **nothing** — the folder is created
but stays empty. The reason is in the log.
 
---

## Choosing a pair
 
Take rank 1. Four columns decide it:
 
| Column | Want |
|---|---|
| `Offtarget_Isoforms_Amplified` | 0 |
| `Repeat_Derived` | no |
| `BLAST_Offgene_Amplicons` | 0 |
| `BLAST_Self_Amplicons` | 0 |
 
Then three that shape how you run it:
 
- **`Anchor_Class`** — prefer `chimeric_region`. A `base_difference` anchor
  discriminates on one base, which is allele-specific-PCR territory.
- **`Spans_Exon_Junction`** — `yes` means genomic DNA can't amplify. If `no`, use
  DNase and a no-RT control.
- **`Min_3p_Mismatches_vs_Any_Isoform`** — higher is safer, max 5. **Treat a 0 with
  suspicion rather than as a verdict**: the seeded search sometimes reports a phantom
  site where brute force finds no binding position at all. Run `check_primers.py`
  before dismissing a pair that otherwise looks good.
`Self_3p_Compl_F` / `_R` are Primer3's `SELF_END`, the same number Primer-BLAST
reports as "Self 3′ complementarity". Reported, not filtered. A 6 means the primer's
3′ end can prime off itself — worth a no-template control.
 
---

## When no pair exists
 
Most targets in a dense annotation have none, and the reason decides whether anything
can be done.
 
| Status | Meaning | Fixable? |
|---|---|---|
| `strict subset` | the whole transcript sits inside one named sibling | **No** — not by any method |
| `piecewise cover` | every stretch is in *some* sibling, none in all of them | Maybe, with a combination assay |
| `no_candidates` | unique sequence exists but no position can anchor a 3′ end, or nothing met the Tm/GC/amplicon limits | Sometimes — widen the constraints |
| `no_sequence` | Ensembl returned 500s or timed out | Yes — rerun |
 
The sibling named in a `strict subset` message is often **not** one you were
interested in. On ARPC5L the blocker was ARPC5L-207, a model that may not be
expressed in the tissue at all. That is the cost of designing against the full
annotation.
 
**A pattern worth knowing:** in every two-isoform pair tested so far, the short or
non-canonical isoform works and the long canonical one fails — TFRC-214 over
TFRC-201, MET-204 over MET-202, LRRTM4-206 over LRRTM4-201. A canonical transcript
contains all the exons and all the standard junctions, so the short variants are
built from pieces of it, but it has nothing they lack.
 
---

## Checking a pair
 
`check_primers.py` is a second, independent implementation. Where isosolo uses
k-mer seeding for speed, this checks **every position by brute force** — when the
two disagree, this one is right.
 
```python
GENE    = "MET"
FORWARD = "TTCAGTTACCGTGAAGATCCC"
REVERSE = "GCGATGTTGACATGCCTAATAA"
```
 
```bash
python check_primers.py
```
 
```
isoform   length   result               forward primer            reverse primer
MET-204     6722   AMPLIFIES  74 bp     pos 2601, 0 mm (0 at 3')  pos 2653, 0 mm (0 at 3')
MET-203     3389   no                   pos 2746, 0 mm (0 at 3')  pos 822, 8 mm (1 at 3')
...
SPECIFIC: only MET-204 amplifies (74 bp).
```
 
`no site` is the strongest possible result — no plausible binding position at all.
 
`No transcript amplifies` usually means a primer was typed on the wrong strand. Both
go in 5′→3′ exactly as ordered; the reverse is complemented by the script.
 
One gene only. Genes with many isoforms are slow — TFRC has 61 and takes a few
minutes. On a cluster, run it on a compute node.
 
---

## How it works
 
### 1. Fetch
 
The target plus every sibling isoform, from Ensembl REST.
 
### 2. Find unique sequence
 
Two detectors, because isoforms differ in two ways.
 
**`chimeric_region`** — slide a 41 bp window; if that exact string appears in no
sibling, the position is unique. This catches junctions: at a junction unique to your
target, the window's left half exists in the siblings and its right half exists in the
siblings, but never adjacent. The window is a chimera of two shared pieces.
 
**`base_difference`** — for single-base differences, 41 bp is too coarse. This locates
the equivalent region in each sibling with a 15 bp window at ≥90% identity, then
compares that one base. Catches shifted splice sites and short indels.
 
### 3. Decide where a primer may END
 
**This is the step that distinguishes isosolo from a generic designer.**
 
A primer's 3′ end is the only part that discriminates — the polymerase grabs that end
and builds outward, so a mismatch there stops extension while a mismatch near the 5′
end is almost free. (This is why 5′ tails and adapters amplify normally.)
 
Sitting inside a unique region is not enough. A base is flagged `chimeric_region` when
the 41 bp window **around** it is absent from the siblings, which is true for positions
up to 20 nt from the actual difference. A primer ending there carries the
discriminating base mid-primer, where the polymerase ignores it.
 
So a position may be a 3′ anchor only if, for every sibling, at least one of the
primer's last 5 bases differs from what that sibling carries there
(`REQUIRE_DISCRIMINATING_ANCHOR`).
 
On MET-204 this cut 39 "unique" positions to 19, and turned a pair that amplified 8 of
11 isoforms into one that amplified none.
 
### 4. Build pairs
 
Every length 18–24 nt at every allowed anchor, with every partner giving a 70–150 bp
product. Filtered on Tm 56–63 °C, GC 35–65%, Tm difference ≤5 °C, and hairpin/dimer
ΔG above −5 kcal/mol.
 
Both orientations are tried. Note that a reverse primer binds the sense strand
backwards, so its 3′-terminal base is the **leftmost** base of its footprint on the
cDNA.
 
### 5. Four specificity checks
 
Each catches what the others cannot.
 
| Check | Scope |
|---|---|
| Sibling in-silico PCR | other isoforms of this gene |
| Repeat filter | 3′ k-mer frequency across the transcriptome — catches Alu |
| BLAST forward × reverse | any other human transcript |
| BLAST forward × forward, reverse × reverse | one primer amplifying alone |
 
The last one matters: EMBOSS PrimerSearch, which BLAST replaced here, only ever pairs
primer A with primer B — a primer binding one template twice in opposite orientations
is structurally invisible to it. That is how an Alu primer passed every within-gene
check on a first TFRC design.
 
### 6. Rank
 
Repeats last, then off-gene hits, then sibling hits, then 3′ margin, then anchor class,
then junction-spanning, then thermodynamics. Specificity dominates chemistry.
 
---
 
## Parameters
 
All at the top of `isosolo.py`.
 
| Setting | Default | Source |
|---|---|---|
| `TERM_WINDOW` | 5 | Stadhouders 2010; same window Primer-BLAST uses |
| `TERM_MM_TO_BLOCK` | 2 | Kwok 1990 — one terminal mismatch is often tolerated, two block |
| `PB_MAX_TOTAL_MM` | 5 | Primer-BLAST ignores a site at 6+; Christopherson 1997 |
| `PRIMER_BLAST_COMPATIBLE` | True | use the two rules above, so surviving pairs also pass Primer-BLAST |
| `REQUIRE_DISCRIMINATING_ANCHOR` | True | the fix for window smearing (step 3) |
| `CHIMERIC_WINDOW` | 41 bp | inherited |
| `ALIGN_WINDOW` / `ALIGN_IDENTITY` | 15 bp / 0.90 | inherited; tolerates exactly one mismatch |
| `MIN` / `MAX_AMPLICON` | 70 / 150 bp | Bio-Rad ddPCR guidance is 75–200 |
| `PRIMER_TM_MIN` / `MAX` | 56 / 63 °C | Primer-BLAST defaults |
| GC range | 35–65% | looser than Bio-Rad's 40–60%, because the anchor position is fixed |
| `REPEAT_MAX_TRANSCRIPTS` | 50 | ours — Alu scores in the hundreds, a real primer 1–20 |
| `MAX_OFFTARGET_PRODUCT` | 1200 bp | ours |
| `GC_CLAMP_PARTNER` | True | G or C at the 3′ end of the **partner** primer only |
| `GC_CLAMP_BOTH` | False | tested and rejected — dropped the margin from 3/5 to 0/5 |
| `MAX_SELF_3P` | None | report Primer3 `SELF_END` without filtering, as Primer-BLAST does |
| `TOP_N` | 10 | raise to 50 to surface genuinely different anchors |
| `RUN_PRIMERSEARCH` | False | redundant once `BLAST_DB` is set |
 
The GC clamp can only apply to the partner primer. The anchored primer's 3′ base is
whichever base discriminates from the siblings — not ours to choose.
 
---
 
## Before ordering
 
1. Open `*_assay.html` and confirm the anchored primer sits **on** the highlighted
   unique region, not next to it. Five seconds, and it catches the worst failure mode.
2. Check the four columns above are zero.
3. Run `check_primers.py` on the chosen pair.
4. Run NCBI [Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast) —
   *Homo sapiens*, `refseq_rna`. It queries live NCBI rather than your local GENCODE
   snapshot.
5. At the bench: temperature gradient 58–66 °C, then a gel. No-template control for
   primer-dimer. DNase and no-RT if the amplicon doesn't span a junction.
### Reading a Primer-BLAST hit
 
Three questions settle most of them: are **both** primers present, is the product
short enough to amplify, and is the 3′ end clean? Three noes and it isn't a real
off-target — a 2 kb reverse-with-reverse product with a mismatched terminal base is
noise.
 
Note also that isosolo screens against **GENCODE** and Primer-BLAST against
**RefSeq**. They won't always agree even with identical thresholds. TFRC-214 has no
RefSeq equivalent at all, so Primer-BLAST reports "no target templates found" for a
perfectly good pair.
 
---
 
## Limitations
 
- **Annotation-bound.** An isoform absent from Ensembl is invisible to the screen.
- **Two databases.** GENCODE here, RefSeq at Primer-BLAST. See above.
- **Dense annotation blocks designs.** HNRNPC has 140 annotated isoforms, TFRC 61.
  Designing against all of them is conservative and means many targets return
  `no_unique_region`.
- **The top-10 collapses.** Several reported pairs are often the same assay with one
  or two bases added at the 5′ end. Raising `TOP_N` helps; dedup by anchor position
  would be better.
- **The margin number is pessimistic.** It reports 0 for pairs that brute force shows
  have no binding site at all.
- **These are predictions.** In-silico PCR models sequence, not polymerase, cycling
  conditions or template structure. Wet-lab validation is still required.
---
 
## References
 
- Kwok S et al. (1990) *Nucleic Acids Res* 18:999–1005 — effects of primer-template
  mismatches; one 3′-terminal mismatch is often tolerated, two block.
- Christopherson C et al. (1997) *Nucleic Acids Res* 25:654–658 — 2–4 mismatches
  barely affect yield, 5 reduces it ~22-fold.
- Stadhouders R et al. (2010) *J Mol Diagn* 12:109–117 — the 5 nt 3′ window.
- Ye J et al. (2012) *BMC Bioinformatics* 13:134 — Primer-BLAST; the off-target rule
  isosolo matches.
- Britten RJ et al. (1988) *PNAS* 85:4770–4774 — Alu is over 10% of the human genome.
- Song Y-C et al. (2023) *New Phytologist* 240:928–939 — QuantAS; independently
  arrives at junction anchoring, and the source of the combination-assay idea.
- McPherson MJ & Møller SG (2000) *PCR*, BIOS Scientific — the 5′ end does not
  initiate extension.
Built on primer3-py (SantaLucia 1998 thermodynamics), NCBI BLAST+, the Ensembl REST
API and GENCODE v50.
 
## License
 
MIT
