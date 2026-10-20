# isosolo

Design PCR primer pairs that amplify **one transcript isoform** and not its siblings.

---

# Getting Started

A complete guide for someone who has just downloaded this folder and never run it before.

## Before you begin

**What you need:**

- A computer running macOS, Linux, or Windows with WSL
- An internet connection (the tool fetches sequences from Ensembl)
- About 15 minutes for first-time setup

**What you don't need:**

- Any prior Python knowledge
- A high-performance cluster — a laptop is fine
- Your own sequence files

**What this tool does:** You give it a gene name and one transcript ID. It gives you a PCR primer pair that amplifies that transcript and not the other versions of the same gene.

## Part 1 — Setup (once)

### 1.1 Open a terminal

- **macOS:** press `Cmd + Space`, type "Terminal", press Enter
- **Linux:** `Ctrl + Alt + T`
- **Windows:** open "Ubuntu" or "WSL" from the Start menu

A window with a text prompt appears. Every command below is typed there, one line at a time, pressing Enter after each.

### 1.2 Go to the downloaded folder

If you downloaded a ZIP from GitHub, unzip it first by double-clicking. Then:

```bash
cd ~/Downloads/isosolo
```

Adjust the path if you put it elsewhere. Check you're in the right place:

```bash
ls
```

You should see `pyproject.toml`, `README.md`, and a folder called `src`. If you see something else, you're in the wrong directory, or the real folder is nested one level deeper — try `cd isosolo` again.

### 1.3 Check you have Python

```bash
python3 --version
```

You need **3.9 or newer**. If you get "command not found," install Python from [python.org](https://www.python.org/downloads/) first.

### 1.4 Install

```bash
pip install -e .
```

This reads `pyproject.toml` and pulls in the two things isosolo needs: `primer3-py` for melting-temperature calculations and `requests` for fetching sequences.

Takes a minute or two. Lots of text scrolls past — normal.

**If `pip` isn't found:** try `python3 -m pip install -e .`

**If you get a permissions error:** use a virtual environment instead (recommended anyway):

```bash
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -e .
```

With a venv you must run `source venv/bin/activate` each time you open a new terminal.

### 1.5 Confirm it worked

```bash
isosolo --version
```

```
isosolo 0.1.0
```

If you see that, you're ready.

## Part 2 — Your first run

### 2.1 Look at a gene

Start here, always. It's fast and tells you whether everything is working.

```bash
isosolo check --gene LRRTM4
```

```
LRRTM4: 8 transcript(s)
transcript            name                  length
ENST00000409911       LRRTM4-205              3063
ENST00000409884       LRRTM4-204              3198
ENST00000409093       LRRTM4-202              3120
ENST00000409088       LRRTM4-201              3616
ENST00000409282       LRRTM4-203              2633
ENST00000456154       LRRTM4-206               485
ENST00000491242       LRRTM4-207               596
ENST00001115845       LRRTM4-208              3966
```

Every transcript (isoform) of the gene, with its ID, its name, and its length in bases.

**Pick the one you want to target** and copy its transcript ID from the left column.

*If this hangs or errors, you have no internet connection, or Ensembl is temporarily down.*

### 2.2 Design primers

```bash
isosolo design \
  --gene LRRTM4 \
  --transcript ENST00000456154 \
  --outdir ~/Desktop/my_primers \
  --no-primersearch \
  -v
```

> **Typing long commands:** the `\` at the end of a line means "continues on the next line." You can type it all on one line without the backslashes if you prefer.

Breaking it down:

| Part | Meaning |
|---|---|
| `design` | the thing you want to do |
| `--gene LRRTM4` | which gene |
| `--transcript ENST00000456154` | which isoform to target |
| `--outdir ~/Desktop/my_primers` | where to put the results |
| `--no-primersearch` | skip the slow optional checks for now |
| `-v` | show detailed progress |

Takes 1–3 minutes depending on transcript length.

### 2.3 Read the output

```
======================================================================
LRRTM4 | ENST00000456154
target LRRTM4-206 (ENST00000456154), 485 nt
1 exon-exon junction(s)
7 sibling isoform(s)
uniqueness detection 0.3s
3 region(s) -> 1 interval(s), 36 bp unique
41356 thermodynamically valid candidate(s) in 52.3s
10/10 reported candidates clean on every check
======================================================================
batch complete: {'ok': 1}
```

Line by line:

- **`485 nt`** — your target is 485 bases long
- **`7 sibling isoform(s)`** — the 7 other versions of this gene it must avoid
- **`36 bp unique`** — **the most important number.** This much of your target exists nowhere in the siblings. Zero would mean no assay is possible.
- **`41356 candidates`** — possible primer pairs that passed the chemistry checks
- **`10/10 clean`** — of the 10 reported, all amplify your target and none of the siblings

## Part 3 — Getting your primers

### 3.1 Find the results

```bash
cd ~/Desktop/my_primers/LRRTM4-206
ls
```

```
LRRTM4-206_assay.csv        ← your primers
LRRTM4-206_assay.html       ← visual check
LRRTM4-206_regions.html     ← what was unique
LRRTM4-206_regions.tsv
LRRTM4-206_regions.txt
LRRTM4-206_positions.csv
```

### 3.2 Open the CSV

Double-click `LRRTM4-206_assay.csv` to open it in Excel or Numbers. Each row is a primer pair, **best first**.

Or from the terminal:

```bash
column -s, -t < LRRTM4-206_assay.csv | less -S
```

Press `q` to quit, arrow keys to scroll.

### 3.3 The columns that matter

Ignore most of them. These six decide whether a pair is good:

| Column | What you want | Why |
|---|---|---|
| `Sibling_Amplicons` | **0** | how many other isoforms this pair would amplify |
| `Anchor_Class` | `chimeric_region` | the strong kind of uniqueness; `base_difference` is weaker |
| `Min_3p_Mismatches_vs_Any_Isoform` | high (max 5) | how different the discriminating end is from the siblings |
| `Spans_Exon_Junction` | `yes` | means contaminating genomic DNA can't amplify |
| `Amplicon_Size_bp` | 70–150 for ddPCR | product length |
| `PS_Offgene_Amplicons` | **0** | other genes hit (blank if you used `--no-primersearch`) |

Then take `Forward_Primer_Sequence_5to3` and `Reverse_Primer_Sequence_5to3` from the best row.

### 3.4 Look at it

Double-click `LRRTM4-206_assay.html`. It opens in your browser and shows, for each candidate, two aligned sequence tracks:

- **top:** the full transcript with unique regions highlighted
- **bottom:** where the primers sit

Check that the anchored primer is sitting *on* a highlighted region, not next to it. This is the fastest sanity check and it takes five seconds.

## Part 4 — What the four outcomes mean

The run reports one of these:

**`ok`** — primers found.

**`no_unique_region`** — your target's sequence exists entirely inside another isoform. **No primer pair can tell them apart.** This isn't a bug or a failure of the tool; it's a fact about the biology. It happens when one isoform is just a shortened version of another. Pick a different target.

**`no_candidates`** — unique sequence exists, but no primer pair met the temperature, GC, or size requirements. Try relaxing them:

```bash
--min-amplicon 60 --max-amplicon 250 --tm-min 54 --tm-max 65
```

**`no_sequence`** — couldn't fetch from Ensembl. Check your gene name and transcript ID are correct and that you have internet.

## Part 5 — The optional extra checks

The run above only checked your target against **other isoforms of the same gene**. It did not check whether your primers hit some completely different gene.

Two ways to cover that.

### 5.1 NCBI Primer-BLAST (easiest, no setup)

1. Go to [ncbi.nlm.nih.gov/tools/primer-blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
2. Paste your forward sequence into "Forward primer"
3. Paste your reverse sequence into "Reverse primer"
4. Set Organism to *Homo sapiens*, Database to `refseq_rna`
5. Click "Get Primers"

You want it to report hits only on your target transcript.

This is the most reliable check and needs nothing installed.

### 5.2 Local PrimerSearch (optional, advanced)

Requires EMBOSS and a transcriptome FASTA file:

```bash
conda install -c bioconda emboss
```

Then:

```bash
isosolo design --gene LRRTM4 --transcript ENST00000456154 \
  --outdir ~/Desktop/my_primers \
  --transcriptome /path/to/gencode.transcripts.fa \
  --primersearch-mismatch 10 \
  -v
```

**Be aware: this is slow** — roughly 60 seconds per primer pair against a full transcriptome, and it checks 30 pairs by default. Budget half an hour.

To speed it up, reduce how many pairs get checked:

```bash
sed -i 's/^PRIMERSEARCH_SCREEN_TOP = 30/PRIMERSEARCH_SCREEN_TOP = 8/' src/isosolo/config.py
```

## Part 6 — Several targets at once

Make a file called `targets.tsv`. Two columns, separated by a **tab**:

```
gene	transcript
LRRTM4	ENST00000456154
BIRC3	ENST00000532808
GLUL	ENST00000331872
```

Then:

```bash
isosolo batch --targets targets.tsv --outdir ~/Desktop/batch_run -v
```

Each target gets its own folder, plus a `summary.csv` at the top level with one row per target. If one gene fails, the rest keep going.

## Part 7 — How it actually works

Optional reading, but it helps you judge the output.

### The problem

You want to amplify isoform 206 and nothing else. But isoforms share most of their sequence — they're different arrangements of the same exons. So the question is: **what part of 206 exists nowhere else?**

### Step 1 — Find unique sequence

Two detectors, because isoforms differ in two different ways.

**Detector A — is this stretch new?**

Slide a 41-base window along your target. At each position, search for that exact string in every sibling. If it's found nowhere, that position is unique.

This catches **junctions**. Suppose your target splices exon 1 directly to exon 3. The siblings have exon 1 and they have exon 3, but with exon 2 wedged between. So a window straddling your junction is half one shared piece, half another shared piece — and that *combination* appears nowhere. The window is a chimera. Hence the label `chimeric_region`.

**Detector B — is this one letter different?**

For single-base differences, Detector A is too coarse. Detector B first finds where the equivalent region sits in each sibling (sliding a 15-base window, keeping matches ≥90% identical), then compares that single base.

Catches splice sites shifted by a few bases, and short insertions or deletions.

### Step 2 — Anchor the 3′ end

A primer has two ends. The **3′ end** is where the enzyme starts copying.

A primer that mismatches at its 5′ end still works. A primer that mismatches in its last 2–3 bases usually doesn't work at all.

So the discriminating base must be the 3′ terminal base. Merely *overlapping* the unique region isn't enough — if the 3′ end sits in shared sequence, the primer works on every isoform.

### Step 3 — Check it doesn't amplify anything else

Simulates PCR against each sibling. A primer is assumed to work on a sibling only if its last 5 bases match perfectly and total mismatches stay low. If both primers would bind at a plausible distance, that's a predicted off-target and the pair is penalised.

### Step 4 — Rank

Specificity first, chemistry second. A thermodynamically perfect pair that hits another isoform loses to a merely decent pair that doesn't.

## Part 8 — Worked example

**LRRTM4-206** is 485 bases: exon 1 (382 bp) joined to exon 2 (103 bp). On the genome those exons sit 72,000 bases apart.

Neither exon is unique:

- exon 1 sits entirely inside LRRTM4-207's longer first exon
- exon 2 sits entirely inside the large exon shared by five other isoforms

A primer anywhere in exon 1 also hits 207. Anywhere in exon 2, it hits five others.

**But the join is unique.** Only in 206 do those two particular pieces sit adjacent. Detector A found a 36 bp zone centred on the junction.

The result:

```
Forward   GCTCACACTGAGAGCTGT          binds positions 317-334  (shared sequence)
Reverse   GCGTAATTAAATGGAAACCTGGG     binds positions 378-400  (across the join)
Amplicon  84 bp
```

The reverse primer needs 4 bases of exon 1 followed immediately by 19 bases of exon 2. In every other isoform, those 4 bases are followed by something completely different — so it can't bind.

Two bonuses fall out of this:

- **Genomic DNA can't amplify.** On the genome the primers sit 72,000 bases apart. PCR can't bridge that, so DNA contamination produces nothing.
- **PrimerSearch against the full human transcriptome returned exactly one hit** — the target, zero mismatches, 84 bp.

## Part 9 — Common problems

**`command not found: isosolo`**
The install didn't complete, or you're in a new terminal without the virtual environment active. Run `source venv/bin/activate`, or reinstall with `pip install -e .`

**`gene XXXX not found`**
Check the spelling. Use the official HGNC symbol (`LRRTM4`, not `lrrtm4` or an alias).

**It hangs at "fetching"**
No internet, or Ensembl is down. Wait and retry. Sequences are cached, so a retry skips what already downloaded.

**`0 bp unique`**
Your target is contained inside another isoform. No tool can design primers for it. Pick a different target.

**`no_candidates` despite unique sequence found**
The unique region is too small or awkwardly placed for the size and temperature constraints. Relax them with `--min-amplicon 60 --max-amplicon 250`.

**It's very slow**
Long transcripts with many siblings take longer. Skip `--transcriptome` and `--genome` for a first pass — they're the slow parts by far.

## Part 10 — Before you order primers

1. Open `*_assay.html` and confirm the anchored primer sits on a highlighted unique region.
2. Check `Sibling_Amplicons` is **0** on your chosen row.
3. Run NCBI Primer-BLAST on the pair (Part 5.1).
4. Save the folder — `*_positions.csv` is your audit trail if anyone asks why a position was called unique.

## What this tool cannot do

- **It only knows annotated isoforms.** If a transcript isn't in Ensembl, isosolo can't avoid it.
- **It's a prediction.** Primers must still be validated at the bench.
- **A `base_difference` anchor is weak.** A single-base difference gives limited discrimination. Prefer `chimeric_region`.
- **It won't design gene-level primers.** If you want to amplify *all* isoforms at once, that's the opposite problem and a different tool.

---

# Operating Guide

Reference documentation for routine use.

## 1. What it does

Designs PCR primer pairs that amplify **one transcript isoform** and not the other isoforms of the same gene.

Most primer tools solve the opposite problem. Primer3 designs on a sequence you hand it and has no concept of isoforms. Primer-BLAST checks genome-wide specificity but happily returns a pair hitting every transcript of your gene. IsoPrimer deliberately maximises coverage across isoforms.

isosolo has an **exclusion** objective: amplify the target, exclude the siblings.

## 2. Installation

```bash
cd /scratch/chikatipalli.t
tar -xzf isosolo.tar.gz
cd isosolo

conda activate isoprimer       # already has primer3 and EMBOSS
pip install -e .

isosolo --version              # isosolo 0.1.0
```

EMBOSS is only needed for the optional transcriptome and genome checks. The core design runs without it.

## 3. Inspect a gene first

```bash
isosolo check --gene LRRTM4
```

```
LRRTM4: 8 transcript(s)
transcript            name                  length
ENST00000409911       LRRTM4-205              3063
ENST00000456154       LRRTM4-206               485
...
```

Confirms the gene resolves, lists every isoform with its cDNA length, and verifies Ensembl is reachable. Always run this before designing.

## 4. Design

```bash
isosolo design \
  --gene LRRTM4 \
  --transcript ENST00000456154 \
  --outdir /scratch/chikatipalli.t/isosolo_out \
  --no-primersearch \
  -v
```

| Flag | Meaning |
|---|---|
| `--gene` | gene symbol |
| `--transcript` | Ensembl transcript ID, version optional |
| `--outdir` | where results go |
| `--no-primersearch` | skip the external checks for a fast first pass |
| `-v` | verbose log |

Runtime: roughly 1–2 minutes for a short transcript.

## 5. What happens inside, step by step

### Step 1 — Fetch sequences

Pulls the target's cDNA plus every sibling isoform's cDNA from Ensembl REST. Cached under `~/.ensembl_cache`, so reruns are instant.

`cdna` is the spliced mature transcript with UTRs — what your reverse transcriptase actually sees.

### Step 2 — Find what's unique to the target

Two independent detectors, because isoforms differ in two different ways.

**Detector A — structural uniqueness (41 bp window)**

At every position, take the 41-base window centred there and ask whether that exact string exists in any sibling. Not found anywhere → unique.

This catches junctions. At a junction unique to your target, the window's left half exists in the siblings and its right half exists in the siblings, but never adjacent. The window is a *chimera* of two shared pieces, and the combination is new. Hence the class name `chimeric_region`.

Also catches unique exons, alternative first exons, and alternative 3′ ends.

**Detector B — base-level uniqueness (15 bp window, 90% identity)**

For single-base differences, Detector A is too coarse. B first *locates* the equivalent region in each sibling by sliding a 15 bp window and keeping the best alignment ≥90% identical, then compares that one base.

The 90% floor matters: without it, a random 15-mer elsewhere in a 4 kb sibling could carry the right base by chance and falsely clear your position.

Catches splice acceptors shifted a few bases, short indels.

**Output:** contiguous intervals, each labelled `chimeric_region` or `base_difference`. Chimeric is the stronger anchor.

### Step 3 — Anchor a primer's 3′ end inside a unique interval

PCR specificity is concentrated at the 3′ terminus, where polymerase begins extending. A primer mismatched at its 5′ end works fine; one mismatched in its last 2–3 bases usually won't extend at all.

So "overlaps the unique region" is not enough — the **terminal base** has to be in it.

Both orientations are tried:

- `fwd-anchored` — forward primer's 3′ end in the interval, reverse downstream
- `rev-anchored` — reverse primer's 3′ end in the interval, forward upstream

Note: a reverse primer binds the sense strand backwards, so its 3′-terminal base corresponds to the **leftmost** base of its footprint on the cDNA.

Each candidate is then filtered on GC 35–65%, Tm 56–63 °C, Tm balance ≤5 °C, amplicon 70–150 bp, and hairpin/dimer ΔG.

### Step 4 — Screen off-targets

**Sibling isoforms (always runs).** Simulates PCR against every other transcript of the gene. A primer is assumed to extend on a sibling only if its 3′-terminal 5 bases match perfectly *and* total mismatches are ≤3. If both primers would bind, in the right orientation, at a plausible distance, that's a predicted off-target product.

**Transcriptome (optional, `--transcriptome`).** EMBOSS PrimerSearch against all human transcripts. Catches paralogs and pseudogenes — things the sibling check is completely blind to.

**Genome (optional, `--genome`).** PrimerSearch against GRCh38. Catches amplification from genomic DNA contamination.

The two kinds of check do different jobs:

| | Python screen | PrimerSearch |
|---|---|---|
| Scope | sibling isoforms | whole transcriptome + genome |
| 3′-end aware | **yes** | no — flat mismatch % |
| Catches paralogs | no | **yes** |
| Catches gDNA products | no | **yes** |

### Step 5 — Rank

Specificity dominates chemistry. A thermodynamically perfect pair that hits a paralog loses to a merely acceptable pair that doesn't.

```
off-gene amplicons      (PrimerSearch, transcriptome-wide)
sibling amplicons       (Python screen)
3′ mismatch margin      (larger is better)
anchor class            (chimeric_region over base_difference)
junction-spanning       (avoids gDNA amplification)
genomic amplicons       (PrimerSearch, genome-wide)
thermodynamics          (dimers, hairpins, size, Tm balance)
```

## 6. Reading the log

```
LRRTM4 | ENST00000456154
target LRRTM4-206 (ENST00000456154), 485 nt
1 exon-exon junction(s)
7 sibling isoform(s)
uniqueness detection 0.3s
3 region(s) -> 1 interval(s), 36 bp unique          ← the key line
41356 thermodynamically valid candidate(s) in 52.3s
10/10 reported candidates clean on every check
```

**`N bp unique`** is the number that matters. Zero means no assay is possible. A few dozen means one anchor point. Hundreds means plenty of options.

## 7. Output files

| File | Contents |
|---|---|
| `*_assay.csv` | ranked primer pairs, ready to order |
| `*_assay.html` | dual sequence tracks — unique regions above, primer placement below |
| `*_regions.html` | full cDNA with unique regions colour-coded |
| `*_regions.tsv` | each region with 50 bp of flanking context |
| `*_positions.csv` | per-base audit: every position against every sibling |
| `summary.csv` | one row per target (batch mode) |

**Key CSV columns:**

| Column | Want |
|---|---|
| `Anchor_Class` | `chimeric_region` |
| `Anchor_Orientation` | either |
| `Sibling_Amplicons` | **0** |
| `Min_3p_Mismatches_vs_Any_Isoform` | high (max = TERM_WINDOW, default 5) |
| `Spans_Exon_Junction` | `yes` (blocks gDNA) |
| `Amplicon_Size_bp` | 70–150 for ddPCR |
| `PS_Offgene_Amplicons` | **0** |

The `*_positions.csv` audit trail is worth knowing about. If you ask "why was position 1,847 called unique," it tells you: sibling 201 aligned at 93% with a C, sibling 202 at 96% with a G, sibling 203 had no equivalent region. Every call is traceable.

## 8. Four outcomes

`status` in `summary.csv`:

- **`ok`** — candidates found
- **`no_unique_region`** — the target's sequence is contained within a sibling. **No primer pair can distinguish them.** This is the alternative-TSS / alternative-polyA case and it's a biological verdict, not a failure.
- **`no_candidates`** — unique sequence exists but nothing met the Tm/GC/amplicon constraints. Widen them.
- **`no_sequence`** — fetch failed.

## 9. Adding the external checks

PrimerSearch takes about **60 seconds per primer pair** against the full transcriptome, and it processes pairs sequentially. Screening 30 pairs = 30 minutes. Reduce first:

```bash
sed -i 's/^PRIMERSEARCH_SCREEN_TOP = 30/PRIMERSEARCH_SCREEN_TOP = 8/' src/isosolo/config.py
```

Then:

```bash
isosolo design --gene LRRTM4 --transcript ENST00000456154 \
  --outdir /scratch/chikatipalli.t/isosolo_out \
  --transcriptome /scratch/chikatipalli.t/refs/gencode.v50.transcripts.clean.fa \
  --primersearch-mismatch 10 \
  -v
```

Use 10% rather than the default 20% — at 20%, a 20 nt primer tolerates 4 mismatches anywhere, producing large amounts of meaningless output.

**Faster alternative:** run the design with `--no-primersearch`, then check only your chosen pair manually:

```bash
printf "mypair\tFORWARD_SEQ\tREVERSE_SEQ\n" > ps_in.txt

primersearch -seqall /path/transcriptome.fa -infile ps_in.txt \
  -mismatchpercent 10 -outfile ps_tx.txt -auto

grep "Sequence:" ps_tx.txt | sort -u
```

One hit, your target → clean.

The genome pass is much slower (GRCh38 is ~20× the transcriptome). Run it on a compute node.

## 10. Batch mode

`targets.tsv`:

```
gene	transcript
LRRTM4	ENST00000456154
BIRC3	ENST00000532808
GLUL	ENST00000331872
```

```bash
isosolo batch --targets targets.tsv \
  --outdir /scratch/chikatipalli.t/isosolo_batch -v
```

Writes `summary.csv` with one row per target. Exceptions in one gene don't kill the batch.

## 11. Configuration

Three levels, increasing precedence: `src/isosolo/config.py` → `--config file.yaml` → CLI flag.

```bash
--min-amplicon 70 --max-amplicon 150     # ddPCR sizing
--tm-min 56 --tm-max 63 --tm-balance 5
--chimeric-window 41                      # structural detection
--align-window 15 --align-identity 0.90   # base-level detection
--term-window 5 --max-mismatches 3        # specificity model
--top-n 10
```

`isosolo design --help` lists everything.

## 12. Worked example — LRRTM4-206

**The problem.** 206 is two exons, 382 bp + 103 bp. Neither exon belongs to it alone:

- exon 1 sits entirely inside LRRTM4-207's longer first exon
- exon 2 sits entirely inside the large exon shared by 201, 202, 204, 205, 208

A primer anywhere in exon 1 also hits 207. Anywhere in exon 2, it hits five others.

**What's unique.** The join. Only in 206 do those two specific pieces sit adjacent. Detector A found a 36 bp zone centred on the junction at cDNA position 382.

**The design.**

```
Forward   GCTCACACTGAGAGCTGT          binds 317-334   (exon 1, shared)
Reverse   GCGTAATTAAATGGAAACCTGGG     binds 378-400   (across the join)
Amplicon  84 bp
```

The reverse primer needs 4 bases of exon 1 followed immediately by 19 bases of exon 2. In every other isoform those 4 bases are followed by something different, so it can't bind.

**Validation.**

- 7/7 sibling isoforms excluded, 3′ mismatch margin 4/5
- PrimerSearch, full GENCODE v50 transcriptome, 10% mismatch: one hit, `ENST00000456154.1`, 0 mismatches, 84 bp
- Junction spans a 72 kb intron → genomic DNA can't amplify

**Note:** IsoPrimer cannot design for this target. 206 is TSL 4 and IsoPrimer filters to TSL ≤ 3 before Primer3 runs.

## 13. Known limitations

- **Annotation-bound.** An isoform absent from Ensembl is invisible.
- **PrimerSearch is not 3′-aware.** It applies one mismatch percentage across the whole primer, so it over-calls off-targets that wouldn't extend. Check *where* reported mismatches fall.
- **Single-base anchors are weak.** A `base_difference` anchor gives roughly allele-specific-PCR discrimination, perhaps 1–2 logs. Prefer `chimeric_region`.
- **Candidate generation is inefficient.** A 36 bp interval produced 41,356 candidates in 52 s, most differing by one base. Needs deduplication by 3′ anchor position before the 166-gene panel is practical.
- **PrimerSearch is not batched or parallelised.** One sequential call. IsoPrimer splits into 4 parallel batches; isosolo should too.
- **Wet-lab validation still required.** This is in-silico prediction.

## 14. Before ordering

1. Open `*_assay.html` and confirm the anchored primer sits on the unique region, not drifting into shared sequence.
2. Run PrimerSearch against the transcriptome on your chosen pair.
3. Run NCBI Primer-BLAST against *Homo sapiens* refseq_rna — it queries live NCBI rather than your local snapshot, and it is 3′-end aware.
4. Save the validation output alongside the design.

## 15. Testing checklist

Worth running against targets with known structure to confirm behaviour:

| Case | Expected |
|---|---|
| Target with a unique junction | `chimeric_region` anchor, `spans_junction=yes` |
| Target with a unique exon | `chimeric_region`, possibly `spans_junction=no` |
| Target that is a strict truncation of a sibling | `no_unique_region` |
| Single-isoform gene | large unique region, many candidates |
| Target with only a few bases different | `base_difference` anchor, low mismatch margin |
| Same gene run twice | identical output (cache makes the second instant) |

---

## License

MIT
