#!/usr/bin/env python3
"""
isosolo - isoform-specific PCR primer design
============================================

Designs PCR primer pairs that amplify ONE transcript isoform and not the
other isoforms of the same gene.

Sequences come from the Ensembl REST API, so the basic run needs nothing
downloaded. You need an internet connection the first time you run a
gene; after that its sequences are cached in ~/.ensembl_cache and reruns
are instant.

    pip install primer3-py requests

OPTIONAL - transcriptome- and genome-wide specificity
    The sibling check below only looks inside the target's own gene. To
    also catch paralogs, pseudogenes and genomic-DNA products, install
    EMBOSS and point the script at reference FASTA files:

        conda install -c bioconda emboss

        transcriptome  https://www.gencodegenes.org/human/
                       "Transcript sequences" (~400 MB gzipped)
        genome         same page, "Genome sequence" (~900 MB gzipped)

    Then set TRANSCRIPTOME_FASTA / GENOME_FASTA below. Budget about a
    minute per primer pair for the transcriptome and three for the
    genome, which is why only the top few pairs are checked. Leave them
    as None and the script skips the step with a note.

Edit the TARGET_JOBS block below, then press Run.

HOW IT WORKS
    1. Fetch the target transcript and every sibling isoform of the gene.
    2. Find stretches of the target that appear in no sibling.
         chimeric_region  - a 41 bp window absent from every sibling.
                            Catches novel junctions, unique exons,
                            alternative first and last exons.
         base_difference  - a single base that differs after locating the
                            equivalent region by local alignment.
                            Catches shifted splice sites and short indels.
    3. Place a primer so its 3'-TERMINAL BASE sits inside one of those
       stretches. The 3' end is where polymerase starts extending, so a
       primer that merely overlaps a unique region still works on every
       sibling.
    4. Simulate PCR against every sibling. A primer is assumed to extend
       on a sibling only if its last 5 bases match perfectly and total
       mismatches stay low. Pairs where both primers bind at a plausible
       distance are flagged as off-target.
    5. OPTIONAL: EMBOSS PrimerSearch against the whole transcriptome
       and the genome. This is the complement of step 4 - it looks
       outside the gene, but applies a flat mismatch percentage rather
       than weighting the 3' end, so neither check replaces the other.
    6. Rank by specificity first, chemistry second.

OUTPUT, per target, in a folder named after the isoform:
    *_assay.csv       ranked primer pairs, ready to order
    *_assay.html      unique regions and primers on paired sequence tracks
    *_regions.html    full cDNA, unique regions colour-coded
    *_regions.tsv     each region with 50 bp of flanking context
    *_positions.csv   per-base audit: every position against every sibling

If it reports NO UNIQUE REGIONS, every part of the target also exists in
another isoform and no primer pair can tell them apart. That is a fact
about the transcripts, not a failure of the search. Try another isoform.

Prefer anchor_class = chimeric_region. A base_difference anchor
discriminates on one base, which is allele-specific-PCR territory (1-2
logs) and needs careful annealing optimisation at the bench.

Without the optional FASTA files, specificity here is WITHIN THE GENE
only. Either way, run the chosen pair through NCBI Primer-BLAST (Homo
sapiens, refseq_rna) before ordering - it queries live NCBI rather than a
local snapshot, and it weights 3' mismatches properly, which PrimerSearch
does not. This is in-silico prediction; wet-lab validation is still
required.
"""

import csv
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict

try:
    import primer3
except ImportError:
    sys.exit("isosolo needs primer3-py:\n    pip install primer3-py")

try:
    import requests
except ImportError:
    sys.exit("isosolo needs requests:\n    pip install requests")


# ==========================================
# CONFIGURATION
# ==========================================
TARGET_JOBS = [
    {"gene": "BIRC3", "BIRC3-201": "ENST00000263464"},
    {"gene": "BIRC3", "BIRC3-206": "ENST00000532808"}
]

OUTPUT_ROOT = os.path.expanduser("~/Desktop/BIRC3_unique_isoform_primer_design")

MAX_RETRIES = 5
TOP_N = 10
headers = {"Content-Type": "application/json"}
CACHE_DIR = os.path.expanduser("~/.ensembl_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# Amplicon (ddPCR-specific; NCBI defaults to 70-1000)
MIN_AMPLICON = 70
MAX_AMPLICON = 150
CONTEXT_FLANK = 50

# Uniqueness detection
CHIMERIC_WINDOW = 41
ALIGN_WINDOW = 15
ALIGN_IDENTITY = 0.90

# True  = a position is non-unique if ANY alignment above ALIGN_IDENTITY
#         carries the same base (conservative, the safer default)
# False = only the single best-scoring alignment is consulted
STRICT_BASE_DIFFERENCE = True

# Primer length
FORWARD_PRIMER_LEN_RANGE = (18, 24)
REVERSE_PRIMER_LEN_RANGE = (18, 24)

# Primer Tm constraints (match NCBI Primer-BLAST)
PRIMER_TM_MIN = 56.0
PRIMER_TM_MAX = 63.0
PRIMER_TM_TARGET = 60.0

# Hard limit on the Tm difference between forward and reverse primers
TM_BALANCE_MAX = 5.0

# Minimum primer bases overlapping a chimeric_region interval.
# base_difference intervals are usually 1-2 bp, so they are exempt: there
# the requirement is that the 3'-terminal base sits on the differing base.
MIN_PRIMER_OVERLAP = 10

# --- specificity model ---
# A primer is assumed to extend on an off-target only if its 3'-terminal
# window matches perfectly and total mismatches stay low.
TERM_WINDOW = 5                # 3' bases that must match to extend
MAX_TOTAL_MM = 3               # total mismatches tolerated across primer
MAX_OFFTARGET_PRODUCT = 1200   # bp; longer predicted products ignored
SEED_LEN = 6                   # k-mer seed for the mismatch search

# Running the screen on every thermodynamically valid candidate is
# wasteful; screen this many best-by-thermo instead.
SPECIFICITY_SCREEN_TOP = 400

# Prefer junction-spanning amplicons (genomic DNA cannot amplify them)
PREFER_JUNCTION_SPANNING = True

# --- OPTIONAL: EMBOSS PrimerSearch, transcriptome- and genome-wide ---
# Leave both FASTA paths as None to skip this step entirely.
# The screen above is within-gene and 3'-end aware; PrimerSearch is the
# complement - it sees the whole transcriptome and the genome, so it
# catches paralogs, pseudogenes and gDNA products that a sibling-only
# check cannot. It is NOT 3'-end aware, applying one mismatch percentage
# across the whole primer, so it over-reports. Both are kept because
# neither replaces the other.
TRANSCRIPTOME_FASTA = None
GENOME_FASTA = None
# TRANSCRIPTOME_FASTA = os.path.expanduser("~/ref/gencode.v47.transcripts.fa")
# GENOME_FASTA = os.path.expanduser("~/ref/GRCh38.primary_assembly.genome.fa")

PRIMERSEARCH_BIN = "primersearch"
# 10%, not the EMBOSS-ish default of 20. At 20% a 20 nt primer tolerates
# 4 mismatches anywhere, which never amplifies, and the output fills with
# hits that mean nothing.
PRIMERSEARCH_MISMATCH_PCT = 10
# Measured at ~55 s per pair against a full human transcriptome, and the
# pairs are processed sequentially - so 30 pairs is half an hour. Five
# keeps a single-gene run to about five minutes.
PRIMERSEARCH_SCREEN_TOP = 5
PRIMERSEARCH_KEEP_FILES = True     # keep working files for inspection


log = logging.getLogger("isosolo")
logging.basicConfig(level=logging.INFO, format="%(message)s",
                    stream=sys.stdout)


# ==========================================
# UTILITIES
# ==========================================
def revcomp(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def gc_percent(seq):
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in ("G", "C"))
    return round((gc / len(seq)) * 100, 2)


def tm(seq):
    try:
        return round(primer3.calc_tm(seq), 2)
    except Exception:
        return 0.0


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

# Names below describe what primer3-py actually computes. The previous
# version labelled calc_hairpin as SELF_END_TH and calc_end_stability as
# COMPL_END_TH; neither is correct.
def homodimer_dg(seq):
    """Self-complementarity, any position (~Primer3 SELF_ANY_TH)."""
    try:
        return round(primer3.calc_homodimer(seq).dg / 1000.0, 2)
    except Exception:
        return 0.0


def hairpin_dg(seq):
    """Hairpin stability (~Primer3 HAIRPIN_TH)."""
    try:
        return round(primer3.calc_hairpin(seq).dg / 1000.0, 2)
    except Exception:
        return 0.0


def heterodimer_dg(fwd, rev):
    """Primer-dimer, any position (~Primer3 COMPL_ANY_TH)."""
    try:
        return round(primer3.calc_heterodimer(fwd, rev).dg / 1000.0, 2)
    except Exception:
        return 0.0


def end_stability_dg(fwd, rev):
    """3'-end stability of fwd against rev (Primer3 END_STABILITY)."""
    try:
        return round(primer3.calc_end_stability(fwd, rev).dg / 1000.0, 2)
    except Exception:
        return 0.0

def build_kmer_index(seq, k):
    idx = defaultdict(list)
    for i in range(len(seq) - k + 1):
        idx[seq[i:i + k]].append(i)
    return idx


def build_kmer_set(seqs, k):
    s = set()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            s.add(seq[i:i + k])
    return s


# ==========================================
# SEQUENCE RETRIEVAL - Ensembl REST, with an on-disk cache
# ==========================================
headers = {"Content-Type": "application/json"}
os.makedirs(CACHE_DIR, exist_ok=True)

def _cache_path(key):
    safe = "".join(c if c.isalnum() else "_" for c in key)
    return os.path.join(CACHE_DIR, f"{safe}.json")


def fetch_json(url, description="data", cache_key=None):
    if cache_key:
        cp = _cache_path(cache_key)
        if os.path.exists(cp):
            try:
                with open(cp) as f:
                    return json.load(f)
            except Exception:
                pass
    backoff = 1.5
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            res = requests.get(url, headers=headers, timeout=30)
            if res.ok:
                data = res.json()
                if cache_key:
                    try:
                        with open(_cache_path(cache_key), "w") as f:
                            json.dump(data, f)
                    except Exception:
                        pass
                return data
            if res.status_code in (429, 500, 502, 503, 504):
                wait = backoff + random.uniform(0, 1.0)
                print(f"[!] {description} attempt {attempt}/{MAX_RETRIES} "
                      f"-> {res.status_code}; backing off {wait:.1f}s")
                time.sleep(wait)
                backoff *= 1.8
            else:
                print(f"[!] {description} -> {res.status_code}")
                return None
        except requests.exceptions.RequestException:
            wait = backoff + random.uniform(0, 1.0)
            print(f"[!] {description} attempt {attempt}/{MAX_RETRIES} "
                  f"network error; backing off {wait:.1f}s")
            time.sleep(wait)
            backoff *= 1.8
    print(f"[!] {description}: exhausted retries.")
    return None


_GENE_CACHE = {}


def fetch_gene_transcripts(gene):
    if gene in _GENE_CACHE:
        return _GENE_CACHE[gene]
    sym_url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
               f"?content-type=application/json")
    gene_info = fetch_json(sym_url, f"symbol lookup {gene}",
                           cache_key=f"symbol_{gene}")
    if not gene_info or "id" not in gene_info:
        return None
    gene_id = gene_info["id"]
    exp_url = (f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1;"
               f"content-type=application/json")
    expanded = fetch_json(exp_url, f"expanded gene {gene_id}",
                          cache_key=f"gene_expand_{gene_id}")
    if not expanded or "Transcript" not in expanded:
        return None
    transcripts = expanded["Transcript"]
    names = {}
    for t in transcripts:
        tid = t["id"].split(".")[0]
        names[tid] = t.get("display_name") or tid
    _GENE_CACHE[gene] = {"transcripts": transcripts, "names": names}
    return _GENE_CACHE[gene]


def fetch_exon_junctions(tid):
    """cDNA offsets of exon-exon junctions for a transcript.

    Returns a sorted list of integers; offset j means the junction falls
    between cDNA positions j-1 and j. Empty list on failure or for a
    single-exon transcript.
    """
    url = (f"https://rest.ensembl.org/lookup/id/{tid}?expand=1;"
           f"content-type=application/json")
    data = fetch_json(url, f"exons {tid}", cache_key=f"tx_expand_{tid}")
    if not data or "Exon" not in data:
        return []
    exons = data["Exon"]
    strand = data.get("strand", 1)
    # Ensembl lists exons in genomic order; transcript order is reversed
    # on the minus strand.
    ordered = sorted(exons, key=lambda e: e["start"],
                     reverse=(strand == -1))
    junctions = []
    cum = 0
    for e in ordered[:-1]:
        cum += (e["end"] - e["start"] + 1)
        junctions.append(cum)
    return junctions


def fetch_ncbi_cdna(transcript_id):
    base_id = transcript_id.split(".")[0]
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           f"db=nuccore&id={base_id}&rettype=fasta&retmode=text")
    try:
        res = requests.get(url, timeout=15)
        if res.ok and ">" in res.text:
            lines = res.text.strip().splitlines()
            return "".join(l.strip() for l in lines if not l.startswith(">"))
    except Exception:
        pass
    return None


def _load_cached_cdna(tid):
    cp = _cache_path(f"cdna_{tid}.txt")
    if not os.path.exists(cp):
        return None
    try:
        with open(cp) as f:
            raw = f.read().strip()
    except Exception:
        return None
    if not raw:
        return None
    if raw[0] == '"' and raw[-1] == '"':
        try:
            return json.loads(raw)
        except Exception:
            pass
    return raw


def fetch_cdna(tid):
    s = _load_cached_cdna(tid)
    if s:
        return s
    url = (f"https://rest.ensembl.org/sequence/id/{tid}"
           f"?content-type=application/json;type=cdna")
    data = fetch_json(url, f"cDNA {tid}", cache_key=f"cdna_{tid}")
    seq = data.get("seq", "") if data else ""
    if not seq:
        seq = fetch_ncbi_cdna(tid) or ""
    if seq:
        try:
            with open(_cache_path(f"cdna_{tid}.txt"), "w") as f:
                f.write(seq)
        except Exception:
            pass
    return seq


# ==========================================
# UNIQUENESS DETECTION
# ==========================================
def chimeric_region_flags(target_cdna, other_cdnas, window=None):
    """Position is flagged when its surrounding window appears verbatim
    in no off-target. Identical semantics to the original; the set
    lookup replaces an O(n*m*L) substring scan."""
    window = CHIMERIC_WINDOW if window is None else window
    t = target_cdna.upper()
    others = [o.upper() for o in other_cdnas if o]
    n = len(t)
    half = window // 2
    flags = [False] * n

    full_kmers = build_kmer_set(others, window) if others else set()

    for p in range(n):
        ws = max(0, p - half)
        we = min(n, p + half + 1)
        w = t[ws:we]
        if len(w) == window:
            shared = w in full_kmers
        else:
            # truncated edge window: fall back to substring search
            shared = any(w in o for o in others)
        flags[p] = not shared
    return flags


def base_difference_flags(target_cdna, other_cdnas,
                          window=None, identity=None, strict=None):
    """Position is flagged when no off-target carries the same base in a
    high-identity local alignment.

    Acceleration: at identity 0.90 over a 15 bp window only one mismatch
    is tolerated, so by pigeonhole at least one of the window's two
    half-length seeds must match exactly. Seeding on those two k-mers
    turns the full scan into a handful of verifications per position.
    """
    window = ALIGN_WINDOW if window is None else window
    identity = ALIGN_IDENTITY if identity is None else identity
    strict = STRICT_BASE_DIFFERENCE if strict is None else strict
    t = target_cdna.upper()
    others = [o.upper() for o in other_cdnas if o]
    n = len(t)
    half = window // 2
    unique = [True] * n
    aligned = [{} for _ in range(n)]

    seed_k = max(4, half)
    indices = [build_kmer_index(o, seed_k) for o in others]

    max_mm = int(window - (identity * window) + 1e-9)
    seeded_ok = max_mm <= 1   # pigeonhole guarantee only holds here

    for p in range(n):
        ws = max(0, p - half)
        we = min(n, p + half + 1)
        w = t[ws:we]
        wlen = len(w)
        pos_in_window = p - ws

        for idx, o in enumerate(others):
            if len(o) < wlen:
                continue

            # ---- collect candidate offsets -------------------------
            if seeded_ok and wlen == window:
                offsets = set()
                for seed_off in (0, wlen - seed_k):
                    seed = w[seed_off:seed_off + seed_k]
                    for hit in indices[idx].get(seed, ()):
                        start = hit - seed_off
                        if 0 <= start <= len(o) - wlen:
                            offsets.add(start)
            else:
                offsets = range(len(o) - wlen + 1)

            # ---- evaluate ------------------------------------------
            best_ident = 0.0
            best_base = None
            same_base_hit = False
            for j in offsets:
                sub = o[j:j + wlen]
                matches = 0
                for a, b in zip(w, sub):
                    if a == b:
                        matches += 1
                ident = matches / wlen
                if ident > best_ident:
                    best_ident = ident
                    best_base = sub[pos_in_window]
                if ident >= identity and sub[pos_in_window] == t[p]:
                    same_base_hit = True
                    if not strict:
                        break

            aligned[p][idx] = (best_ident, best_base)

            if strict:
                if same_base_hit:
                    unique[p] = False
            else:
                if best_ident >= identity and best_base == t[p]:
                    unique[p] = False

    return unique, aligned


def merge_runs(chimeric_flags, base_flags):
    n = len(chimeric_flags)
    classes = [""] * n
    for i in range(n):
        if base_flags[i]:
            classes[i] = "base_difference"
        elif chimeric_flags[i]:
            classes[i] = "chimeric_region"

    regions = []
    i = 0
    while i < n:
        cur = classes[i]
        if cur == "":
            i += 1
            continue
        s = i
        while i < n and classes[i] == cur:
            i += 1
        regions.append({
            "start": s,
            "end": i,
            "length": i - s,
            "definition": cur,
        })
    return regions


def merge_into_uniqueness_intervals(regions):
    """Merge overlapping/adjacent regions, keeping the strongest class.

    Returns (start, end, cls) tuples. 'chimeric_region' wins over
    'base_difference' when a merged interval contains both, because a
    chimeric anchor gives far better discrimination.
    """
    if not regions:
        return []
    sorted_regions = sorted(regions, key=lambda r: r["start"])
    merged = []
    cur_s = sorted_regions[0]["start"]
    cur_e = sorted_regions[0]["end"]
    cur_c = {sorted_regions[0]["definition"]}
    for r in sorted_regions[1:]:
        if r["start"] <= cur_e:
            cur_e = max(cur_e, r["end"])
            cur_c.add(r["definition"])
        else:
            merged.append((cur_s, cur_e, _class_of(cur_c)))
            cur_s, cur_e, cur_c = r["start"], r["end"], {r["definition"]}
    merged.append((cur_s, cur_e, _class_of(cur_c)))
    return merged


def _class_of(classes):
    return ("chimeric_region" if "chimeric_region" in classes
            else "base_difference")


# ==========================================
# ASSAY DESIGN
# ==========================================

def _overlap_ok(p_start, p_end, u_s, u_e, cls):
    """chimeric intervals need real overlap; base_difference intervals
    are often 1 bp, so the 3'-end anchoring check carries them."""
    if cls == "base_difference":
        return True
    ov = max(0, min(p_end, u_e) - max(p_start, u_s))
    return ov >= MIN_PRIMER_OVERLAP


def design_assay(target_cdna, uniqueness_intervals,
                 min_amp=None, max_amp=None):
    min_amp = MIN_AMPLICON if min_amp is None else min_amp
    max_amp = MAX_AMPLICON if max_amp is None else max_amp
    t = target_cdna.upper()
    n = len(t)
    if not uniqueness_intervals:
        return []

    fmin, fmax = FORWARD_PRIMER_LEN_RANGE
    rmin, rmax = REVERSE_PRIMER_LEN_RANGE

    raw = []

    for (u_s, u_e, cls) in uniqueness_intervals:
        # --------------------------------------------------------------
        # Orientation 1: FORWARD primer. Its 3'-terminal base is the
        # rightmost base of the footprint, so index f3 must lie inside
        # the unique interval.
        # --------------------------------------------------------------
        for f3 in range(u_s, u_e):
            f_end = f3 + 1                      # exclusive
            for f_len in range(fmin, fmax + 1):
                f_start = f_end - f_len
                if f_start < 0:
                    continue
                fwd = t[f_start:f_end]
                if not _overlap_ok(f_start, f_end, u_s, u_e, cls):
                    continue
                f_gc = gc_percent(fwd)
                if not (35 <= f_gc <= 65):
                    continue
                f_tm = tm(fwd)
                if not (PRIMER_TM_MIN <= f_tm <= PRIMER_TM_MAX):
                    continue

                r_start_lo = max(f_end, f_start + min_amp - rmax)
                r_start_hi = min(n, f_start + max_amp - rmin + 1)
                for r_start in range(r_start_lo, r_start_hi):
                    for r_len in range(rmin, rmax + 1):
                        r_end = r_start + r_len
                        if r_end > n:
                            continue
                        prod = r_end - f_start
                        if not (min_amp <= prod <= max_amp):
                            continue
                        rev = revcomp(t[r_start:r_end])
                        c = _finish_candidate(
                            fwd, f_start, f_end, f_tm, f_gc, f_len,
                            rev, r_start, r_end, r_len,
                            prod, "fwd-anchored", f3, cls)
                        if c:
                            raw.append(c)

        # --------------------------------------------------------------
        # Orientation 2: REVERSE primer. FIXED. The reverse primer's
        # 3'-terminal base is the complement of the LEFTMOST base of its
        # sense-strand footprint, so it is r_start - not r_end - that
        # must lie inside the unique interval.
        # --------------------------------------------------------------
        for r3 in range(u_s, u_e):
            r_start = r3
            for r_len in range(rmin, rmax + 1):
                r_end = r_start + r_len
                if r_end > n:
                    continue
                rev = revcomp(t[r_start:r_end])
                if not _overlap_ok(r_start, r_end, u_s, u_e, cls):
                    continue
                r_gc = gc_percent(rev)
                if not (35 <= r_gc <= 65):
                    continue
                r_tm = tm(rev)
                if not (PRIMER_TM_MIN <= r_tm <= PRIMER_TM_MAX):
                    continue

                f_start_lo = max(0, r_end - max_amp)
                f_start_hi = max(0, r_end - min_amp + 1)
                for f_start in range(f_start_lo, f_start_hi):
                    for f_len in range(fmin, fmax + 1):
                        f_end = f_start + f_len
                        if f_end > r_start:
                            continue
                        prod = r_end - f_start
                        if not (min_amp <= prod <= max_amp):
                            continue
                        fwd = t[f_start:f_end]
                        f_gc = gc_percent(fwd)
                        if not (35 <= f_gc <= 65):
                            continue
                        f_tm = tm(fwd)
                        if not (PRIMER_TM_MIN <= f_tm <= PRIMER_TM_MAX):
                            continue
                        c = _finish_candidate(
                            fwd, f_start, f_end, f_tm, f_gc, f_len,
                            rev, r_start, r_end, r_len,
                            prod, "rev-anchored", r3, cls)
                        if c:
                            raw.append(c)

    seen = set()
    deduped = []
    for c in raw:
        key = (c["fwd"], c["rev"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(c)
    return deduped


def _finish_candidate(fwd, f_start, f_end, f_tm, f_gc, f_len,
                      rev, r_start, r_end, r_len,
                      prod, orientation, anchor_pos, cls):
    r_gc = gc_percent(rev)
    if not (35 <= r_gc <= 65):
        return None
    r_tm = tm(rev)
    if not (PRIMER_TM_MIN <= r_tm <= PRIMER_TM_MAX):
        return None
    if abs(f_tm - r_tm) > TM_BALANCE_MAX:
        return None

    homo_f = homodimer_dg(fwd)
    hairpin_f = hairpin_dg(fwd)
    homo_r = homodimer_dg(rev)
    hairpin_r = hairpin_dg(rev)
    hetero = heterodimer_dg(fwd, rev)
    end_stab = end_stability_dg(fwd, rev)

    if hetero < -5.0:
        return None
    if hairpin_f < -5.0 or hairpin_r < -5.0:
        return None

    return {
        "orientation": orientation,
        "anchor_class": cls,
        "anchor_position": anchor_pos,
        "fwd": fwd, "fwd_start": f_start, "fwd_end": f_end,
        "fwd_tm": f_tm, "fwd_gc": f_gc, "fwd_len": f_len,
        "rev": rev, "rev_start": r_start, "rev_end": r_end,
        "rev_tm": r_tm, "rev_gc": r_gc, "rev_len": r_len,
        "product_size": prod,
        "homo_f": homo_f, "hairpin_f": hairpin_f,
        "homo_r": homo_r, "hairpin_r": hairpin_r,
        "hetero": hetero, "end_stab": end_stab,
    }


def thermo_key(c):
    return (
        -c["hetero"],
        -(c["hairpin_f"] + c["hairpin_r"]),
        -(c["homo_f"] + c["homo_r"]),
        c["product_size"],
        abs(c["fwd_tm"] - c["rev_tm"]),
        abs(c["fwd_tm"] - PRIMER_TM_TARGET) +
        abs(c["rev_tm"] - PRIMER_TM_TARGET),
    )


def rank_candidates(candidates, junctions):
    for c in candidates:
        c["spans_junction"] = any(c["fwd_start"] < j < c["rev_end"]
                                  for j in junctions)

    def key(c):
        junction_penalty = 0
        if PREFER_JUNCTION_SPANNING and junctions:
            junction_penalty = 0 if c["spans_junction"] else 1
        class_penalty = 0 if c["anchor_class"] == "chimeric_region" else 1
        # PrimerSearch fields are absent until annotate_with_primersearch
        # has run, so this key works for both ranking passes.
        offgene = len(c.get("ps_offgene") or ())
        genome = len(c.get("ps_genome_hits") or ())
        return (
            offgene,                       # paralog/pseudogene amplicons
            c["n_offtarget"],              # sibling isoform amplicons
            -c["min_3p_mismatches"],       # discrimination margin
            class_penalty,                 # prefer chimeric anchors
            junction_penalty,              # prefer no-gDNA amplicons
            genome,                        # gDNA products
        ) + thermo_key(c)

    return sorted(candidates, key=key)


# ==========================================
# SPECIFICITY - in-silico PCR against the sibling isoforms
# ==========================================
def best_binding(footprint, three_prime_at_right, template, index):
    """Most dangerous binding site of a primer on an off-target.

    footprint             primer sequence as it maps onto the sense
                          strand of the template
    three_prime_at_right  True for a forward primer, False for a reverse
                          primer (whose 3' end maps to the LEFT edge)

    Returns (term_mm, total_mm, position) for the site that would be
    easiest to extend, or None if no seeded site exists.
    """
    L = len(footprint)
    if len(template) < L:
        return None

    # Three seeds: pigeonhole covers up to 2 mismatches across the primer.
    seed_offsets = {0, max(0, (L - SEED_LEN) // 2), L - SEED_LEN}
    offsets = set()
    for so in seed_offsets:
        seed = footprint[so:so + SEED_LEN]
        for hit in index.get(seed, ()):
            start = hit - so
            if 0 <= start <= len(template) - L:
                offsets.add(start)
    if not offsets:
        return None

    best = None
    for j in offsets:
        sub = template[j:j + L]
        total_mm = 0
        for a, b in zip(footprint, sub):
            if a != b:
                total_mm += 1
        if three_prime_at_right:
            term = sum(1 for a, b in zip(footprint[-TERM_WINDOW:],
                                         sub[-TERM_WINDOW:]) if a != b)
        else:
            term = sum(1 for a, b in zip(footprint[:TERM_WINDOW],
                                         sub[:TERM_WINDOW]) if a != b)
        key = (term, total_mm)
        if best is None or key < (best[0], best[1]):
            best = (term, total_mm, j)
    return best


def screen_specificity(cand, target_cdna, off_templates, off_names,
                       off_indices, binding_cache):
    """Attach off-target binding and in-silico PCR results to a candidate."""
    f_fp = cand["fwd"]                       # forward footprint = sense
    r_fp = revcomp(cand["rev"])              # reverse footprint = sense

    offtargets = []
    worst_discrimination = 99

    for k, tmpl in enumerate(off_templates):
        fk = (f_fp, k, True)
        if fk not in binding_cache:
            binding_cache[fk] = best_binding(f_fp, True, tmpl, off_indices[k])
        fb = binding_cache[fk]

        rk = (r_fp, k, False)
        if rk not in binding_cache:
            binding_cache[rk] = best_binding(r_fp, False, tmpl,
                                             off_indices[k])
        rb = binding_cache[rk]

        f_term = fb[0] if fb else TERM_WINDOW
        r_term = rb[0] if rb else TERM_WINDOW
        # The better-discriminating primer is what protects the assay.
        discrimination = max(f_term, r_term)
        worst_discrimination = min(worst_discrimination, discrimination)

        amplifies = False
        product = None
        if fb and rb:
            f_term_mm, f_total, f_pos = fb
            r_term_mm, r_total, r_pos = rb
            if (f_term_mm == 0 and r_term_mm == 0
                    and f_total <= MAX_TOTAL_MM
                    and r_total <= MAX_TOTAL_MM
                    and r_pos + len(r_fp) > f_pos):
                size = (r_pos + len(r_fp)) - f_pos
                if 0 < size <= MAX_OFFTARGET_PRODUCT:
                    amplifies = True
                    product = size

        if amplifies:
            offtargets.append({
                "isoform": off_names[k],
                "product_size": product,
                "fwd_total_mm": fb[1],
                "rev_total_mm": rb[1],
            })

    cand["offtarget_amplicons"] = offtargets
    cand["n_offtarget"] = len(offtargets)
    cand["min_3p_mismatches"] = (worst_discrimination
                                 if off_templates else TERM_WINDOW)
    return cand


# ==========================================
# PRIMERSEARCH - transcriptome- and genome-wide (optional)
# ==========================================
def primersearch_available(binary=None):
    """Resolve the binary at call time, not import time, so editing the
    module constant (or setting it from a wrapper) actually takes."""
    b = binary or PRIMERSEARCH_BIN
    return shutil.which(b) is not None or (
        os.path.isfile(b) and os.access(b, os.X_OK))


def write_primersearch_input(candidates, path):
    """EMBOSS -infile format: name <tab> forward <tab> reverse.

    Names must be whitespace-free; we use the candidate index so results
    map straight back.
    """
    with open(path, "w") as f:
        for i, c in enumerate(candidates):
            f.write(f"cand{i}\t{c['fwd']}\t{c['rev']}\n")
    return path


def parse_primersearch(path):
    """Parse primersearch -outfile into {name: [(seq_id, length), ...]}.

    Real EMBOSS output looks like:

        Primer name cand0
        Amplimer 1
        \tSequence: ENST00000456328.2
        \tHomo sapiens some description
        \tGCTGACT... hits forward strand at 100 with 0 mismatches
        \tCCAGTCA... hits reverse strand at [250] with 0 mismatches
        \tAmplimer length: 150 bp
    """
    results = defaultdict(list)
    name = None
    seq_id = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Primer name"):
                name = line.split("Primer name", 1)[1].strip()
                seq_id = None
            elif line.startswith("Amplimer ") and "length" not in line:
                seq_id = None
            elif line.startswith("Sequence:"):
                seq_id = line.split("Sequence:", 1)[1].strip()
                # EMBOSS sometimes prefixes a database tag: "EMBL:HSFAU"
                if ":" in seq_id:
                    seq_id = seq_id.split(":")[-1].strip()
                seq_id = seq_id.split()[0] if seq_id else seq_id
            elif line.startswith("Amplimer length:"):
                m = re.search(r"Amplimer length:\s*(\d+)", line)
                length = int(m.group(1)) if m else None
                if name is not None and seq_id is not None:
                    results[name].append((seq_id, length))
    return dict(results)


def run_primersearch(candidates, fasta_path, work_dir, tag,
                     mismatch_pct=None, binary=None):
    """Run one primersearch pass. Returns {cand_name: [(seq_id, len)]}.

    Returns None (rather than {}) when the run could not be performed,
    so callers can tell 'no hits' from 'not checked'.
    """
    binary = binary or PRIMERSEARCH_BIN
    mismatch_pct = (PRIMERSEARCH_MISMATCH_PCT if mismatch_pct is None
                    else mismatch_pct)
    if not fasta_path or not os.path.exists(fasta_path):
        return None
    if not primersearch_available(binary):
        return None

    infile = os.path.join(work_dir, f"primersearch_{tag}_in.txt")
    outfile = os.path.join(work_dir, f"primersearch_{tag}_out.txt")
    write_primersearch_input(candidates, infile)

    cmd = [binary,
           "-seqall", fasta_path,
           "-infile", infile,
           "-mismatchpercent", str(mismatch_pct),
           "-outfile", outfile,
           "-auto"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True,
                       timeout=3600)
    except subprocess.CalledProcessError as e:
        print(f"[!] primersearch ({tag}) failed: "
              f"{(e.stderr or '').strip()[:300]}")
        return None
    except subprocess.TimeoutExpired:
        print(f"[!] primersearch ({tag}) timed out")
        return None

    parsed = parse_primersearch(outfile)
    if not PRIMERSEARCH_KEEP_FILES:
        for p in (infile, outfile):
            try:
                os.remove(p)
            except OSError:
                pass
    return parsed


def annotate_with_primersearch(candidates, work_dir,
                               target_id, sibling_ids):
    """Attach transcriptome-wide and genome-wide PrimerSearch results.

    Transcriptome hits are split three ways:
      - the target itself          expected, ignored
      - a sibling isoform          already covered by the Python screen,
                                   recorded for cross-checking
      - anything else              NEW information: a paralog,
                                   pseudogene or unrelated transcript
    """
    def _strip(x):
        return x.split("|")[0].split(".")[0] if x else x

    target_base = _strip(target_id)
    sibling_base = {_strip(s) for s in sibling_ids}

    for c in candidates:
        c["ps_checked"] = False
        c["ps_offgene"] = []
        c["ps_sibling"] = []
        c["ps_genome_hits"] = None

    tx = run_primersearch(candidates, TRANSCRIPTOME_FASTA, work_dir,
                          "transcriptome")
    if tx is not None:
        for i, c in enumerate(candidates):
            c["ps_checked"] = True
            for seq_id, length in tx.get(f"cand{i}", []):
                base = _strip(seq_id)
                pretty = seq_id.split("|")[0]
                if base == target_base:
                    continue
                if base in sibling_base:
                    c["ps_sibling"].append((pretty, length))
                else:
                    c["ps_offgene"].append((pretty, length))

    gn = run_primersearch(candidates, GENOME_FASTA, work_dir, "genome")
    if gn is not None:
        for i, c in enumerate(candidates):
            c["ps_genome_hits"] = gn.get(f"cand{i}", [])

    return candidates


def crosscheck_primersearch(candidates):
    """Report where PrimerSearch and the Python screen disagree.

    They should mostly agree on sibling isoforms. Disagreement is
    informative, not an error: PrimerSearch uses a flat mismatch
    percentage, the Python screen requires a perfect 3'-terminal window.
    A PrimerSearch-only hit is usually a product that would not actually
    extend; a Python-only hit is usually a 3'-clean primer whose total
    mismatch count exceeds the PrimerSearch threshold.
    """
    ps_only = py_only = agree = 0
    for c in candidates:
        if not c.get("ps_checked"):
            continue
        ps = {s.split("|")[0].split(".")[0] for s, _ in c["ps_sibling"]}
        py = {o["isoform"] for o in c.get("offtarget_amplicons", [])}
        if ps and py:
            agree += 1
        elif ps and not py:
            ps_only += 1
        elif py and not ps:
            py_only += 1
    return {"agree": agree, "primersearch_only": ps_only,
            "python_only": py_only}


# ==========================================
# REPORTS
# ==========================================
CSS_COMMON = """
body { font-family: Menlo, Consolas, monospace; font-size: 13px;
       padding: 20px; color: #222; }
h1, h2 { color: #333; }
h2 { margin-top: 40px; }
.ruler { color: #aaa; font-size: 10px; margin-bottom: 6px;
         white-space: nowrap; overflow: hidden; }
.seq { background: #f8f8f8; padding: 16px 20px; line-height: 2.0;
       word-break: break-all; border-radius: 4px; font-size: 13px;
       letter-spacing: 0.5px; }
.base_difference { color: #1a1a1a; font-weight: bold;
                   background: #bfbfbf; padding: 3px 2px;
                   border-radius: 2px;
                   border-bottom: 3px dotted #555; }
.chimeric_region { color: #1a1a1a; font-weight: bold;
                   background: #e8e8e8; padding: 3px 2px;
                   border-radius: 2px;
                   border-bottom: 3px solid #999; }
.fwd { background: #c8e6ff; color: #003a70; font-weight: bold;
       padding: 3px 2px; border-bottom: 3px solid #0055aa; }
.rev { background: #d6f5d6; color: #0a4d0a; font-weight: bold;
       padding: 3px 2px; border-bottom: 3px solid #228b22; }
.candidate { margin-bottom: 32px; padding: 16px; background: #fafafa;
             border-radius: 6px; border-left: 4px solid #4a90e2; }
.specific { border-left-color: #2e7d32; }
.nonspecific { border-left-color: #c62828; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 10px;
         font-size: 11px; font-weight: bold; margin-left: 8px; }
.badge-ok { background: #c8e6c9; color: #1b5e20; }
.badge-bad { background: #ffcdd2; color: #b71c1c; }
table { border-collapse: collapse; margin-top: 12px;
        font-size: 12px; }
table td, table th { padding: 6px 12px; border: 1px solid #ccc;
                     text-align: left; }
table th { background: #f0f0f0; }
table tr.summary td { background: #f0f4f9; font-weight: 600;
                      border-top: 2px solid #b0bec5; }
table tr.summary td b { color: #111; }
code { background: #f0f0f0; padding: 2px 6px; border-radius: 3px; }
.track-label { color: #666; font-size: 11px; margin-bottom: 2px;
               margin-top: 8px; }
"""


def _build_ruler(n):
    out = ["<div class='ruler'>"]
    for i in range(0, n, 50):
        out.append(f"|{i}" + "&nbsp;" * 41)
    out.append("</div>")
    return "".join(out)


def _render_span_sequence(target_cdna, classes):
    out = []
    i = 0
    n = len(target_cdna)
    while i < n:
        cur = classes[i]
        j = i
        while j < n and classes[j] == cur:
            j += 1
        chunk = target_cdna[i:j]
        if cur:
            out.append(f"<span class='{cur}'>{chunk}</span>")
        else:
            out.append(f"<span style='color:#888;'>{chunk}</span>")
        i = j
    return "".join(out)

def write_unified_report(out_dir, display_name, target_base_id,
                         target_cdna, regions, aligned, other_ids):
    safe = safe_filename(display_name)

    tsv_path = os.path.join(out_dir, f"{safe}_regions.tsv")
    with open(tsv_path, "w") as f:
        f.write(f"# {display_name} ({target_base_id})\n")
        f.write("region_start\tregion_end\tlength\tdefinition\t"
                "upstream_50bp\tregion_sequence\tdownstream_50bp\n")
        for r in regions:
            s, e = r["start"], r["end"]
            f.write(f"{s}\t{e}\t{r['length']}\t{r['definition']}\t"
                    f"{target_cdna[max(0, s-CONTEXT_FLANK):s]}\t"
                    f"{target_cdna[s:e]}\t"
                    f"{target_cdna[e:e+CONTEXT_FLANK]}\n")

    n = len(target_cdna)
    classes = [""] * n
    for r in regions:
        for i in range(r["start"], r["end"]):
            classes[i] = r["definition"]

    html_path = os.path.join(out_dir, f"{safe}_regions.html")
    with open(html_path, "w") as f:
        f.write("<!DOCTYPE html><html><head><meta charset='utf-8'>")
        f.write(f"<title>Unique regions in {display_name}</title>")
        f.write("<style>" + CSS_COMMON + "</style></head><body>")
        f.write(f"<h1>{display_name} "
                f"<span style='color:#888; font-weight:normal;'>"
                f"({target_base_id})</span></h1>")
        f.write(f"<p>Total length: {n} nt. "
                f"{len(regions)} region(s) found.</p>")
        f.write("<p><b>Legend:</b> "
                "<span class='base_difference'>&nbsp;base difference"
                "&nbsp;</span> "
                "(dotted underline, dark grey) &mdash; no off-target shares "
                "the aligned base. "
                "<span class='chimeric_region'>&nbsp;chimeric region"
                "&nbsp;</span> "
                "(solid underline, light grey) &mdash; no off-target "
                f"contains this {CHIMERIC_WINDOW} bp string verbatim."
                "</p>")
        f.write(_build_ruler(n))
        f.write("<div class='seq'>")
        f.write(_render_span_sequence(target_cdna, classes))
        f.write("</div>")
        f.write("<h2>Region list</h2>")
        f.write("<table>")
        f.write("<tr><th>Start</th><th>End</th><th>Length</th>"
                "<th>Definition</th><th>Sequence</th></tr>")
        for r in regions:
            f.write(f"<tr><td>{r['start']}</td><td>{r['end']}</td>"
                    f"<td>{r['length']}</td>"
                    f"<td><span class='{r['definition']}'>"
                    f"{r['definition']}</span></td>"
                    f"<td><code>{target_cdna[r['start']:r['end']]}"
                    f"</code></td></tr>")
        f.write("</table>")
        f.write("</body></html>")

    txt_path = os.path.join(out_dir, f"{safe}_regions.txt")
    with open(txt_path, "w") as f:
        f.write(f"{display_name} ({target_base_id})\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Total length: {n} nt. "
                f"{len(regions)} region(s) found.\n\n")
        for r in regions:
            s, e = r["start"], r["end"]
            f.write(f"{r['definition']:20s} "
                    f"cDNA {s}-{e} ({r['length']} bp)\n")
            f.write(f"  sequence: {target_cdna[s:e]}\n\n")

    pos_path = os.path.join(out_dir, f"{safe}_positions.csv")
    with open(pos_path, "w", newline="") as f:
        w = csv.writer(f)
        header = (["position", "base", "definition"] +
                  [f"ident_vs_{nm}" for nm in other_ids] +
                  [f"base_vs_{nm}" for nm in other_ids])
        w.writerow(header)
        for p in range(n):
            row = [p, target_cdna[p], classes[p] or ""]
            for idx in range(len(other_ids)):
                ident, _ = aligned[p].get(idx, (None, None))
                row.append(f"{ident:.3f}" if ident is not None else "")
            for idx in range(len(other_ids)):
                _, base = aligned[p].get(idx, (None, None))
                row.append(base if base else "")
            w.writerow(row)

    return tsv_path, html_path, txt_path, pos_path


def write_assay_reports(out_dir, display_name, target_base_id,
                        gene_name, target_cdna, candidates, regions,
                        junctions):
    if not candidates:
        return None, None
    safe = safe_filename(display_name)

    csv_path = os.path.join(out_dir, f"{safe}_assay.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Gene", "Transcript_ID", "Isoform_Name", "Rank",
            "Forward_Primer_Sequence_5to3",
            "Reverse_Primer_Sequence_5to3",
            "Forward_Primer_Length", "Reverse_Primer_Length",
            "Amplicon_Size_bp",
            "Forward_Primer_Position", "Reverse_Primer_Position",
            # --- specificity ---
            "Anchor_Class", "Anchor_Orientation", "Anchor_Position",
            "Offtarget_Isoforms_Amplified",
            "Offtarget_Detail",
            "Min_3p_Mismatches_vs_Any_Isoform",
            "Spans_Exon_Junction",
            # --- EMBOSS PrimerSearch (blank when not run) ---
            "PrimerSearch_Run", "PS_Offgene_Amplicons",
            "PS_Offgene_Detail", "PS_Sibling_Amplicons",
            "PS_Genomic_Amplicons",
            # --- thermodynamics ---
            "Tm_F_C", "Tm_R_C", "Tm_Difference_C",
            "GC_F_percent", "GC_R_percent",
            "HOMODIMER_DG_F_kcal_mol", "HAIRPIN_DG_F_kcal_mol",
            "HOMODIMER_DG_R_kcal_mol", "HAIRPIN_DG_R_kcal_mol",
            "HETERODIMER_DG_kcal_mol", "END_STABILITY_DG_kcal_mol",
        ])
        for rank, c in enumerate(candidates, 1):
            tm_diff = round(abs(c["fwd_tm"] - c["rev_tm"]), 2)
            detail = "; ".join(
                f"{o['isoform']}:{o['product_size']}bp"
                f"(mmF{o['fwd_total_mm']}/mmR{o['rev_total_mm']})"
                for o in c["offtarget_amplicons"]
            ) or "none"
            w.writerow([
                gene_name, target_base_id, display_name, rank,
                c["fwd"], c["rev"], c["fwd_len"], c["rev_len"],
                c["product_size"],
                f"{c['fwd_start']}-{c['fwd_end']}",
                f"{c['rev_start']}-{c['rev_end']}",
                c["anchor_class"], c["orientation"], c["anchor_position"],
                c["n_offtarget"], detail,
                c["min_3p_mismatches"],
                "yes" if c.get("spans_junction") else "no",
                "yes" if c.get("ps_checked") else "no",
                len(c.get("ps_offgene") or ()),
                "; ".join(f"{sid}:{ln}bp"
                          for sid, ln in (c.get("ps_offgene") or ()))
                or "none",
                len(c.get("ps_sibling") or ()),
                ("n/a" if c.get("ps_genome_hits") is None
                 else len(c["ps_genome_hits"])),
                c["fwd_tm"], c["rev_tm"], tm_diff,
                c["fwd_gc"], c["rev_gc"],
                c["homo_f"], c["hairpin_f"],
                c["homo_r"], c["hairpin_r"],
                c["hetero"], c["end_stab"],
            ])

    html_path = os.path.join(out_dir, f"{safe}_assay.html")
    n = len(target_cdna)

    region_classes = [""] * n
    for r in regions:
        for i in range(r["start"], r["end"]):
            region_classes[i] = r["definition"]

    with open(html_path, "w") as f:
        f.write("<!DOCTYPE html><html><head><meta charset='utf-8'>")
        f.write(f"<title>Assay designs for {display_name}</title>")
        f.write("<style>" + CSS_COMMON + "</style></head><body>")
        f.write(f"<h1>{display_name} "
                f"<span style='color:#888; font-weight:normal;'>"
                f"({target_base_id})</span> &mdash; primer pair designs</h1>")
        f.write(f"<p>Total length: {n} nt. "
                f"{len(candidates)} candidate(s). Ranked by predicted "
                f"isoform specificity (in-silico PCR against sibling "
                f"isoforms), then 3'-terminal discrimination margin, "
                f"then thermodynamics.</p>")
        f.write(f"<p><b>Specificity model:</b> a primer is assumed to "
                f"extend on an off-target isoform when its 3'-terminal "
                f"{TERM_WINDOW} bases match perfectly and total "
                f"mismatches are &le; {MAX_TOTAL_MM}. This is a "
                f"within-gene check and it is 3'-end aware. EMBOSS "
                f"PrimerSearch, when reference FASTA files are supplied, "
                f"adds the complement &mdash; whole transcriptome and "
                f"genome at a flat "
                f"{PRIMERSEARCH_MISMATCH_PCT}% mismatch allowance, not "
                f"3'-end aware &mdash; so paralogs, pseudogenes and gDNA "
                f"products are caught too.</p>")
        f.write("<p><b>Legend:</b> "
                "<span class='base_difference'>&nbsp;base difference"
                "&nbsp;</span> "
                "<span class='chimeric_region'>&nbsp;chimeric region"
                "&nbsp;</span> "
                "<span class='fwd'>&nbsp;forward primer&nbsp;</span> "
                "<span class='rev'>&nbsp;reverse primer&nbsp;</span></p>")

        for idx, c in enumerate(candidates, 1):
            ok = (c["n_offtarget"] == 0
                  and not c.get("ps_offgene"))
            f.write(f"<div class='candidate "
                    f"{'specific' if ok else 'nonspecific'}'>")
            if ok:
                badge = "<span class='badge badge-ok'>specific</span>"
            else:
                bits = []
                if c["n_offtarget"]:
                    bits.append(f"{c['n_offtarget']} sibling isoform(s)")
                if c.get("ps_offgene"):
                    bits.append(f"{len(c['ps_offgene'])} off-gene "
                                f"transcript(s)")
                badge = (f"<span class='badge badge-bad'>amplifies "
                         f"{' + '.join(bits)}</span>")
            f.write(f"<h2>Rank {idx}{badge}</h2>")

            f.write("<div class='track-label'>"
                    "Reference: unique regions</div>")
            f.write(_build_ruler(n))
            f.write("<div class='seq'>")
            f.write(_render_span_sequence(target_cdna, region_classes))
            f.write("</div>")

            cls = [""] * n
            for i in range(c["fwd_start"], c["fwd_end"]):
                cls[i] = "fwd"
            for i in range(c["rev_start"], c["rev_end"]):
                cls[i] = "rev"

            f.write("<div class='track-label'>"
                    "Assay: forward / reverse</div>")
            f.write(_build_ruler(n))
            f.write("<div class='seq'>")
            f.write(_render_span_sequence(target_cdna, cls))
            f.write("</div>")

            tm_diff = round(abs(c["fwd_tm"] - c["rev_tm"]), 2)
            f.write("<table>")
            f.write(
                "<tr>"
                "<th>Oligo</th><th>Sequence (5'&rarr;3')</th>"
                "<th>Position</th><th>Length</th><th>Tm (&deg;C)</th>"
                "<th>GC%</th><th>Homodimer &Delta;G</th>"
                "<th>Hairpin &Delta;G</th><th>Heterodimer &Delta;G</th>"
                "<th>3' end stability</th></tr>"
            )
            for label, seq, st, en, ln, t_, g_, ho, ha in (
                ("Forward primer", c["fwd"], c["fwd_start"], c["fwd_end"],
                 c["fwd_len"], c["fwd_tm"], c["fwd_gc"],
                 c["homo_f"], c["hairpin_f"]),
                ("Reverse primer", c["rev"], c["rev_start"], c["rev_end"],
                 c["rev_len"], c["rev_tm"], c["rev_gc"],
                 c["homo_r"], c["hairpin_r"]),
            ):
                f.write(
                    f"<tr><td>{label}</td><td><code>{seq}</code></td>"
                    f"<td>{st}&ndash;{en}</td><td>{ln}</td><td>{t_}</td>"
                    f"<td>{g_}</td><td>{ho}</td><td>{ha}</td>"
                    f"<td>{c['hetero']}</td><td>{c['end_stab']}</td></tr>"
                )
            f.write(
                f"<tr class='summary'><td colspan='10'>"
                f"<b>Amplicon:</b> {c['product_size']} bp &nbsp;|&nbsp; "
                f"<b>Tm difference:</b> {tm_diff} &deg;C "
                f"(limit {TM_BALANCE_MAX}) &nbsp;|&nbsp; "
                f"<b>Anchor:</b> {c['anchor_class']} at position "
                f"{c['anchor_position']} ({c['orientation']}) "
                f"&nbsp;|&nbsp; <b>Spans exon junction:</b> "
                f"{'yes' if c.get('spans_junction') else 'no'}</td></tr>"
            )
            if c["offtarget_amplicons"]:
                detail = ", ".join(
                    f"{o['isoform']} ({o['product_size']} bp)"
                    for o in c["offtarget_amplicons"])
                f.write(
                    f"<tr class='summary'><td colspan='10'>"
                    f"<b>Predicted off-target amplicons:</b> {detail}"
                    f"</td></tr>")
            else:
                f.write(
                    f"<tr class='summary'><td colspan='10'>"
                    f"<b>No sibling isoform amplifies.</b> Smallest "
                    f"3'-terminal mismatch margin against any isoform: "
                    f"{c['min_3p_mismatches']} of {TERM_WINDOW} bases."
                    f"</td></tr>")
            if c.get("ps_checked"):
                og = c.get("ps_offgene") or []
                gh = c.get("ps_genome_hits")
                og_txt = (", ".join(f"{sid} ({ln} bp)" for sid, ln in og)
                          if og else "none")
                if gh is None:
                    gh_txt = "not checked"
                elif gh:
                    short = [ln for _, ln in gh if ln and ln < 1000]
                    gh_txt = (f"{len(gh)} reported, shortest "
                              f"{min(ln for _, ln in gh if ln)} bp"
                              + (" - CHECK, a short product is a real "
                                 "gDNA risk" if short else
                                 " - all far too long to amplify"))
                else:
                    gh_txt = "none"
                f.write(
                    f"<tr class='summary'><td colspan='10'>"
                    f"<b>PrimerSearch, transcriptome-wide:</b> off-gene "
                    f"amplicons &mdash; {og_txt} &nbsp;|&nbsp; "
                    f"<b>genomic:</b> {gh_txt}</td></tr>")
            else:
                f.write(
                    "<tr class='summary'><td colspan='10'>"
                    "<b>Within-gene check only</b> &mdash; PrimerSearch "
                    "was not run. Confirm the chosen pair with NCBI "
                    "Primer-BLAST before ordering.</td></tr>")
            f.write(
                "<tr class='summary'><td colspan='10'>&Delta;G values "
                "are primer3-py thermodynamic calculations. More "
                "negative = stronger structure = worse.</td></tr>")
            f.write("</table></div>")

        f.write("</body></html>")

    return csv_path, html_path


# ==========================================
# PIPELINE
# ==========================================
    
log = logging.getLogger("isosolo")


class TargetResult:
    """Outcome of one target. ``status`` is the thing to branch on."""

    OK = "ok"
    NO_UNIQUE_REGION = "no_unique_region"   # target is a sequence subset
    NO_CANDIDATES = "no_candidates"         # unique region, bad chemistry
    NO_SEQUENCE = "no_sequence"             # fetch failed

    def __init__(self, gene, transcript_id, status, **kw):
        self.gene = gene
        self.transcript_id = transcript_id
        self.status = status
        self.display_name = kw.get("display_name", transcript_id)
        self.candidates = kw.get("candidates", [])
        self.intervals = kw.get("intervals", [])
        self.n_siblings = kw.get("n_siblings", 0)
        self.outputs = kw.get("outputs", {})
        self.message = kw.get("message", "")

    @property
    def best(self):
        return self.candidates[0] if self.candidates else None

    def summary_row(self):
        c = self.best
        return {
            "gene": self.gene,
            "transcript": self.transcript_id,
            "isoform": self.display_name,
            "status": self.status,
            "n_siblings": self.n_siblings,
            "n_unique_intervals": len(self.intervals),
            "unique_bp": sum(e - s for s, e, _ in self.intervals),
            "anchor_classes": ",".join(sorted({c for _, _, c
                                               in self.intervals})) or "",
            "n_candidates": len(self.candidates),
            "best_forward": c["fwd"] if c else "",
            "best_reverse": c["rev"] if c else "",
            "best_amplicon_bp": c["product_size"] if c else "",
            "best_anchor_class": c["anchor_class"] if c else "",
            "sibling_amplicons": c["n_offtarget"] if c else "",
            "offgene_amplicons": (len(c.get("ps_offgene") or ())
                                  if c else ""),
            "spans_junction": (("yes" if c.get("spans_junction") else "no")
                               if c else ""),
            "message": self.message,
        }


def run_target(gene, transcript_id, output_root):
    """Design isoform-specific primers for one transcript."""
    target_base_id = transcript_id.split(".")[0]
    log.info("=" * 70)
    log.info("%s | %s", gene, transcript_id)

    gene_data = fetch_gene_transcripts(gene)
    if not gene_data:
        return TargetResult(gene, transcript_id, TargetResult.NO_SEQUENCE,
                            message=f"gene {gene} not found")

    names = gene_data["names"]
    display_name = names.get(target_base_id, target_base_id)
    iso_dir = ensure_dir(os.path.join(output_root,
                                      safe_filename(display_name)))

    target_cdna = fetch_cdna(target_base_id)
    if not target_cdna:
        return TargetResult(gene, transcript_id, TargetResult.NO_SEQUENCE,
                            display_name=display_name,
                            message="no cDNA for target")
    target_cdna = target_cdna.upper()
    log.info("target %s (%s), %d nt", display_name, target_base_id,
             len(target_cdna))

    junctions = fetch_exon_junctions(target_base_id)
    log.info("%d exon-exon junction(s)", len(junctions))

    other_cdnas, other_ids = [], []
    for t in gene_data["transcripts"]:
        tid = t["id"].split(".")[0]
        if tid == target_base_id:
            continue
        s = fetch_cdna(tid)
        if s:
            other_cdnas.append(s.upper())
            other_ids.append(names.get(tid, tid))
    log.info("%d sibling isoform(s)", len(other_cdnas))

    # ---- uniqueness -------------------------------------------------
    t0 = time.time()
    chim = chimeric_region_flags(target_cdna, other_cdnas)
    base, aligned = base_difference_flags(target_cdna, other_cdnas)
    log.info("uniqueness detection %.1fs", time.time() - t0)

    regions = merge_runs(chim, base)
    intervals = merge_into_uniqueness_intervals(regions)

    if not regions:
        return TargetResult(
            gene, transcript_id, TargetResult.NO_UNIQUE_REGION,
            display_name=display_name, n_siblings=len(other_cdnas),
            message=("no sequence unique to this transcript - it is a "
                     "subset of a sibling, so no primer pair can "
                     "distinguish them"))

    log.info("%d region(s) -> %d interval(s), %d bp unique",
             len(regions), len(intervals),
             sum(e - s for s, e, _ in intervals))

    outputs = {}
    paths = write_unified_report(iso_dir, display_name, target_base_id,
                                 target_cdna, regions, aligned, other_ids)
    outputs["regions_tsv"], outputs["regions_html"] = paths[0], paths[1]
    outputs["regions_txt"], outputs["positions_csv"] = paths[2], paths[3]

    # ---- design -----------------------------------------------------
    t0 = time.time()
    candidates = design_assay(target_cdna, intervals)
    log.info("%d thermodynamically valid candidate(s) in %.1fs",
             len(candidates), time.time() - t0)
    if not candidates:
        return TargetResult(
            gene, transcript_id, TargetResult.NO_CANDIDATES,
            display_name=display_name, intervals=intervals,
            n_siblings=len(other_cdnas), outputs=outputs,
            message=("unique sequence exists but no primer pair met the "
                     "Tm/GC/amplicon constraints - try widening them"))

    # ---- within-gene specificity ------------------------------------
    candidates.sort(key=thermo_key)
    screen = candidates[:SPECIFICITY_SCREEN_TOP]
    off_indices = [build_kmer_index(o, SEED_LEN)
                   for o in other_cdnas]
    cache = {}
    for c in screen:
        screen_specificity(c, target_cdna, other_cdnas, other_ids,
                           off_indices, cache)
    ranked = rank_candidates(screen, junctions)

    # ---- optional: transcriptome / genome specificity ---------------
    shortlist = ranked[:max(TOP_N, PRIMERSEARCH_SCREEN_TOP)]
    if TRANSCRIPTOME_FASTA or GENOME_FASTA:
        if not primersearch_available():
            log.warning("  '%s' not on PATH - skipping the "
                        "transcriptome/genome check "
                        "(conda install -c bioconda emboss)",
                        PRIMERSEARCH_BIN)
        else:
            shortlist = shortlist[:PRIMERSEARCH_SCREEN_TOP]
            log.info("  PrimerSearch on the top %d pair(s) at %d%% "
                     "mismatch - about a minute each...",
                     len(shortlist), PRIMERSEARCH_MISMATCH_PCT)
            t0 = time.time()
            annotate_with_primersearch(shortlist, iso_dir,
                                       target_base_id, other_ids)
            log.info("  PrimerSearch %.0fs; %d/%d hit a transcript "
                     "outside %s; cross-check vs the sibling screen: %s",
                     time.time() - t0,
                     sum(1 for c in shortlist if c["ps_offgene"]),
                     len(shortlist), gene,
                     crosscheck_primersearch(shortlist))
            shortlist = rank_candidates(shortlist, junctions)

    final = shortlist[:TOP_N]
    n_clean = sum(1 for c in final
                  if c["n_offtarget"] == 0 and not c.get("ps_offgene"))
    log.info("%d/%d reported pairs amplify nothing else", n_clean,
             len(final))

    acsv, ahtml = write_assay_reports(iso_dir, display_name,
                                      target_base_id, gene, target_cdna,
                                      final, regions, junctions)
    outputs["assay_csv"], outputs["assay_html"] = acsv, ahtml

    return TargetResult(gene, transcript_id, TargetResult.OK,
                        display_name=display_name, candidates=final,
                        intervals=intervals, n_siblings=len(other_cdnas),
                        outputs=outputs)

# ==========================================
# MAIN
# ==========================================
if os.path.exists(OUTPUT_ROOT):
    shutil.rmtree(OUTPUT_ROOT)
ensure_dir(OUTPUT_ROOT)
print(f"[*] Output root: {OUTPUT_ROOT}")

_results = []
for _job in TARGET_JOBS:
    print()
    try:
        _r = run_target(_job["gene"], _job["transcript"], OUTPUT_ROOT)
    except Exception as _exc:
        log.exception("%s | %s failed", _job["gene"], _job["transcript"])
        continue
    _results.append(_r)

    if _r.status != TargetResult.OK:
        print(f"[!] {_r.status}: {_r.message}")
        continue

    _c = _r.best
    _d = round(abs(_c["fwd_tm"] - _c["rev_tm"]), 2)
    print(f"\n    Best pair:")
    print(f"      Fwd      : {_c['fwd']}  "
          f"(Tm {_c['fwd_tm']}, GC {_c['fwd_gc']}%)")
    print(f"      Rev      : {_c['rev']}  "
          f"(Tm {_c['rev_tm']}, GC {_c['rev_gc']}%)")
    print(f"      Amplicon : {_c['product_size']} bp  |  Tm diff {_d} C")
    print(f"      Anchor   : {_c['anchor_class']} "
          f"({_c['orientation']}) at cDNA {_c['anchor_position']}")
    print(f"      Siblings amplified : {_c['n_offtarget']}")
    if _c.get("ps_checked"):
        print(f"      Off-gene amplicons : "
              f"{len(_c.get('ps_offgene') or ())}")
    print(f"      3' mismatch margin : {_c['min_3p_mismatches']}"
          f"/{TERM_WINDOW}")
    print(f"      Spans exon junction: "
          f"{'yes' if _c.get('spans_junction') else 'no'}")
    for _k in ("assay_csv", "assay_html", "regions_html"):
        if _r.outputs.get(_k):
            print(f"      {_r.outputs[_k]}")

# cross-target summary, only when there is more than one
if len(_results) > 1:
    _sum = os.path.join(OUTPUT_ROOT, "summary.csv")
    _rows = [r.summary_row() for r in _results]
    with open(_sum, "w", newline="") as _f:
        _w = csv.DictWriter(_f, fieldnames=list(_rows[0]))
        _w.writeheader()
        _w.writerows(_rows)
    print(f"\n[+] Summary: {_sum}")

print("\n" + "=" * 70)
print(f"[done] Output: {OUTPUT_ROOT}")
print("Open the *_assay.html file in a browser first.")
print("Then run the chosen pair through NCBI Primer-BLAST.")
print("=" * 70)
