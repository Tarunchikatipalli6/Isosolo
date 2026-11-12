#!/usr/bin/env python3
"""
check_primers.py - which transcripts does this primer pair amplify?
====================================================================

Give it a gene and a forward/reverse primer pair. It fetches every
transcript of that gene from Ensembl and reports, for each one:

    AMPLIFIES   both primers can prime, facing each other, sane product
    no          at least one primer cannot prime there

For every transcript it also shows WHERE each primer binds best and
how many mismatches sit there - counting the 3'-terminal ones
separately, since those are what stop a primer extending.

    pip install requests

Edit the block below, then run.

WHAT IT IS FOR
    Independent confirmation of a pair isosolo chose, or a check on a
    pair from anywhere else. It searches every position by brute
    force, where isosolo uses k-mer seeding for speed - so when the
    two disagree, this one is the more thorough.

    One gene only. It cannot see other genes; use NCBI Primer-BLAST
    for that.

THE RULE
    Matches NCBI Primer-BLAST defaults: a site primes unless it
    carries 2+ mismatches in the last 5 bases or 6+ in total.

READING THE OUTPUT
    pos 363, 0 mm (0 at 3')   binds at base 363, no mismatches, none
                              of them in the last 5 bases
    no site                   no plausible binding position at all -
                              the strongest possible result

    SPECIFIC                  one isoform amplifies. What you want.
    NOT SPECIFIC              names the others it also hits.
    No transcript amplifies   usually a primer typed on the wrong
                              strand. Both go in 5' to 3' exactly as
                              ordered; the reverse is complemented
                              here, not by you.

Genes with many isoforms are slow - TFRC has 61 and takes a few
minutes. On a cluster, run it on a compute node; login nodes kill it.
"""

# =====================================================================
# ===  EDIT THIS, THEN PRESS RUN  =====================================
# =====================================================================

GENE = "HNRNPC"

FORWARD = "CAAACAAGCAGCTGGAGTTG"
REVERSE = "AGTCATCCTCGCCATTGGC"

# A primer is assumed to extend on a template only if it carries FEWER
# THAN TERM_MM_TO_BLOCK mismatches within the last TERM_WINDOW bases at
# the 3' end, and no more than MAX_TOTAL_MM mismatches overall.
#
# The 5 nt window: Stadhouders et al. 2010 J Mol Diagn 12:109-117;
#   same window used by Primer-BLAST (Ye et al. 2012 BMC Bioinf 13:134).
# Two mismatches to block: Kwok et al. 1990 Nucleic Acids Res 18:999-1005
#   found a single 3'-terminal mismatch is frequently tolerated while two
#   within the last four bases reduce yield dramatically.
# Five total: NCBI Primer-BLAST ignores a site only at 6+ total
#   mismatches (Ye et al. 2012, default "mismatch threshold to ignore
#   targets" = 6). Christopherson et al. 1997 Nucleic Acids Res
#   25:654-658 found 2-4 mismatches had little effect on yield and 5
#   reduced it ~22-fold - reduced, not abolished. Set to 5 so this
#   tool and Primer-BLAST agree on what counts as an off-target;
#   at 3 they disagreed on a real MET-204 site.
TERM_WINDOW = 5
TERM_MM_TO_BLOCK = 2
MAX_TOTAL_MM = 5
MAX_PRODUCT = 2000        # bp; anything longer is not a PCR product

# =====================================================================

import json
import os
import sys

try:
    import requests
except ImportError:
    sys.exit("needs requests:\n    pip install requests")

CACHE_DIR = os.path.expanduser("~/.ensembl_cache")
os.makedirs(CACHE_DIR, exist_ok=True)
HEADERS = {"Content-Type": "application/json"}


def revcomp(s):
    c = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(c.get(b, "N") for b in reversed(s.upper()))


def _cache(key):
    safe = "".join(ch if ch.isalnum() else "_" for ch in key)
    return os.path.join(CACHE_DIR, f"{safe}.json")


def _get(url, key):
    p = _cache(key)
    if os.path.exists(p):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            pass
    r = requests.get(url, headers=HEADERS, timeout=30)
    if not r.ok:
        return None
    d = r.json()
    try:
        with open(p, "w") as f:
            json.dump(d, f)
    except Exception:
        pass
    return d


def transcripts_of(gene):
    g = _get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/"
             f"{gene}?content-type=application/json", f"symbol_{gene}")
    if not g or "id" not in g:
        sys.exit(f"gene {gene} not found")
    e = _get(f"https://rest.ensembl.org/lookup/id/{g['id']}?expand=1;"
             f"content-type=application/json", f"gene_expand_{g['id']}")
    out = []
    for t in e.get("Transcript", []):
        tid = t["id"].split(".")[0]
        out.append((tid, t.get("display_name") or tid))
    return out


def cdna(tid):
    p = _cache(f"cdna_{tid}.txt")
    if os.path.exists(p):
        s = open(p).read().strip()
        if s:
            return s.strip('"')
    d = _get(f"https://rest.ensembl.org/sequence/id/{tid}"
             f"?content-type=application/json;type=cdna", f"cdna_{tid}")
    s = (d or {}).get("seq", "")
    if s:
        open(p, "w").write(s)
    return s


def best_site(footprint, three_prime_right, template):
    """Where would this primer bind most easily on this template?

    Returns (terminal_mismatches, total_mismatches, position) for the
    site that is easiest to extend from, or None.
    """
    L, best = len(footprint), None
    for j in range(len(template) - L + 1):
        sub = template[j:j + L]
        total = sum(1 for a, b in zip(footprint, sub) if a != b)
        if total > 8:                      # nowhere near; skip
            continue
        if three_prime_right:
            term = sum(1 for a, b in zip(footprint[-TERM_WINDOW:],
                                         sub[-TERM_WINDOW:]) if a != b)
        else:
            term = sum(1 for a, b in zip(footprint[:TERM_WINDOW],
                                         sub[:TERM_WINDOW]) if a != b)
        # "Best" means most likely to EXTEND, not merely the cleanest
        # 3' end. A site with term=0 but 9 total mismatches cannot
        # prime; one with term=1 and 5 total can. Sorting on term
        # alone picked the former, rejected it, and hid the latter.
        primes = term < TERM_MM_TO_BLOCK and total <= MAX_TOTAL_MM
        key = (0 if primes else 1, term, total)
        cur = ((0 if (best[0] < TERM_MM_TO_BLOCK
                      and best[1] <= MAX_TOTAL_MM) else 1),
               best[0], best[1]) if best else None
        if best is None or key < cur:
            best = (term, total, j)
    return best


def main():
    fwd = FORWARD.upper().strip()
    rev = REVERSE.upper().strip()
    f_fp = fwd                       # forward footprint, sense strand
    r_fp = revcomp(rev)              # reverse footprint, sense strand

    print(f"Gene      : {GENE}")
    print(f"Forward   : {fwd}  ({len(fwd)} nt)")
    print(f"Reverse   : {rev}  ({len(rev)} nt)")
    print(f"Rule      : extends if fewer than {TERM_MM_TO_BLOCK} "
          f"mismatches in the last {TERM_WINDOW} bases and total "
          f"mismatches <= {MAX_TOTAL_MM}")
    print()

    txs = transcripts_of(GENE)
    print(f"{len(txs)} transcript(s) in Ensembl\n")

    hits = []
    rows = []
    for tid, name in txs:
        seq = (cdna(tid) or "").upper()
        if not seq:
            rows.append((name, tid, len(seq), "no sequence", "", ""))
            continue

        fb = best_site(f_fp, True, seq)
        rb = best_site(r_fp, False, seq)

        def fmt(b):
            if not b:
                return "no site"
            return f"pos {b[2]}, {b[1]} mm ({b[0]} at 3')"

        amplifies = False
        product = ""
        if fb and rb:
            if (fb[0] < TERM_MM_TO_BLOCK and rb[0] < TERM_MM_TO_BLOCK
                    and fb[1] <= MAX_TOTAL_MM and rb[1] <= MAX_TOTAL_MM):
                size = (rb[2] + len(r_fp)) - fb[2]
                if 0 < size <= MAX_PRODUCT:
                    amplifies = True
                    product = f"{size} bp"
        if amplifies:
            hits.append((name, product))
        rows.append((name, tid, len(seq),
                     "AMPLIFIES" if amplifies else "no",
                     fmt(fb), fmt(rb)))
        if amplifies:
            rows[-1] = rows[-1][:3] + (f"AMPLIFIES  {product}",
                                       fmt(fb), fmt(rb))

    w = max(len(r[0]) for r in rows) + 2
    print(f"{'isoform':<{w}}{'length':>7}   {'result':<22}"
          f"{'forward primer':<30}{'reverse primer'}")
    print("-" * (w + 7 + 3 + 22 + 30 + 28))
    for name, tid, ln, res, f_, r_ in rows:
        print(f"{name:<{w}}{ln:>7}   {res:<22}{f_:<30}{r_}")

    print()
    if not hits:
        print("No transcript of this gene amplifies. Check the primer "
              "sequences - one of them may be from the wrong strand.")
    elif len(hits) == 1:
        print(f"SPECIFIC: only {hits[0][0]} amplifies ({hits[0][1]}).")
    else:
        print(f"NOT SPECIFIC: {len(hits)} isoforms amplify -")
        for n, p in hits:
            print(f"    {n}  {p}")

    print("\nmm = mismatches. The count in brackets is how many of them "
          f"fall in the last {TERM_WINDOW} bases at the 3' end, which is"
          "\nwhere polymerase starts - a mismatch there usually stops "
          "the primer working at all.")
    print("\nThis covers one gene only. For other genes, use NCBI "
          "Primer-BLAST.")


if __name__ == "__main__":
    main()
