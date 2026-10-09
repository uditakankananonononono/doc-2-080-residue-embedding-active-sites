# DOC-2-080 "Functional motifs without alignment": does a per-residue ESM-2 representation locate annotated enzyme active sites in held-out Pfam families better than local amino-acid context? (frozen protocol, lock-1)

Written 2026-10-09 IST and committed BEFORE the active-site data was downloaded, sampled or embedded. The ledger gives only the title and source line for DOC-2-080 (no spec text), so this is one narrow, checkable test of "motifs without alignment".

## Why this is not a re-ask
- DOC-2-073/077/009-R3/R4: zero-shot variant-effect rho vs features. DOC-2-072: family-split EC classification. DOC-2-078: agreement as correctness flag. DOC-2-079: dark vs characterized 1-NN family recovery. None predicts residue-level function.
- Prior art: protein language models are known to carry residue-level functional signal (catalytic and binding sites). No novelty claimed. The question is whether this holds for a 35M model on family-held-out data, and by how much over a strong local-context baseline.

## Data (frozen)
- UniProtKB reviewed, length 100-400, single Pfam id, annotated active site (ft_act_site; 1 to 5 active-site residues), no X/B/Z/U/O residues (acquire_act.py; release and md5 in DATA_HASHES.tsv).
- Keep families with >= 6 such proteins; up to 200 families (seed 80); exactly 6 proteins per family (build_080.py). Expect up to 1,200 proteins.
- Residues: all annotated active-site residues (positives) plus 15 random non-active residues per protein (negatives, seed 80). Prevalence in the scored set is therefore about 1 in 8, set by sampling and not natural prevalence (about 1 in 150); AUPRC values are relative to that sampled prevalence only.
- Labels include UniProt annotations propagated by similarity or rule (ECO:0000250/0000255) as well as experimental ones; evidence level is not filtered.

## Features (per residue)
- E: ESM-2 35M (esm2_t12_35M_UR50D) last-layer residue embedding (480 dims), full-length forward pass, fp32, CPU (embed_res.py).
- B1: amino-acid identity (20) + relative position. B2 (primary comparator): amino-acid identity in a +/-3 window (140 dims, zero padded). B2 has fewer dims than E (140 vs 480); noted as a limit.

## Method
- Standardize then L2 logistic regression (C = 0.1, max_iter 3000) per feature set. Families sorted, assigned round-robin to 5 folds (no family spans folds); pooled out-of-fold scores.
- Metrics: AUPRC and AUROC on pooled residues. CI: family-cluster bootstrap, 2,000 resamples, seed 12345, percentile 95%.

## Gates
- G1 (control): E trained on within-protein permuted labels has out-of-fold AUPRC < prevalence + 0.03. Failure = INVALID, protocol stop.
- G2: AUPRC(E) - AUPRC(B2) >= +0.10 and CI lower bound > +0.03.
- G3 (reported only, cannot change the label): absolute E AUPRC with CI and E AUROC.
- Label (mechanical): INVALID if G1 fails; EMBEDDING-BEATS-LOCAL-CONTEXT if G2 passes; HONEST NEGATIVE otherwise.

## Limits stated up front
- Active site is a narrow, enzyme-biased target; "motif" here means the annotated catalytic residue, not a sequence motif. Pfam family holdout does not remove clan-level homology between folds.
- Negatives are subsampled, so absolute AUPRC is not a deployment number. Similarity-propagated labels may be recovered through family-level patterns.
- Only 6 proteins per family with at least 6 qualifying; family counts are fixed by the rule and may be fewer than 200; counts were not examined before lock-1 (lesson from DOC-2-079: the realised n will be reported against expectation).
- One model size, last layer only, linear probe. B2 is a local-window baseline, not the strongest sequence baseline. Single seed, CPU only. No simulated data in results; analysis_080.py smoke test was on synthetic noise only, no record saved.
- Single run; crash fixes are dated AMENDMENT-N.md files committed before outcomes exist.
