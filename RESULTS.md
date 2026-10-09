# DOC-2-080 RESULTS - review pending (independent gate has not cleared; no claim is final)

Label (set mechanically by analysis_080.py): **EMBEDDING-BEATS-LOCAL-CONTEXT**. G1 pass, G2 pass. Prior art: residue-level functional signal in protein language models is known; this is a replication-style measurement with a 35M model, not a novelty claim.

| item | value |
|---|---|
| data | UniProt 2026_03, 65,746 reviewed entries with active-site annotation fetched; 200 families x 6 proteins = 1,200 proteins (meets the protocol's expectation) |
| scored residues | 19,966 (1,966 active-site positives + sampled negatives), sampled prevalence 0.0985 (natural prevalence is about 1 in 150) |
| G1 control (within-protein permuted labels, E OOF AUPRC) | 0.0989 vs prevalence 0.0985; threshold < 0.1285; pass |
| AUPRC (pooled OOF, family-held-out) | E (ESM-2 35M residue embedding) **0.751**, B2 (amino-acid +/-3 window) 0.517, B1 (identity + relative position) 0.419 |
| AUROC | E 0.949, B2 0.908, B1 0.894 |
| G2: AUPRC(E) minus AUPRC(B2) | **+0.234**, 95% family-bootstrap CI [0.192, 0.274]; needed >= 0.10 and CI lb > 0.03: pass |
| G3 (reported only) | E AUPRC 0.751, CI [0.712, 0.787]; E AUROC 0.949 |

## What this shows and does not show
- In held-out Pfam families (5-fold, no family split across folds), a linear probe on 35M-parameter ESM-2 residue embeddings ranks annotated active-site residues far above a strong local-context baseline (AUPRC 0.751 vs 0.517 at sampled prevalence 0.10).
- AUPRC values are relative to the sampled prevalence (about 10%) and are not deployment numbers; at natural prevalence precision would be much lower. AUROC (0.949 vs 0.908) is prevalence-independent in expectation and shows a smaller gap.
- Pfam family holdout does not remove clan-level homology between folds; related families can share catalytic motifs, which would help E and B2 differently. Not measured.
- Labels include UniProt annotations propagated by sequence similarity or rule (ECO:0000250/0000255), not only experimental ones; evidence level was not filtered. A probe may partly learn the propagation pattern rather than catalytic chemistry. Not tested.
- B2 has 140 dims and E has 480; B2 is a local-window baseline, not the strongest sequence baseline (no alignment-based or conservation baseline, no larger-window or learned-sequence baseline). The result shows E beats a local window, not that E beats alignment-based motif tools.
- Narrow target: only annotated active-site residues in enzymes (1 to 5 per protein, 6 proteins per family, 200 families). "Motif" here means the catalytic residue, not a sequence motif. No claim about binding sites or other motifs.
- One model size, last layer only, linear probe, single seed, CPU only.

## Disclosures
- Run once. analysis_080.py md5 0cc61504a26d8757ffdaff768d4ce775 equals the lock-1 file; tag_tree_check.txt records the lock-1 tag tree (5 files). run_log.txt: analysis start UTC 17:18:31, numpy 2.2.6, pandas 2.3.3, sklearn 1.7.2, torch 2.14.1+cpu, transformers 5.19.0, python 3.10.12. run_log.txt lacks the exit and end-time lines: the builder shell hit its time limit after python had printed RESULT_JSON and written results.json (run_note.txt). No rerun.
- input_md5.txt holds md5 of act_raw.tsv, proteins.tsv, residues.tsv, res_emb.npy and the five scripts. embed_res.py ran as locked; no amendments were needed.
- The synthetic smoke test is a builder statement; no record saved.
- Family counts were not examined before lock-1; the realised 200 families and 1,200 proteins matched the expectation (unlike DOC-2-079).
- Large files (act_raw.tsv, proteins.tsv, residues.tsv, res_emb.npy) are in the Drive folder; md5s in input_md5.txt.
