"""DOC-2-080: build proteins.tsv and residues.tsv from act_raw.tsv (frozen sampling rules in PROTOCOL.md)."""
import io, re, numpy as np, pandas as pd
d = pd.read_csv("act_raw.tsv", sep="\t"); d.columns = ["accession", "pfam", "length", "sequence", "act"]
d["pf"] = d.pfam.fillna("").str.strip(";").str.split(";"); d = d[d.pf.str.len() == 1].copy(); d["family"] = d.pf.str[0]
d = d[~d.sequence.str.contains("[XBZUO]")].copy()
def pos(s): return sorted({int(m) - 1 for m in re.findall(r"ACT_SITE (\d+);", str(s))})
d["act_pos"] = d.act.apply(pos); d = d[d.act_pos.apply(len).between(1, 5)].copy()
d = d[d.apply(lambda r: all(0 <= p < len(r.sequence) for p in r.act_pos), axis=1)].sort_values("accession")
cnt = d.family.value_counts(); fams = sorted(cnt[cnt >= 6].index); rng = np.random.default_rng(80)
if len(fams) > 200: fams = sorted(rng.choice(fams, 200, replace=False))
out = []
for f in fams:
    g = d[d.family == f]; out.append(g.iloc[np.sort(rng.choice(len(g), 6, replace=False))])
P = pd.concat(out).reset_index(drop=True); P["pid"] = P.index
R = []
for r in P.itertuples():
    L = len(r.sequence); neg = [i for i in range(L) if i not in set(r.act_pos)]; neg = list(rng.choice(neg, min(15, len(neg)), replace=False))
    for i in sorted(set(r.act_pos) | set(int(x) for x in neg)): R.append((r.pid, r.accession, r.family, i, int(i in r.act_pos)))
R = pd.DataFrame(R, columns=["pid", "accession", "family", "pos", "label"]); P[["pid", "accession", "family", "length", "sequence"]].to_csv("proteins.tsv", sep="\t", index=False); R.to_csv("residues.tsv", sep="\t", index=False)
print("BUILD_DONE families", len(fams), "proteins", len(P), "residues", len(R), "positives", int(R.label.sum()), "prevalence", round(R.label.mean(), 4))
