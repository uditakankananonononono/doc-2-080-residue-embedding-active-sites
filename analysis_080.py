"""DOC-2-080 frozen analysis (PROTOCOL.md lock-1). Run once. Needs proteins.tsv, residues.tsv, res_emb.npy."""
import json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import average_precision_score as ap, roc_auc_score as auc
rng = np.random.default_rng(12345)
P = pd.read_csv("proteins.tsv", sep="\t"); R = pd.read_csv("residues.tsv", sep="\t"); E = np.load("res_emb.npy"); y = R.label.values; fam = R.family.values; n = len(R)
seqs = P.set_index("pid").sequence.to_dict(); AA = "ACDEFGHIKLMNPQRSTVWY"; ix = {a: i for i, a in enumerate(AA)}
def window(pid, pos, w):
    s = seqs[pid]; v = np.zeros((2 * w + 1) * 20, dtype=np.float32)
    for o in range(-w, w + 1):
        j = pos + o
        if 0 <= j < len(s): v[(o + w) * 20 + ix[s[j]]] = 1
    return v
B1 = np.array([np.r_[window(p, q, 0), q / len(seqs[p])] for p, q in zip(R.pid, R.pos)]); B2 = np.array([window(p, q, 3) for p, q in zip(R.pid, R.pos)])
uf = sorted(set(fam)); fold_of = {f: i % 5 for i, f in enumerate(uf)}; fid = np.array([fold_of[f] for f in fam])
def oof(X, yy):
    s = np.zeros(n)
    for k in range(5):
        a, b = fid != k, fid == k; m = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=3000)).fit(X[a], yy[a]); s[b] = m.decision_function(X[b])
    return s
S = {"E": oof(E, y), "B1": oof(B1, y), "B2": oof(B2, y)}
yperm = y.copy()
for pid in set(R.pid): m = (R.pid == pid).values; yperm[m] = rng.permutation(y[m])
sperm = oof(E, yperm)
prev = float(y.mean()); res = {"n_residues": n, "n_proteins": len(P), "n_families": len(uf), "prevalence": prev, "AUPRC": {k: float(ap(y, v)) for k, v in S.items()}, "AUROC": {k: float(auc(y, v)) for k, v in S.items()}}
res["G1"] = dict(perm_AUPRC_E=float(ap(yperm, sperm)), prevalence=prev, pass_=bool(ap(yperm, sperm) < prev + 0.03))
members = {u: np.where(fam == u)[0] for u in uf}; bs = {"gap": [], "e": []}
for _ in range(2000):
    idx = np.concatenate([members[u] for u in rng.choice(uf, len(uf))]); yy = y[idx]
    if yy.sum() == 0: continue
    a, b = ap(yy, S["E"][idx]), ap(yy, S["B2"][idx]); bs["gap"].append(a - b); bs["e"].append(a)
ci = lambda k: [float(x) for x in np.percentile(bs[k], [2.5, 97.5])]
gap = res["AUPRC"]["E"] - res["AUPRC"]["B2"]
res["G2"] = dict(delta_AUPRC_E_minus_B2=gap, ci=ci("gap"), pass_=bool(gap >= 0.10 and ci("gap")[0] > 0.03))
res["G3_reported"] = dict(E_AUPRC=res["AUPRC"]["E"], E_AUPRC_ci=ci("e"), E_AUROC=res["AUROC"]["E"])
res["LABEL"] = "INVALID" if not res["G1"]["pass_"] else ("EMBEDDING-BEATS-LOCAL-CONTEXT" if res["G2"]["pass_"] else "HONEST NEGATIVE")
print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))
