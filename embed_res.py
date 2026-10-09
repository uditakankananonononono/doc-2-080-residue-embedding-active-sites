"""DOC-2-080: per-residue ESM-2 35M last-layer embeddings for the rows of residues.tsv -> res_emb.npy (row order = residues.tsv)."""
import numpy as np, pandas as pd, torch, time
from transformers import AutoTokenizer, EsmModel
torch.set_num_threads(2); P = pd.read_csv("proteins.tsv", sep="\t"); R = pd.read_csv("residues.tsv", sep="\t")
tok = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D"); m = EsmModel.from_pretrained("facebook/esm2_t12_35M_UR50D").eval()
emb = np.zeros((len(R), m.config.hidden_size), dtype=np.float32); rows = R.groupby("pid").indices; t0 = time.time()
with torch.no_grad():
    for k, pid in enumerate(P.pid):
        b = tok(P.sequence.iloc[pid], return_tensors="pt"); h = m(**b).last_hidden_state[0]
        idx = rows[pid]; emb[idx] = h[R.pos.values[idx] + 1].numpy()
        if k % 100 == 0: print(k, len(P), round(time.time() - t0), flush=True)
np.save("res_emb.npy", emb); print("EMBED_DONE", emb.shape)
