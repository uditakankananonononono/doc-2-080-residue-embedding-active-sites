"""DOC-2-080: fetch reviewed UniProt entries with annotated active sites, Pfam, length 100-400. Writes act_raw.tsv, DATA_HASHES.tsv."""
import time, hashlib, requests
Q = "reviewed:true AND ft_act_site:* AND length:[100 TO 400] AND xref:pfam-*"
url = "https://rest.uniprot.org/uniprotkb/search"; params = dict(query=Q, fields="accession,xref_pfam,length,sequence,ft_act_site", format="tsv", size=500)
rows, rel, first = [], None, True
while url:
    for attempt in range(6):
        try:
            r = requests.get(url, params=params if first else None, timeout=120); r.raise_for_status(); break
        except requests.exceptions.RequestException as ex:
            print("retry", attempt, type(ex).__name__, flush=True); time.sleep(5 * (attempt + 1))
    else: raise SystemExit("acquisition failed")
    first = False; rel = r.headers.get("x-uniprot-release", rel); lines = [x for x in r.text.split("\n") if x]; rows += lines[1:] if rows else lines
    url = r.links.get("next", {}).get("url"); print("page", len(rows), flush=True)
raw = "\n".join(rows) + "\n"; open("act_raw.tsv", "w").write(raw)
open("DATA_HASHES.tsv", "w").write(f"file\tmd5\tuniprot_release\tn_rows\nact_raw.tsv\t{hashlib.md5(raw.encode()).hexdigest()}\t{rel}\t{len(rows)-1}\n")
print("ACT_DONE release", rel, "rows", len(rows) - 1)
