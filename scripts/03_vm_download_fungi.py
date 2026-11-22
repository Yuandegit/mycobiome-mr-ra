#!/usr/bin/env python3
"""Download all 115 genus-level fungi manhattan JSONs (Westlake pheweb).
Sets a browser User-Agent (default urllib UA is blocked with 403)."""
import json, os, time, urllib.parse, urllib.request

ROOT = "/home/user/mycobiome_joint_mr"
OUT = os.path.join(ROOT, "data/gwas_microbe/fungi_genus")
os.makedirs(OUT, exist_ok=True)
BASE = "https://omics.lab.westlake.edu.cn/data/fungi/api/manhattan/pheno/"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")

with open(os.path.join(ROOT, "data/genus_phenocodes.txt")) as f:
    codes = [l.strip() for l in f if l.strip()]

mapping = {}
for i, code in enumerate(codes):
    idx = f"{i:03d}"
    path = os.path.join(OUT, idx + ".json")
    mapping[idx] = code
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        continue
    url = BASE + urllib.parse.quote(code, safe="") + ".json"
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                               "Referer": "https://omics.lab.westlake.edu.cn/"})
    ok = False
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            with open(path, "wb") as fh:
                fh.write(data)
            ok = True
            break
        except Exception as e:
            time.sleep(2 + attempt * 2)
    print(f"{idx}  {'OK' if ok else 'FAIL'}  {len(data) if ok else 0} bytes  {code[:60]}")

with open(os.path.join(OUT, "mapping.json"), "w") as f:
    json.dump(mapping, f, indent=1)
print("saved mapping.json; total", len(mapping))
