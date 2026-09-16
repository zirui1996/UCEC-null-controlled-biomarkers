"""§23 HPA IHC for core-4 panel proteins (PEG3, RUNX1T1).
HAND2-AS1 is lncRNA (no protein); OSR2 IHC reliability = Uncertain -> not used.
Parses HPA XML, reports normal endometrium + endometrial cancer staining for
every antibody (all reported, not just favourable), downloads one image per
(antibody, group, staining level)."""
import json, os, re, subprocess, sys, urllib.request, xml.etree.ElementTree as ET
from collections import Counter

PROJ = os.environ.get("UCEC_PROJ", os.getcwd())
GENES = {"PEG3": "ENSG00000198300", "RUNX1T1": "ENSG00000079102"}
OUT_IMG = f"{PROJ}/results/protein/images_panel"
os.makedirs(OUT_IMG, exist_ok=True)

def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return
    subprocess.run(["curl", "-s", "-m", "120", "-o", dest, url], check=True)

summary, manifest = [], []
for gene, ens in GENES.items():
    xml = f"{PROJ}/data_raw/HPA/{ens}.xml"
    os.makedirs(os.path.dirname(xml), exist_ok=True)
    fetch(f"https://www.proteinatlas.org/{ens}.xml", xml)
    root = ET.parse(xml).getroot()
    # IHC reliability (tissue section)
    rel = [(v.get("type"), v.text) for v in root.iter("verification")]
    print(f"\n===== {gene} ({ens}) =====")
    print("  reliability/verification:", rel if rel else "(none in XML)")
    for ab in root.iter("antibody"):
        abid = ab.get("id")
        for te in ab.findall("tissueExpression"):
            at = te.get("assayType")
            if at not in ("tissue", "cancer"):
                continue
            for d in te.findall("data"):
                t = d.find("tissue")
                if t is None or "ndometri" not in (t.text or ""):
                    continue
                tissue = t.text
                group = "normal" if at == "tissue" else "cancer"
                cells = {c.findtext("cellType"): [(l.get("type"), l.get("count"), l.text)
                         for l in c.findall("level")] for c in d.findall("tissueCell")}
                pats = d.findall("patient")
                lv = Counter(p.findtext("level[@type='staining']") for p in pats)
                print(f"  {abid} | {group:6s} | {tissue:20s} | n_patients={len(pats)} | per-patient staining={dict(lv)}")
                for ct, v in cells.items():
                    print(f"      {ct}: {v}")
                summary.append(dict(gene=gene, antibody=abid, group=group, tissue=tissue,
                                    n_patients=len(pats), patient_staining=dict(lv), cells=cells))
                # one image per staining level
                seen = set()
                for p in pats:
                    s = p.findtext("level[@type='staining']") or "NA"
                    if s in seen:
                        continue
                    url = p.findtext(".//imageUrl")
                    if not url:
                        continue
                    seen.add(s)
                    tag = re.sub(r"\W", "", tissue)
                    fn = f"{gene}_{abid}_{group}_{tag}_{s.replace(' ', '')}.jpg"
                    fetch(url, f"{OUT_IMG}/{fn}")
                    manifest.append(dict(gene=gene, antibody=abid, group=group, tissue=tissue,
                                         staining=s, location=p.findtext("location"),
                                         url=url, file=f"results/protein/images_panel/{fn}"))

json.dump(summary, open(f"{PROJ}/results/protein/hpa_panel_summary.json", "w"), indent=1)
json.dump(manifest, open(f"{PROJ}/results/protein/HPA_panel_image_manifest.json", "w"), indent=1)
print(f"\n{len(manifest)} images -> {OUT_IMG}")
print("===== 35_hpa_panel DONE =====")
