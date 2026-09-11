#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import pandas as pd
import requests


URL = "https://api.platform.opentargets.org/api/v4/graphql"
TARGETS = {
    "GPR25": "ENSG00000170128",
    "MST1": "ENSG00000173531",
    "FADS1": "ENSG00000149485",
    "TMEM258": "ENSG00000134825",
}
REGIONS = {
    "block154": ("1", 200_134_006, 201_067_952),
    "block464": ("3", 47_588_462, 50_387_742),
    "block1671": ("11", 60_515_106, 61_717_117),
}

Q_TARGET = """
query($ensemblId: String!, $pageIdx: Int!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    credibleSets(page: {index: $pageIdx, size: 500}) {
      count
      rows {
        studyLocusId studyId studyType chromosome position qtlGeneId isTransQtl
        variant { id rsIds }
      }
    }
  }
}
"""

Q_COLOC = """
query($studyLocusId: String!, $pageIdx: Int!) {
  credibleSet(studyLocusId: $studyLocusId) {
    studyLocusId
    colocalisation(page: {index: $pageIdx, size: 500}) {
      count
      rows {
        h3 h4 clpp colocalisationMethod rightStudyType numberColocalisingVariants
        otherStudyLocus {
          studyId studyType studyLocusId chromosome position qtlGeneId
          variant { id rsIds }
        }
      }
    }
  }
}
"""

Q_STUDY = """
query($studyId: String!) {
  study(studyId: $studyId) {
    id traitFromSource nSamples nCases nControls pubmedId publicationFirstAuthor publicationDate
  }
}
"""


def post(query: str, variables: dict, retries: int = 3) -> dict:
    last = None
    for attempt in range(retries):
        try:
            response = requests.post(URL, json={"query": query, "variables": variables}, timeout=90)
            response.raise_for_status()
            payload = response.json()
            if payload.get("errors"):
                raise RuntimeError(payload["errors"][0]["message"])
            return payload["data"]
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(str(last))


def region_for(chromosome: str, position: int) -> str | None:
    for locus, (chrom, start, stop) in REGIONS.items():
        if str(chromosome) == chrom and start - 1_000_000 <= int(position) <= stop + 1_000_000:
            return locus
    return None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    out = args.root / "results" / "opentargets"
    raw = args.root / "data" / "opentargets"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)

    query_log = []
    rows = []
    study_ids = set()
    for symbol, ensembl_id in TARGETS.items():
        all_sets = []
        page = 0
        total = None
        while True:
            data = post(Q_TARGET, {"ensemblId": ensembl_id, "pageIdx": page})
            target = data.get("target")
            if not target:
                break
            credible = target["credibleSets"]
            total = credible["count"]
            page_rows = credible["rows"]
            all_sets.extend(page_rows)
            if not page_rows or len(all_sets) >= total:
                break
            page += 1
            time.sleep(0.4)
        pqtl_sets = [item for item in all_sets if item.get("studyType") == "pqtl"]
        query_log.append({
            "gene_symbol": symbol, "ensembl_id": ensembl_id,
            "all_target_credible_sets": len(all_sets), "pqtl_credible_sets": len(pqtl_sets),
            "target_pages": page + 1,
        })
        with open(raw / f"{symbol}_pqtl_credible_sets.json", "w") as handle:
            json.dump(pqtl_sets, handle, indent=2)

        for index, qtl_set in enumerate(pqtl_sets):
            coloc_rows = []
            coloc_page = 0
            coloc_total = None
            while True:
                data = post(Q_COLOC, {
                    "studyLocusId": qtl_set["studyLocusId"], "pageIdx": coloc_page,
                })
                credible = data.get("credibleSet")
                if not credible:
                    break
                coloc = credible["colocalisation"]
                coloc_total = coloc["count"]
                page_rows = coloc["rows"]
                coloc_rows.extend(page_rows)
                if not page_rows or len(coloc_rows) >= coloc_total:
                    break
                coloc_page += 1
                time.sleep(0.25)
            for coloc in coloc_rows:
                if coloc.get("rightStudyType") != "gwas":
                    continue
                other = coloc.get("otherStudyLocus") or {}
                locus = region_for(other.get("chromosome", ""), other.get("position", -1))
                if locus is None:
                    continue
                study_ids.add(other.get("studyId", ""))
                qvar = qtl_set.get("variant") or {}
                gvar = other.get("variant") or {}
                rows.append({
                    "locus": locus, "target_gene": symbol, "target_ensembl_id": ensembl_id,
                    "pqtl_study_locus_id": qtl_set.get("studyLocusId", ""),
                    "pqtl_study_id": qtl_set.get("studyId", ""),
                    "pqtl_variant_rsids": ",".join(qvar.get("rsIds") or []),
                    "pqtl_chr": qtl_set.get("chromosome", ""), "pqtl_pos": qtl_set.get("position", ""),
                    "pqtl_is_trans": qtl_set.get("isTransQtl", ""),
                    "gwas_study_locus_id": other.get("studyLocusId", ""),
                    "gwas_study_id": other.get("studyId", ""),
                    "gwas_variant_rsids": ",".join(gvar.get("rsIds") or []),
                    "gwas_chr": other.get("chromosome", ""), "gwas_pos": other.get("position", ""),
                    "h3": coloc.get("h3", ""), "h4": coloc.get("h4", ""),
                    "clpp": coloc.get("clpp", ""), "method": coloc.get("colocalisationMethod", ""),
                    "n_colocalising_variants": coloc.get("numberColocalisingVariants", ""),
                })
            if index % 50 == 0:
                time.sleep(0.3)

    metadata = {}
    for study_id in sorted(value for value in study_ids if value):
        try:
            metadata[study_id] = post(Q_STUDY, {"studyId": study_id}).get("study") or {}
        except Exception as exc:
            metadata[study_id] = {"metadata_error": str(exc)}
    with open(raw / "gwas_study_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)

    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["gwas_trait"] = frame.gwas_study_id.map(
            lambda value: metadata.get(value, {}).get("traitFromSource", "")
        )
        trait = frame.gwas_trait.fillna("").str.lower()
        frame["gwas_family"] = "other"
        frame.loc[trait.str.contains(r"inflammatory bowel|crohn|ulcerative colitis|\bibd\b", regex=True), "gwas_family"] = "IBD"
        frame.loc[trait.str.contains("depress", regex=True), "gwas_family"] = "depression"
        frame["h4"] = pd.to_numeric(frame.h4, errors="coerce")
        frame["h3"] = pd.to_numeric(frame.h3, errors="coerce")
        frame = frame.sort_values(["locus", "target_gene", "gwas_family", "h4"], ascending=[True, True, True, False])
    frame.to_csv(out / "opentargets_pqtl_coloc.tsv", sep="\t", index=False)
    pd.DataFrame(query_log).to_csv(out / "opentargets_query_audit.tsv", sep="\t", index=False)
    manifest = []
    for path in sorted(list(out.glob("*")) + list(raw.glob("*"))):
        if path.is_file():
            manifest.append({"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(out / "opentargets_file_manifest.tsv", sep="\t", index=False)
    print(pd.DataFrame(query_log).to_string(index=False))
    print("region-matched pQTL-GWAS rows:", len(frame))
    if not frame.empty:
        print(frame.loc[frame.gwas_family.ne("other")].to_string(index=False))


if __name__ == "__main__":
    main()
