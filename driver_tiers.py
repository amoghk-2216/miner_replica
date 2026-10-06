'''
Driver gene tiers for the MINER pan-cancer notebooks (used by 02, 06 and 07 of every cancer).

Why: MINER's causal step tests every gene mutated in >= 5% of patients. In heavily mutated cancers (LUAD: 785 genes,
UCEC: 5,027) almost all of those are passengers that track total mutation burden, so their "causal" flows are confounded.
Frequency alone cannot separate drivers from passengers (Lawrence et al. 2013, Nature), and a naive
nonsynonymous / synonymous ratio is biased (Martincorena et al. 2017, Cell). So the gene list is taken from published
driver calls, fetched per cancer from the Open Targets Platform, which carries two citable sources as separate channels:

A gene counts as a driver of a cancer only if IntOGen calls it one in that cancer (statistical evidence of selection in
tumour cohorts of that type, incl. dNdScv; Martinez-Jimenez et al. 2020). The Cancer Gene Census channel (Sondka et al. 2018)
then splits those genes into two tiers:

    known      IntOGen driver call for this cancer AND a catalogued cancer gene (Cancer Gene Census evidence)
    candidate  IntOGen driver call for this cancer, NOT in the Cancer Gene Census (the novel side)

The Cancer Gene Census channel is not used on its own: Open Targets gives it a score for every catalogued cancer gene with
somatic mutations observed in the tumour type (745 genes for breast carcinoma, MUC16 among them), so it is not specific.

Everything else is excluded. The result is cached as driver_genes.csv in the cancer's data folder, so the three
notebooks use exactly the same list; delete that file (or pass refresh=True) to fetch again.
'''

import os
import time
import pandas as pd
import requests

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Open Targets disease search term per cancer; OT_DISEASE_IDS forces a specific ID where the search is ambiguous
OT_QUERIES = {
    "BRCA": "breast carcinoma", "KIRC": "clear cell renal carcinoma", "LUAD": "lung adenocarcinoma",
    "UCEC": "endometrial carcinoma", "THCA": "thyroid carcinoma", "PRAD": "prostate adenocarcinoma",
    "LUSC": "lung squamous cell carcinoma", "HNSC": "head and neck squamous cell carcinoma",
    "COAD": "colon adenocarcinoma", "LGG": "low grade glioma",
    "GBM": "glioblastoma",
}
OT_DISEASE_IDS = {"BRCA": "MONDO_0004989"}

# Minimum IntOGen score (Open Targets scale, 0-1) for a gene to count as a driver. IntOGen calls with a low score include
# genes widely regarded as passengers (in lung adenocarcinoma: CSMD3 0.28, LRP1B 0.57, FAT4 0.17), while the established
# drivers score 0.68 or more in the three cancers inspected (BRCA, KIRC, LUAD). 0.6 sits in that gap; it is a judgement,
# not a published threshold. Set to 0 to accept every IntOGen call. The cache keeps every call, so changing this needs no new download
MIN_INTOGEN_SCORE = 0.6

# Fallback only (general knowledge, NOT a citable source): used when Open Targets cannot be reached and no cache exists
PROVISIONAL_DRIVERS = {
    "LUAD": ["KRAS", "EGFR", "TP53", "STK11", "KEAP1", "NF1", "BRAF", "RBM10", "SMARCA4", "ATM", "RB1", "CDKN2A",
             "ARID1A", "SETD2", "MET", "ERBB2", "U2AF1", "RIT1", "MGA", "PIK3CA"],
    "KIRC": ["VHL", "PBRM1", "SETD2", "BAP1", "KDM5C", "MTOR"],
    "BRCA": ["TP53", "PIK3CA", "GATA3", "CDH1", "MAP3K1", "PTEN", "KMT2C", "NCOR1", "MAP2K4", "RUNX1", "CBFB", "AKT1",
             "NF1", "RB1", "ARID1A", "FOXA1", "TBX3", "ERBB2"],
    "UCEC": ["PTEN", "PIK3CA", "PIK3R1", "ARID1A", "CTNNB1", "TP53", "KRAS", "CTCF", "FBXW7", "PPP2R1A", "POLE", "FGFR2",
             "ARID5B", "RPL22", "KMT2D", "ZFHX3", "SPOP", "NFE2L2", "SOX17", "MSH6", "CCND1", "KMT2B"],
    "THCA": ["BRAF", "NRAS", "HRAS", "KRAS", "EIF1AX", "TG", "PPM1D", "CHEK2"],
}


def _post(query, variables, retries=4):

    '''

    Sends one GraphQL query to Open Targets, retrying on temporary server errors.

    Args: query (str). GraphQL query text
          variables (dict). Query variables
          retries (int). Number of attempts (default 4)

    Returns: dict. The "data" part of the response

    '''

    for attempt in range(retries):
        resp = requests.post(OT_URL, json={"query": query, "variables": variables}, timeout=60)
        if resp.status_code in (429, 502, 503, 504):
            time.sleep(5 * (attempt + 1))
            continue
        resp.raise_for_status()
        body = resp.json()
        if body.get("errors"):
            raise RuntimeError(f"Open Targets query error: {body['errors'][0].get('message')}")
        return body["data"]
    raise RuntimeError("Open Targets did not respond")


def resolve_disease(cancer):

    '''

    Finds the Open Targets disease ID for a cancer (fixed ID if listed in OT_DISEASE_IDS, otherwise the best search hit).

    Args: cancer (str). TCGA project code, e.g. "LUAD"

    Returns: disease_id (str), disease_name (str)

    '''

    q = """query s($q: String!) { search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 10}) { hits { id name } } }"""
    term = OT_QUERIES[cancer]
    hits = _post(q, {"q": term})["search"]["hits"]
    exact = [h for h in hits if h["name"].lower() == term.lower()]
    print(f"Open Targets search hits for '{term}': {[(h['id'], h['name']) for h in hits[:5]]}")

    # Candidate IDs in order of preference: the fixed ID, an exact name match, then the remaining hits
    # A fixed ID can stop existing when Open Targets renames diseases, so each candidate is checked
    candidates = [OT_DISEASE_IDS[cancer]] if cancer in OT_DISEASE_IDS else []
    candidates += [h["id"] for h in exact] + [h["id"] for h in hits]
    check = """query d($efoId: String!) { disease(efoId: $efoId) { id name } }"""
    for disease_id in dict.fromkeys(candidates):
        found = _post(check, {"efoId": disease_id})["disease"]
        if found is not None:
            return found["id"], found["name"]
        print(f"Open Targets has no disease with ID {disease_id}; trying the next candidate")
    raise LookupError(f"no Open Targets disease found for '{term}'")


def fetch_associations(disease_id):

    '''

    Downloads every target associated with a disease in Open Targets, with the overall score and the scores of the
    two driver-related data sources.

    Args: disease_id (str). Open Targets disease ID

    Returns: df, one row per gene. Columns Gene, ot_score (overall, 0-1), cgc_score (Cancer Gene Census), intogen_score

    '''

    q = """query t($efoId: String!, $page: Int!) { disease(efoId: $efoId) { associatedTargets(page: {index: $page, size: 500}) {
             count rows { target { approvedSymbol } score datasourceScores { id score } } } } }"""
    rows, page = [], 0
    while True:
        disease = _post(q, {"efoId": disease_id, "page": page})["disease"]
        if disease is None:
            raise LookupError(f"Open Targets has no disease with ID {disease_id}")
        data = disease["associatedTargets"]
        if not data["rows"]:
            break
        for r in data["rows"]:
            sources = {s["id"]: s["score"] for s in (r.get("datasourceScores") or [])}
            rows.append({"Gene": r["target"]["approvedSymbol"], "ot_score": r["score"],
                         "cgc_score": sources.get("cancer_gene_census", 0.0), "intogen_score": sources.get("intogen", 0.0)})
        page += 1
        time.sleep(0.3)
    return pd.DataFrame(rows).drop_duplicates("Gene")


def _apply_min_score(table):

    '''

    Keeps the driver calls whose IntOGen score reaches MIN_INTOGEN_SCORE.

    Args: table (df). Driver table with columns Gene, Tier, intogen_score, ...

    Returns: df, same columns. Rows with intogen_score >= MIN_INTOGEN_SCORE

    '''

    kept = table[table["intogen_score"] >= MIN_INTOGEN_SCORE].reset_index(drop=True)
    print(f"IntOGen score >= {MIN_INTOGEN_SCORE}: {len(kept)} of {len(table)} driver calls kept: {kept['Tier'].value_counts().to_dict()}")
    return kept


def get_driver_tiers(cancer, data_dir, refresh=False):

    '''

    Returns the driver gene table for a cancer, from the cache if present, otherwise from Open Targets.

    Args: cancer (str). TCGA project code, e.g. "LUAD"
          data_dir (str). The cancer's data folder; driver_genes.csv and open_targets_associations.csv are cached here
          refresh (bool). True ignores the cache and fetches again (default False)

    Returns: df, one row per driver gene. Columns Gene, Tier ("known" / "candidate" / "provisional"), cgc_score,
             intogen_score, ot_score, Source, Disease_ID. Tier "provisional" means Open Targets could not be reached and
             the built-in fallback list was used (not cached, so the next run tries again)

    '''

    cache = os.path.join(data_dir, "driver_genes.csv")
    if os.path.exists(cache) and not refresh:
        table = pd.read_csv(cache)
        print(f"Driver genes for {cancer} loaded from cache: {table['Tier'].value_counts().to_dict()}")
        return _apply_min_score(table)

    try:
        disease_id, disease_name = resolve_disease(cancer)
        assoc = fetch_associations(disease_id)
        assoc.to_csv(os.path.join(data_dir, "open_targets_associations.csv"), index=False)

        # A driver must be called by IntOGen for this cancer; the Cancer Gene Census channel only splits the tiers
        called    = assoc[assoc["intogen_score"] > 0]
        known     = called[called["cgc_score"] > 0].assign(Tier="known", Source="IntOGen + Cancer Gene Census (via Open Targets)")
        candidate = called[called["cgc_score"] == 0].assign(Tier="candidate", Source="IntOGen only (via Open Targets)")
        table = pd.concat([known, candidate]).assign(Disease_ID=disease_id)
        table = table[["Gene", "Tier", "cgc_score", "intogen_score", "ot_score", "Source", "Disease_ID"]]
        table = table.sort_values("intogen_score", ascending=False).reset_index(drop=True)
        if table.empty:
            raise RuntimeError(f"no IntOGen driver calls returned for {disease_id}")
        table.to_csv(cache, index=False)
        print(f"Driver genes for {cancer} from Open Targets disease {disease_id} ({disease_name}): {table['Tier'].value_counts().to_dict()}")
        print("CHECK that this is the right disease; if not, set its ID in OT_DISEASE_IDS in driver_tiers.py and delete driver_genes.csv")
        return _apply_min_score(table)

    except Exception as e:
        genes = PROVISIONAL_DRIVERS.get(cancer)
        print(f"WARNING: Open Targets driver calls unavailable ({type(e).__name__}: {e})")
        if genes is None:
            print(f"No fallback list for {cancer}: no driver filter can be applied")
            return pd.DataFrame(columns=["Gene", "Tier", "cgc_score", "intogen_score", "ot_score", "Source", "Disease_ID"])
        print(f"Using the provisional built-in list for {cancer} ({len(genes)} genes). This is NOT a citable source; rerun with internet access")
        return pd.DataFrame({"Gene": genes, "Tier": "provisional", "cgc_score": float("nan"), "intogen_score": float("nan"),
                             "ot_score": float("nan"), "Source": "provisional built-in list", "Disease_ID": ""})
