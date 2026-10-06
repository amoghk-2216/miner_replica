'''
Settings shared by every notebook: which cancer to run, where the data lives and the per-cancer choices.
Edit this file (or set the environment variables named below); the notebooks themselves need no editing.
'''

import os

REPO_DIR = os.path.dirname(os.path.abspath(__file__))

# TCGA project code of the cancer to run (BRCA, KIRC, LUAD, UCEC, THCA, PRAD, LUSC, HNSC, COAD, LGG, GBM, ...)
CANCER = os.environ.get("MINER_CANCER", "GBM")

# Data root. Each cancer gets <BASE_DIR>/<cancer>/<cancer>_data and <cancer>_output; shared files sit directly in BASE_DIR:
#   mc3.v0.2.8.PUBLIC.maf   TCGA MC3 somatic mutations (all cancers)
#   SI_geneset.gmt          ESTIMATE stromal / immune signatures (optional, used by 07)
BASE_DIR = os.environ.get("MINER_BASE_DIR", os.path.join(REPO_DIR, "data"))

# Folder holding the Gene Ontology files go-basic.obo and goa_human.gaf (shared across cancers, used by 04)
GO_DIR = os.environ.get("MINER_GO_DIR", REPO_DIR)

# TCGA MC3 somatic mutation file (all cancers). By default the copy in the pan-cancer project folder next to this one
MAF_PATH = os.environ.get("MINER_MAF_PATH", os.path.join(REPO_DIR, "..", "causal_inf_new", "mc3.v0.2.8.PUBLIC.maf"))

# Notebook 01 (GBM). Genes whose microarray and RNA-seq measurements correlate less than this, across the patients measured
# on both platforms, are dropped before the platforms are combined. 0 drops only the genes whose two
# measurements correlate negatively (589 of 10,379 in GBM); -1 keeps every gene. The value is a judgement
MIN_CROSS_PLATFORM_R = 0

# Notebook 02. Pathway mutations (as in the paper): a pathway is mutated in a patient if any member gene is
# PATHWAY_PREFIX selects the pathway collection in miner's database ("PID_" = NCI-Nature, the one the paper uses)
# Pathways mutated in more than MAX_PATHWAY_FREQUENCY of patients are dropped (too few wild-type patients to compare with)
INCLUDE_PATHWAY_MUTATIONS = True
PATHWAY_PREFIX = "PID_"
MAX_PATHWAY_FREQUENCY = 0.95

# Notebook 02. False = the published gbmMINER rule: every gene mutated in >= 5% of patients is tested for causal flows
# True = only published driver genes of the cancer are tested (see driver_tiers.py); much faster in heavily mutated cancers
RESTRICT_TO_DRIVERS = False

# Notebook 03. cBioPortal study whose data_clinical_patient.txt is used, placed in <cancer>_data/<study>/
# (colon cancer is combined with rectal in cBioPortal)
CBIOPORTAL_STUDIES = {
    "BRCA": "brca_tcga", "KIRC": "kirc_tcga", "LUAD": "luad_tcga", "UCEC": "ucec_tcga", "LUSC": "lusc_tcga", "GBM": "gbm_tcga",
    "THCA": "thca_tcga_pan_can_atlas_2018", "PRAD": "prad_tcga_pan_can_atlas_2018", "HNSC": "hnsc_tcga_pan_can_atlas_2018",
    "COAD": "coadread_tcga_pan_can_atlas_2018", "LGG": "lgg_tcga_pan_can_atlas_2018",
}

# Notebook 03. Survival endpoint per cancer; cancers not listed use overall survival ("OS")
# Progression-free survival ("PFS") for cancers with too few deaths; it is only in the PanCancer Atlas clinical files
ENDPOINTS = {"THCA": "PFS", "PRAD": "PFS"}

# Notebooks 04, 05, 07. Independent cohort per cancer: cancer -> (name, folder inside <cancer>_data)
# A cancer that is not listed is run without external validation and its results are marked as such
VALIDATION_COHORTS = {
    "GBM":  ("GRAVENDEEL", "gbm_gravendeel"),
    "BRCA": ("METABRIC", "brca_metabric"),
    "KIRC": ("EMTAB1980", "kirc_emtab1980"),
    "LUAD": ("GSE72094", "luad_gse72094"),
    "UCEC": ("CPTAC", "ucec_cptac_2020"),
}
