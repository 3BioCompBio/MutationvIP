import polars as pl


CYTOBAND_COLOR = {
    "gneg": (1.0, 1.0, 1.0),
    "gpos25": (0.6, 0.6, 0.6),
    "gpos50": (0.4, 0.4, 0.4),
    "gpos75": (0.2, 0.2, 0.2),
    "gpos100": (0.0, 0.0, 0.0),
    "acen": (0.8, 0.4, 0.4),
    "gvar": (0.8, 0.8, 0.8),
    "stalk": (0.9, 0.9, 0.9),
}

CHROMOSOMES = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
]

DBS = {"Ssc": "cosmic", "Ssn": "somamutdb", "Sg": "clinvar"}

DISEASE_MAP = {
    "thymus": "UNMAPPED",
    "oesophagus": "Esophagus",
    "skin": "Skin",
    "gastrointestinal_tract_(site_indeterminate)": "UNMAPPED",
    "genital_tract": "UNMAPPED",
    "fallopian_tube": "UNMAPPED",
    "upper_aerodigestive_tract": "UNMAPPED",
    "kidney": "Kidney",
    "soft_tissue": "Skeletal muscle",
    "adrenal_gland": "Adrenal gland",
    "central_nervous_system": "Brain",
    "vulva": "UNMAPPED",
    "placenta": "Placenta",
    "large_intestine": "Colon",
    "biliary_tract": "UNMAPPED",
    "pancreas": "Pancreas",
    "salivary_gland": "UNMAPPED",
    "eye": "UNMAPPED",
    "breast": "Breast",
    "pituitary": "UNMAPPED",
    "ovary": "UNMAPPED",
    "endometrium": "Endometrium",
    "haematopoietic_and_lymphoid_tissue": "Blood",
    "thyroid": "Thyroid",
    "autonomic_ganglia": "UNMAPPED",
    "lung": "Lung",
    "pleura": "UNMAPPED",
    "parathyroid": "UNMAPPED",
    "peritoneum": "UNMAPPED",
    "liver": "Liver",
    "cervix": "UNMAPPED",
    "penis": "UNMAPPED",
    "urinary_tract": "Bladder",
    "bone": "Bone marrow",
    "NS": "UNMAPPED",
    "small_intestine": "Small intestine",
    "prostate": "Prostate",
    "meninges": "UNMAPPED",
    "vagina": "UNMAPPED",
    "stomach": "Stomach",
    "testis": "Testicle",
}


TYPE_ORDER = [
    "nonsence",
    "missence",
    "silent",
    "splice_region",
    "intronic",
    "UTR",
    "intergenic",
]
TYPE_ORDER = [
    "silent",
    "missence",
    "nonsence",
    "intronic",
    "splice_region",
    "UTR",
    "intergenic",
]
DF_ORDER = ["Ssc", "Ssn", "Sg"]
IMPACT_ORDER = ["Unknown", "Benign", "Pathogenic"]

GENOMIC_CONTEXT_PROBA = pl.DataFrame(
    {
        "type": [
            "intergenic",
            "transcript",
            "intronic",
            "UTR",
            "missence",
            "nonsence",
            "silent",
        ],
        "proba_type": [0.65, 0.35, 0.315, 0.01, 0.0183, 0.0018, 0.0055],
    }
)

CODING_PROBA = pl.DataFrame(
    {"type": ["missence", "nonsence", "silent"], "proba_type": [0.73, 0.05, 0.22]}
)
