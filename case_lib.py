"""Utility functions and mock data for generating synthetic genetic cases."""

import random
import uuid
from datetime import datetime

# Mock list of genetic disorders with associated HPO terms and inheritance
MOCK_GENETIC_DISORDERS = [
    {
        "name": "Cystic Fibrosis",
        "hpo": [
            "HP:0001004 (Meconium ileus)",
            "HP:0002094 (Recurrent respiratory infections)",
            "HP:0001999 (Pancreatic insufficiency)",
            "HP:0000855 (Growth retardation)",
        ],
        "inheritance": "Autosomal Recessive",
        "genes": ["CFTR"],
    },
    {
        "name": "Huntington's Disease",
        "hpo": [
            "HP:0002061 (Chorea)",
            "HP:0002511 (Progressive cognitive decline)",
            "HP:0000716 (Depression)",
        ],
        "inheritance": "Autosomal Dominant",
        "genes": ["HTT"],
    },
    {
        "name": "Duchenne Muscular Dystrophy",
        "hpo": [
            "HP:0003735 (Muscle weakness)",
            "HP:0003731 (Gowers' sign)",
            "HP:0003197 (Calf muscle pseudohypertrophy)",
        ],
        "inheritance": "X-linked",
        "genes": ["DMD"],
    },
    {
        "name": "Sickle Cell Anemia",
        "hpo": [
            "HP:0001903 (Anemia)",
            "HP:0001297 (Splenomegaly)",
            "HP:0002130 (Pain crisis)",
        ],
        "inheritance": "Autosomal Recessive",
        "genes": ["HBB"],
    },
    {
        "name": "Marfan Syndrome",
        "hpo": [
            "HP:0002705 (Aortic dissection)",
            "HP:0001083 (Arachnodactyly)",
            "HP:0000501 (Ectopia lentis)",
        ],
        "inheritance": "Autosomal Dominant",
        "genes": ["FBN1"],
    },
    {
        "name": "Phenylketonuria",
        "hpo": [
            "HP:0001994 (Neonatal vomiting)",
            "HP:0001249 (Intellectual disability)",
            "HP:0000739 (Seizures)",
        ],
        "inheritance": "Autosomal Recessive",
        "genes": ["PAH"],
    },
    {
        "name": "Tay-Sachs Disease",
        "hpo": [
            "HP:0001250 (Severe muscular hypotonia)",
            "HP:0001260 (Cherry red spot)",
            "HP:0001276 (Neurodegeneration)",
        ],
        "inheritance": "Autosomal Recessive",
        "genes": ["HEXA"],
    },
    {
        "name": "Neurofibromatosis Type 1",
        "hpo": [
            "HP:0009733 (Cafe-au-lait spots)",
            "HP:0001661 (Neurofibromas)",
            "HP:0000486 (Lisch nodules)",
        ],
        "inheritance": "Autosomal Dominant",
        "genes": ["NF1"],
    },
    {
        "name": "Hemophilia A",
        "hpo": [
            "HP:0001928 (Prolonged bleeding)",
            "HP:0001981 (Hemarthrosis)",
            "HP:0003010 (Factor VIII deficiency)",
        ],
        "inheritance": "X-linked",
        "genes": ["F8"],
    },
    {
        "name": "Alpha-1 Antitrypsin Deficiency",
        "hpo": [
            "HP:0006531 (Emphysema)",
            "HP:0001404 (Hepatic fibrosis)",
            "HP:0009060 (Neonatal cholestasis)",
        ],
        "inheritance": "Autosomal Co-dominant",
        "genes": ["SERPINA1"],
    },
]

# Mock variants with different classifications for a few genes
MOCK_VARIANTS_DB = {
    "CFTR": [
        {"chr": "chr7", "pos": "117559591", "id": "rs113993960", "ref": "G", "alt": "A", "clnsig": "Pathogenic", "freq": "0.0001"},
        {"chr": "chr7", "pos": "117560416", "id": "rs397509172", "ref": "T", "alt": "A", "clnsig": "Pathogenic", "freq": "0.00005"},
        {"chr": "chr7", "pos": "117559871", "id": "rs123456789", "ref": "C", "alt": "G", "clnsig": "Likely_pathogenic", "freq": "0.00002"},
        {"chr": "chr7", "pos": "117561000", "id": "rs987654321", "ref": "A", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.001"},
        {"chr": "chr7", "pos": "117562000", "id": "rs234567890", "ref": "G", "alt": "C", "clnsig": "Benign", "freq": "0.2"},
        {"chr": "chr7", "pos": "117563000", "id": "rs345678901", "ref": "T", "alt": "C", "clnsig": "Likely_benign", "freq": "0.15"},
    ],
    "HTT": [
        {"chr": "chr4", "pos": "3076600", "id": "rsXXXXXXX", "ref": "CAGCAGCAG", "alt": "CAGCAGCAGCAG", "clnsig": "Pathogenic", "freq": "0.00001"},
        {"chr": "chr4", "pos": "3076700", "id": "rsYYYYYYY", "ref": "C", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.0005"},
        {"chr": "chr4", "pos": "3076800", "id": "rsZZZZZZZ", "ref": "A", "alt": "G", "clnsig": "Benign", "freq": "0.3"},
    ],
    "DMD": [
        {"chr": "chrX", "pos": "31234567", "id": "del_exon50", "ref": "N", "alt": "<DEL>", "clnsig": "Pathogenic", "freq": "0.00001"},
        {"chr": "chrX", "pos": "31234700", "id": "rsDMD123", "ref": "C", "alt": "A", "clnsig": "Likely_pathogenic", "freq": "0.000005"},
        {"chr": "chrX", "pos": "31235000", "id": "rsDMD456", "ref": "G", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.0001"},
    ],
    "HBB": [
        {"chr": "chr11", "pos": "5227000", "id": "rs334", "ref": "A", "alt": "T", "clnsig": "Pathogenic", "freq": "0.01"},
        {"chr": "chr11", "pos": "5227500", "id": "rsHBB123", "ref": "G", "alt": "A", "clnsig": "Likely_benign", "freq": "0.1"},
    ],
    "FBN1": [
        {"chr": "chr15", "pos": "48740000", "id": "rsFBN1_path1", "ref": "C", "alt": "T", "clnsig": "Pathogenic", "freq": "0.00001"},
        {"chr": "chr15", "pos": "48740500", "id": "rsFBN1_VUS1", "ref": "A", "alt": "G", "clnsig": "Uncertain_significance", "freq": "0.00005"},
        {"chr": "chr15", "pos": "48741000", "id": "rsFBN1_benign1", "ref": "T", "alt": "C", "clnsig": "Benign", "freq": "0.05"},
    ],
    "PAH": [
        {"chr": "chr12", "pos": "103234567", "id": "rs5030858", "ref": "T", "alt": "C", "clnsig": "Pathogenic", "freq": "0.0001"},
        {"chr": "chr12", "pos": "103234678", "id": "rs62516034", "ref": "C", "alt": "T", "clnsig": "Likely_pathogenic", "freq": "0.00005"},
        {"chr": "chr12", "pos": "103234789", "id": "rs1800000", "ref": "G", "alt": "A", "clnsig": "Uncertain_significance", "freq": "0.001"},
        {"chr": "chr12", "pos": "103234890", "id": "rs199475", "ref": "A", "alt": "G", "clnsig": "Benign", "freq": "0.05"},
    ],
    "HEXA": [
        {"chr": "chr15", "pos": "72346567", "id": "rs11648", "ref": "C", "alt": "T", "clnsig": "Pathogenic", "freq": "0.0001"},
        {"chr": "chr15", "pos": "72346700", "id": "rs28940870", "ref": "G", "alt": "A", "clnsig": "Likely_pathogenic", "freq": "0.00005"},
        {"chr": "chr15", "pos": "72346800", "id": "rs3780103", "ref": "A", "alt": "G", "clnsig": "Uncertain_significance", "freq": "0.001"},
    ],
    "NF1": [
        {"chr": "chr17", "pos": "29553410", "id": "rs79588645", "ref": "C", "alt": "T", "clnsig": "Pathogenic", "freq": "0.0001"},
        {"chr": "chr17", "pos": "29553520", "id": "rs104894193", "ref": "G", "alt": "A", "clnsig": "Likely_pathogenic", "freq": "0.00005"},
        {"chr": "chr17", "pos": "29553630", "id": "rs55954512", "ref": "A", "alt": "G", "clnsig": "Uncertain_significance", "freq": "0.001"},
        {"chr": "chr17", "pos": "29553740", "id": "rs10692721", "ref": "T", "alt": "C", "clnsig": "Benign", "freq": "0.1"},
    ],
    "F8": [
        {"chr": "chrX", "pos": "154835788", "id": "rs61746189", "ref": "A", "alt": "G", "clnsig": "Pathogenic", "freq": "0.0001"},
        {"chr": "chrX", "pos": "154835999", "id": "rs6048", "ref": "G", "alt": "A", "clnsig": "Likely_pathogenic", "freq": "0.00005"},
        {"chr": "chrX", "pos": "154836222", "id": "rs586178", "ref": "T", "alt": "C", "clnsig": "Uncertain_significance", "freq": "0.001"},
        {"chr": "chrX", "pos": "154837333", "id": "rs1800291", "ref": "C", "alt": "T", "clnsig": "Benign", "freq": "0.1"},
    ],
    "SERPINA1": [
        {"chr": "chr14", "pos": "94844347", "id": "rs28929474", "ref": "G", "alt": "A", "clnsig": "Pathogenic", "freq": "0.002"},
        {"chr": "chr14", "pos": "94845432", "id": "rs17580", "ref": "T", "alt": "C", "clnsig": "Likely_pathogenic", "freq": "0.001"},
        {"chr": "chr14", "pos": "94845500", "id": "rs6647", "ref": "C", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.05"},
        {"chr": "chr14", "pos": "94845789", "id": "rs8004738", "ref": "A", "alt": "G", "clnsig": "Benign", "freq": "0.2"},
    ],
}

# Background variants not tied to disease genes
MOCK_BACKGROUND_VARIANTS = (
    [
        {
            "chr": f"chr{random.randint(1,22)}",
            "pos": str(random.randint(100000, 200000000)),
            "id": f"rsBG_{i}",
            "ref": random.choice(["A", "T", "C", "G"]),
            "alt": random.choice(["A", "T", "C", "G"]),
            "clnsig": "Benign",
            "freq": str(random.uniform(0.05, 0.4)),
        }
        for i in range(50)
    ]
    + [
        {
            "chr": f"chr{random.randint(1,22)}",
            "pos": str(random.randint(100000, 200000000)),
            "id": f"rsLBBG_{i}",
            "ref": random.choice(["A", "T", "C", "G"]),
            "alt": random.choice(["A", "T", "C", "G"]),
            "clnsig": "Likely_benign",
            "freq": str(random.uniform(0.01, 0.05)),
        }
        for i in range(20)
    ]
    + [
        {
            "chr": f"chr{random.randint(1,22)}",
            "pos": str(random.randint(100000, 200000000)),
            "id": f"rsVUSBG_{i}",
            "ref": random.choice(["A", "T", "C", "G"]),
            "alt": random.choice(["A", "T", "C", "G"]),
            "clnsig": "Uncertain_significance",
            "freq": str(random.uniform(0.0001, 0.01)),
        }
        for i in range(10)
    ]
)


def select_genetic_disorder_auto():
    """Return a random disorder from the mock list."""
    return random.choice(MOCK_GENETIC_DISORDERS)


def assemble_variants(
    disorder_info,
    inheritance_model,
    num_pathogenic,
    num_vus,
    num_likely_benign,
    num_benign,
    generation_mode,
    proband_sex,
    background_pool=MOCK_BACKGROUND_VARIANTS,
):
    """Assemble variants for a synthetic case."""
    assembled_variants = []
    disease_genes = disorder_info.get("genes", [])
    all_gene_variants = []
    for gene in disease_genes:
        all_gene_variants.extend(MOCK_VARIANTS_DB.get(gene, []))

    pathogenic_candidates = [
        v for v in all_gene_variants if v["clnsig"] in ["Pathogenic", "Likely_pathogenic"]
    ]

    if inheritance_model == "Autosomal Recessive":
        if len(pathogenic_candidates) >= 2:
            chosen = random.sample(pathogenic_candidates, 2)
            assembled_variants.extend(chosen)
            chosen[0]["gt"] = "0/1"
            chosen[1]["gt"] = "0/1"
        elif len(pathogenic_candidates) == 1:
            chosen = random.sample(pathogenic_candidates, 1)
            assembled_variants.extend(chosen)
            assembled_variants.append(pathogenic_candidates[0])
            chosen[0]["gt"] = "0/1"
            assembled_variants[-1]["gt"] = "0/1"
        else:
            return []
    elif inheritance_model == "Autosomal Dominant":
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            chosen["gt"] = "0/1"
            assembled_variants.append(chosen)
        else:
            return []
    elif inheritance_model == "X-linked":
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            if proband_sex == "Male":
                chosen["gt"] = "1/1"
            else:
                chosen["gt"] = "0/1"
            assembled_variants.append(chosen)
        else:
            return []
    else:
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            chosen["gt"] = "0/1"
            assembled_variants.append(chosen)

    random.shuffle(background_pool)
    existing_positions = {(v["chr"], v["pos"]) for v in assembled_variants}
    filtered_background = [
        v for v in background_pool if (v["chr"], v["pos"]) not in existing_positions
    ]

    def add_variants(count, classification):
        added = 0
        for variant in filtered_background:
            if added >= count:
                break
            if variant["clnsig"] == classification and (
                variant["chr"], variant["pos"]
            ) not in existing_positions:
                variant["gt"] = random.choice(["0/1", "1/0", "0/0"])
                assembled_variants.append(variant)
                existing_positions.add((variant["chr"], variant["pos"]))
                added += 1

    add_variants(num_vus, "Uncertain_significance")
    add_variants(num_likely_benign, "Likely_benign")
    add_variants(num_benign, "Benign")

    if generation_mode == "De Novo/Synthetic (randomized)":
        for var in assembled_variants:
            if var["clnsig"] not in ["Pathogenic", "Likely_pathogenic"]:
                try:
                    original_pos = int(var["pos"])
                    var["pos"] = str(max(1, original_pos + random.randint(-500, 500)))
                except ValueError:
                    pass

    return assembled_variants


def generate_patient_metadata(sex, age, proband_id_base, ethnicity, generate_ped):
    proband_id = proband_id_base if proband_id_base else f"SYNTH_PATIENT_{str(uuid.uuid4())[:8].upper()}"
    return {
        "Sex": sex,
        "Age": age,
        "Proband ID": proband_id,
        "Ethnicity": ethnicity,
        "Family Structure Generated": generate_ped,
    }


def generate_ped_file_content(proband_id, sex, disease_name):
    family_id = f"FAM_{proband_id}"
    ped_content = [
        f"{family_id}\t{proband_id}\t{proband_id}_father\t{proband_id}_mother\t{1 if sex == 'Male' else 2}\t2 # Affected by {disease_name}",
        f"{family_id}\t{proband_id}_father\t0\t0\t1\t1 # Unaffected",
        f"{family_id}\t{proband_id}_mother\t0\t0\t2\t1 # Unaffected",
    ]
    return "\n".join(ped_content)


def generate_vcf_file_content(variants, patient_metadata):
    proband_id = patient_metadata["Proband ID"]
    file_date = datetime.now().strftime("%Y%m%d")
    header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={file_date}",
        "##source=SyntheticPatientCaseGenerator",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth (mock value)\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency (mock value or actual frequency)\">",
        "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{proband_id}",
    ]

    records = []
    for var in variants:
        info = (
            f"DP=100;AF={var.get('freq', '0')};CLNSIG={var.get('clnsig', 'Uncertain_significance').replace(' ', '_')}"
        )
        records.append(
            f"{var.get('chr','.')}\t{var.get('pos','.')}\t{var.get('id','.')}\t{var.get('ref','.')}\t{var.get('alt','.')}\t.\tPASS\t{info}\tGT\t{var.get('gt','0/1')}"
        )

    return "\n".join(header + records)
