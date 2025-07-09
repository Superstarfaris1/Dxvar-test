import streamlit as st
import pandas as pd
import numpy as np
import random
from datetime import datetime
import uuid # For generating unique IDs

# --- Constants and Mock Data (REPLACE WITH YOUR REAL DATA SOURCES) ---

# Mock list of genetic disorders with associated HPO terms and inheritance
MOCK_GENETIC_DISORDERS = [
    {
        "name": "Cystic Fibrosis",
        "hpo": ["HP:0001004 (Meconium ileus)", "HP:0002094 (Recurrent respiratory infections)", "HP:0001999 (Pancreatic insufficiency)", "HP:0000855 (Growth retardation)"],
        "inheritance": "Autosomal Recessive",
        "genes": ["CFTR"]
    },
    {
        "name": "Huntington's Disease",
        "hpo": ["HP:0002061 (Chorea)", "HP:0002511 (Progressive cognitive decline)", "HP:0000716 (Depression)"],
        "inheritance": "Autosomal Dominant",
        "genes": ["HTT"]
    },
    {
        "name": "Duchenne Muscular Dystrophy",
        "hpo": ["HP:0003735 (Muscle weakness)", "HP:0003731 (Gowers' sign)", "HP:0003197 (Calf muscle pseudohypertrophy)"],
        "inheritance": "X-linked",
        "genes": ["DMD"]
    },
    {
        "name": "Sickle Cell Anemia",
        "hpo": ["HP:0001903 (Anemia)", "HP:0001297 (Splenomegaly)", "HP:0002130 (Pain crisis)"],
        "inheritance": "Autosomal Recessive",
        "genes": ["HBB"]
    },
    {
        "name": "Marfan Syndrome",
        "hpo": ["HP:0002705 (Aortic dissection)", "HP:0001083 (Arachnodactyly)", "HP:0000501 (Ectopia lentis)"],
        "inheritance": "Autosomal Dominant",
        "genes": ["FBN1"]
    }
]

# Mock variants with different classifications for a few genes
# In a real scenario, these would come from your ClinVar/gnomAD parsed database
MOCK_VARIANTS_DB = {
    "CFTR": [
        {"chr": "chr7", "pos": "117559591", "id": "rs113993960", "ref": "G", "alt": "A", "clnsig": "Pathogenic", "freq": "0.0001"}, # F508del equivalent
        {"chr": "chr7", "pos": "117560416", "id": "rs397509172", "ref": "T", "alt": "A", "clnsig": "Pathogenic", "freq": "0.00005"}, # W1282X equivalent
        {"chr": "chr7", "pos": "117559871", "id": "rs123456789", "ref": "C", "alt": "G", "clnsig": "Likely_pathogenic", "freq": "0.00002"},
        {"chr": "chr7", "pos": "117561000", "id": "rs987654321", "ref": "A", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.001"},
        {"chr": "chr7", "pos": "117562000", "id": "rs234567890", "ref": "G", "alt": "C", "clnsig": "Benign", "freq": "0.2"},
        {"chr": "chr7", "pos": "117563000", "id": "rs345678901", "ref": "T", "alt": "C", "clnsig": "Likely_benign", "freq": "0.15"},
    ],
    "HTT": [
        {"chr": "chr4", "pos": "3076600", "id": "rsXXXXXXX", "ref": "CAGCAGCAG", "alt": "CAGCAGCAGCAG", "clnsig": "Pathogenic", "freq": "0.00001"}, # HTT expansion
        {"chr": "chr4", "pos": "3076700", "id": "rsYYYYYYY", "ref": "C", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.0005"},
        {"chr": "chr4", "pos": "3076800", "id": "rsZZZZZZZ", "ref": "A", "alt": "G", "clnsig": "Benign", "freq": "0.3"},
    ],
    "DMD": [
        {"chr": "chrX", "pos": "31234567", "id": "del_exon50", "ref": "N", "alt": "<DEL>", "clnsig": "Pathogenic", "freq": "0.00001"}, # Exon deletion
        {"chr": "chrX", "pos": "31234700", "id": "rsDMD123", "ref": "C", "alt": "A", "clnsig": "Likely_pathogenic", "freq": "0.000005"},
        {"chr": "chrX", "pos": "31235000", "id": "rsDMD456", "ref": "G", "alt": "T", "clnsig": "Uncertain_significance", "freq": "0.0001"},
    ],
    "HBB": [
        {"chr": "chr11", "pos": "5227000", "id": "rs334", "ref": "A", "alt": "T", "clnsig": "Pathogenic", "freq": "0.01"}, # SCD variant
        {"chr": "chr11", "pos": "5227500", "id": "rsHBB123", "ref": "G", "alt": "A", "clnsig": "Likely_benign", "freq": "0.1"},
    ],
    "FBN1": [
        {"chr": "chr15", "pos": "48740000", "id": "rsFBN1_path1", "ref": "C", "alt": "T", "clnsig": "Pathogenic", "freq": "0.00001"},
        {"chr": "chr15", "pos": "48740500", "id": "rsFBN1_VUS1", "ref": "A", "alt": "G", "clnsig": "Uncertain_significance", "freq": "0.00005"},
        {"chr": "chr15", "pos": "48741000", "id": "rsFBN1_benign1", "ref": "T", "alt": "C", "clnsig": "Benign", "freq": "0.05"},
    ]
}

# Mock generic background variants (not gene-specific)
MOCK_BACKGROUND_VARIANTS = [
    {"chr": f"chr{random.randint(1,22)}", "pos": str(random.randint(100000, 200000000)), "id": f"rsBG_{i}",
     "ref": random.choice(["A", "T", "C", "G"]), "alt": random.choice(["A", "T", "C", "G"]),
     "clnsig": "Benign", "freq": str(random.uniform(0.05, 0.4))} for i in range(50)
] + [
    {"chr": f"chr{random.randint(1,22)}", "pos": str(random.randint(100000, 200000000)), "id": f"rsLBBG_{i}",
     "ref": random.choice(["A", "T", "C", "G"]), "alt": random.choice(["A", "T", "C", "G"]),
     "clnsig": "Likely_benign", "freq": str(random.uniform(0.01, 0.05))} for i in range(20)
] + [
    {"chr": f"chr{random.randint(1,22)}", "pos": str(random.randint(100000, 200000000)), "id": f"rsVUSBG_{i}",
     "ref": random.choice(["A", "T", "C", "G"]), "alt": random.choice(["A", "T", "C", "G"]),
     "clnsig": "Uncertain_significance", "freq": str(random.uniform(0.0001, 0.01))} for i in range(10)
]


# --- AI Workflow Functions (MOCK IMPLEMENTATIONS) ---

def select_genetic_disorder_ai(mode="auto", user_query=None):
    """
    TODO: REPLACE WITH YOUR ACTUAL GPT-STYLE MODEL LOGIC.
    This mock function randomly selects a disorder or picks based on a query.
    """
    if mode == "auto":
        return random.choice(MOCK_GENETIC_DISORDERS)
    elif mode == "query" and user_query:
        # Simple search for demonstration
        for disorder in MOCK_GENETIC_DISORDERS:
            if user_query.lower() in disorder["name"].lower():
                return disorder
        st.warning(f"Disease '{user_query}' not found in mock database. Selecting a random one.")
        return random.choice(MOCK_GENETIC_DISORDERS)
    return random.choice(MOCK_GENETIC_DISORDERS)

def assemble_variants_mock(
    selected_disease_info, inheritance_model, num_pathogenic_variants,
    num_vus_variants, num_likely_benign_variants, num_benign_variants,
    variant_generation_mode, proband_sex
):
    """
    TODO: REPLACE WITH YOUR ACTUAL CLINVAR/GNOMAD DATABASE QUERYING AND LOGIC.
    This mock function assembles variants based on disease and desired mix.
    """
    assembled_variants = []
    disease_genes = selected_disease_info.get("genes", [])
    all_gene_variants = []
    for gene in disease_genes:
        all_gene_variants.extend(MOCK_VARIANTS_DB.get(gene, []))

    # 1. Add Pathogenic/Likely Pathogenic variants for the selected disease
    pathogenic_candidates = [v for v in all_gene_variants if v["clnsig"] in ["Pathogenic", "Likely_pathogenic"]]

    if inheritance_model == "Autosomal Recessive":
        # Need two pathogenic variants for AR
        if len(pathogenic_candidates) >= 2:
            chosen = random.sample(pathogenic_candidates, 2)
            assembled_variants.extend(chosen)
            # Assign genotypes based on inheritance model (heterozygous for both)
            chosen[0]['gt'] = '0/1'
            chosen[1]['gt'] = '0/1'
        elif len(pathogenic_candidates) == 1:
            # If only one found, duplicate it to simulate compound heterozygous
            chosen = random.sample(pathogenic_candidates, 1)
            assembled_variants.extend(chosen)
            assembled_variants.append(pathogenic_candidates[0])
            chosen[0]['gt'] = '0/1'
            assembled_variants[-1]['gt'] = '0/1'
            st.warning("Not enough distinct pathogenic variants for AR. Simulating compound heterozygous with available.")
        else:
            st.error("No pathogenic variants found for selected AR disease in mock DB.")
            return []
    elif inheritance_model == "Autosomal Dominant":
        # Need at least one pathogenic variant for AD
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            assembled_variants.append(chosen)
            chosen['gt'] = '0/1'
        else:
            st.error("No pathogenic variants found for selected AD disease in mock DB.")
            return []
    elif inheritance_model == "X-linked":
        # For X-linked, male (XY) gets 1, female (XX) can be carrier or affected (2)
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            assembled_variants.append(chosen)
            if proband_sex == "Male":
                chosen['gt'] = '1/1' # Hemizygous
            else: # Female
                chosen['gt'] = '0/1' # Heterozygous (carrier or affected depending on gene dose)
        else:
            st.error("No pathogenic variants found for selected X-linked disease in mock DB.")
            return []
    else: # Other inheritance, default to one pathogenic if any
        if pathogenic_candidates:
            chosen = random.choice(pathogenic_candidates)
            assembled_variants.append(chosen)
            chosen['gt'] = '0/1'
        else:
            st.warning("No specific pathogenic variant logic for this inheritance model, or no pathogenic variants found.")


    # 2. Add other classifications (VUS, LB, B) from background variants
    random.shuffle(MOCK_BACKGROUND_VARIANTS)

    # Filter out any background variants that might conflict with primary disease gene variants
    # Simple check: avoid variants at same position (for this mock)
    existing_positions = {(v['chr'], v['pos']) for v in assembled_variants}
    filtered_background = [
        v for v in MOCK_BACKGROUND_VARIANTS
        if (v['chr'], v['pos']) not in existing_positions
    ]

    def add_variants(count, classification, source_list):
        added = 0
        for variant in source_list:
            if added < count and variant["clnsig"] == classification:
                # Ensure different positions for unique variants
                if (variant['chr'], variant['pos']) not in existing_positions:
                    # Assign a random genotype, mostly heterozygous for background variants
                    # For homozygous benign, add 0/0. For heterozygous, add 0/1 or 1/0
                    variant['gt'] = random.choice(['0/1', '1/0', '0/0'])
                    assembled_variants.append(variant)
                    existing_positions.add((variant['chr'], variant['pos']))
                    added += 1
        return added

    add_variants(num_vus_variants, "Uncertain_significance", filtered_background)
    add_variants(num_likely_benign_variants, "Likely_benign", filtered_background)
    add_variants(num_benign_variants, "Benign", filtered_background)

    # Ensure unique variants (in case of overlaps in mock data or sampling)
    final_variants = []
    seen_ids = set()
    for var in assembled_variants:
        if var['id'] not in seen_ids:
            final_variants.append(var)
            seen_ids.add(var['id'])
        else: # If same ID, ensure position is different or modify ID
            temp_id = var['id']
            count = 1
            while f"{temp_id}_{count}" in seen_ids:
                count += 1
            var['id'] = f"{temp_id}_{count}"
            final_variants.append(var)
            seen_ids.add(var['id'])


    # Simulate de novo/synthetic: slightly randomize positions if chosen
    if variant_generation_mode == "De Novo/Synthetic (randomized)":
        for var in final_variants:
            if var["clnsig"] not in ["Pathogenic", "Likely_pathogenic"]: # Don't randomize core disease variants
                try:
                    original_pos = int(var['pos'])
                    var['pos'] = str(max(1, original_pos + random.randint(-500, 500))) # +- 500 bp
                except ValueError:
                    pass # Keep original if not convertible to int

    return final_variants

def generate_patient_metadata(sex, age, proband_id_base, ethnicity, generate_ped):
    """
    Generates patient metadata.
    """
    proband_id = proband_id_base if proband_id_base else f"SYNTH_PATIENT_{str(uuid.uuid4())[:8].upper()}"
    patient_metadata = {
        "Sex": sex,
        "Age": age,
        "Proband ID": proband_id,
        "Ethnicity": ethnicity,
        "Family Structure Generated": generate_ped
    }
    return patient_metadata

def generate_ped_file_content(proband_id, sex, disease_name):
    """
    Generates simple PED file content (proband, mother, father).
    Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male, 2=female), Phenotype (2=affected)
    """
    family_id = f"FAM_{proband_id}"
    ped_content = [
        f"{family_id}\t{proband_id}\t{proband_id}_father\t{proband_id}_mother\t{ (1 if sex == 'Male' else 2) }\t2 # Affected by {disease_name}",
        f"{family_id}\t{proband_id}_father\t0\t0\t1\t1 # Unaffected",
        f"{family_id}\t{proband_id}_mother\t0\t0\t2\t1 # Unaffected"
    ]
    return "\n".join(ped_content)

def generate_vcf_file_content(variants, patient_metadata):
    """
    Generates VCF file content.
    """
    proband_id = patient_metadata["Proband ID"]
    file_date = datetime.now().strftime("%Y%m%d")

    # VCF Header
    header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={file_date}",
        "##source=SyntheticPatientCaseGenerator",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth (mock value)\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency (mock value or actual frequency)\">",
        "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance of the variant. Possible values: Benign, Likely_benign, Uncertain_significance, Likely_pathogenic, Pathogenic, drug_response, histocompatibility, other, protective, risk_factor, Conflicting_interpretations_of_pathogenicity\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{proband_id}"
    ]

    # VCF Records
    records = []
    for var in variants:
        chrom = var.get("chr", ".")
        pos = var.get("pos", ".")
        id_val = var.get("id", ".")
        ref = var.get("ref", ".")
        alt = var.get("alt", ".")
        qual = "."
        filter_val = "PASS" # Assuming all generated variants pass filters for simplicity
        clnsig = var.get("clnsig", "Uncertain_significance").replace(" ", "_") # VCF requires no spaces
        af = var.get("freq", "0.01") # Use mock frequency
        info = f"DP=100;AF={af};CLNSIG={clnsig}" # Mock DP
        gt = var.get("gt", "0/1") # Genotype from assembly logic

        records.append(f"{chrom}\t{pos}\t{id_val}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\tGT\t{gt}")

    return "\n".join(header + records)


# --- Streamlit Application ---

st.set_page_config(layout="wide", page_title="Synthetic Patient Case Generator")

st.title("🧬 Synthetic Patient Case Generator")
st.markdown("""
    This application automatically generates realistic synthetic patient cases,
    complete with HPO terms, a mix of pathogenic and benign variants,
    and properly formatted VCF and optional PED files.
    **Note:** This application uses mock data for demonstration.
    Real-world usage requires integration with biological databases (ClinVar, gnomAD etc.)
    and advanced AI models.
""")

# --- Sidebar for Global Settings/Information ---
st.sidebar.header("Settings & Information")
st.sidebar.info("This tool is for research and educational purposes only and does not generate real patient data.")
st.sidebar.markdown("---")
st.sidebar.subheader("Contact & Feedback")
st.sidebar.markdown("For questions or suggestions, please reach out.")


# --- Step 1: Disease & Phenotype Selection ---
st.header("Step 1: Disease & Phenotype Selection")
st.markdown("Specify the criteria for the genetic disorder you want to simulate.")

disease_selection_method = st.radio(
    "How do you want to select the disease?",
    ("Automatic (AI-driven)", "Manual Selection")
)

selected_disease_info = None

if disease_selection_method == "Automatic (AI-driven)":
    st.info("The AI model will select a realistic genetic disorder based on internal criteria (mock selection).")
    if st.button("Generate Automatic Disease Suggestion", key="auto_disease_btn"):
        with st.spinner("AI is selecting a disease..."):
            selected_disease_info = select_genetic_disorder_ai(mode="auto")
            st.session_state['selected_disease_info'] = selected_disease_info
            st.success("Disease suggested!")

        if selected_disease_info:
            st.write(f"**Suggested Disease:** {selected_disease_info['name']}")
            st.write(f"**Associated HPO Terms:** {', '.join(selected_disease_info['hpo'])}")
            st.write(f"**Inheritance Model:** {selected_disease_info['inheritance']}")

elif disease_selection_method == "Manual Selection":
    # Allow user to pick from mock list for now
    disease_names = [d["name"] for d in MOCK_GENETIC_DISORDERS]
    selected_disease_name = st.selectbox("Select a genetic disorder:", [""] + disease_names)

    if selected_disease_name:
        selected_disease_info = next((d for d in MOCK_GENETIC_DISORDERS if d["name"] == selected_disease_name), None)
        st.session_state['selected_disease_info'] = selected_disease_info
        if selected_disease_info:
            st.write(f"**Selected Disease:** {selected_disease_info['name']}")
            st.write(f"**Associated HPO Terms:** {', '.join(selected_disease_info['hpo'])}")
            st.write(f"**Inheritance Model:** {selected_disease_info['inheritance']}")
        else:
            st.warning("Please select a disease from the dropdown.")
    else:
        st.info("Please select a disease from the list above.")

# --- Step 2: Variant Assembly ---
st.header("Step 2: Variant Assembly")
st.markdown("Based on the selected disease, variants will be assembled, including pathogenic and background variants.")

if 'selected_disease_info' in st.session_state and st.session_state['selected_disease_info']:
    disease_name_for_variants = st.session_state['selected_disease_info']['name']
    inheritance_model_for_variants = st.session_state['selected_disease_info']['inheritance']
    st.write(f"**Disease for Variant Assembly:** {disease_name_for_variants} (Inheritance: {inheritance_model_for_variants})")

    st.subheader("Variant Classification Distribution:")
    st.markdown("Specify the number of variants for each pathogenicity classification.")

    num_pathogenic = 0
    if inheritance_model_for_variants == "Autosomal Recessive":
        num_pathogenic = st.slider("Number of Pathogenic/Likely Pathogenic variants (for AR, typically 2):", 1, 3, 2)
    elif inheritance_model_for_variants == "Autosomal Dominant" or inheritance_model_for_variants == "X-linked":
        num_pathogenic = st.slider("Number of Pathogenic/Likely Pathogenic variants (for AD/XL, typically 1):", 1, 3, 1)
    else:
        num_pathogenic = st.slider("Number of Pathogenic/Likely Pathogenic variants:", 1, 3, 1)

    num_vus = st.slider("Number of Uncertain Significance (VUS) variants:", 0, 10, 3)
    num_likely_benign = st.slider("Number of Likely Benign variants:", 0, 20, 5)
    num_benign = st.slider("Number of Benign variants:", 0, 30, 10)


    variant_generation_mode = st.radio(
        "Variant Generation Mode:",
        ("Realistic (from mock DB)", "De Novo/Synthetic (randomized positions)")
    )

    if st.button("Assemble Variants", key="assemble_variants_btn"):
        with st.spinner("Assembling variants... (This would query ClinVar/gnomAD in a real app)."):
            # Need to get proband sex for X-linked
            proband_sex_for_variants = st.session_state.get('patient_metadata', {}).get('Sex', 'Male') # Default to Male if not set yet

            assembled_variants = assemble_variants_mock(
                st.session_state['selected_disease_info'],
                inheritance_model_for_variants,
                num_pathogenic,
                num_vus,
                num_likely_benign,
                num_benign,
                variant_generation_mode,
                proband_sex_for_variants
            )
            st.session_state['assembled_variants'] = assembled_variants
            st.success(f"Assembled {len(assembled_variants)} variants!")

            st.subheader("Assembled Variants Overview:")
            if assembled_variants:
                df_variants = pd.DataFrame(assembled_variants)
                df_variants_display = df_variants[['chr', 'pos', 'id', 'ref', 'alt', 'clnsig', 'freq', 'gt']]
                st.dataframe(df_variants_display)
                st.markdown(f"**Total Variants Assembled:** {len(assembled_variants)}")
                st.markdown("**(Note: `gt` is the genotype assigned for the synthetic proband)**")
            else:
                st.warning("No variants were assembled. Please check your selections or mock data.")
else:
    st.info("Please complete 'Step 1: Disease & Phenotype Selection' first.")


# --- Step 3: Patient Metadata ---
st.header("Step 3: Patient Metadata Generation")
st.markdown("Generate demographic and family information for the synthetic patient.")

if 'assembled_variants' in st.session_state and st.session_state['assembled_variants']:
    sex_options = ["Male", "Female", "Other"]
    sex = st.selectbox("Sex:", sex_options, key="proband_sex_input")
    age = st.slider("Age (Years):", 0, 100, 5, key="proband_age_input")
    proband_id_base = st.text_input("Proband ID (leave blank for auto-generate):", value="", key="proband_id_input")
    ethnicity = st.selectbox(
        "Ethnicity (optional):",
        ["Not Specified", "European", "Asian", "African", "Hispanic/Latino", "Middle Eastern", "South Asian", "Other"],
        key="proband_ethnicity_input"
    )
    generate_ped = st.checkbox("Generate Family Structure (PED file)", value=True, key="generate_ped_checkbox")

    if st.button("Generate Patient Metadata", key="generate_metadata_btn"):
        with st.spinner("Generating patient metadata..."):
            patient_metadata = generate_patient_metadata(sex, age, proband_id_base, ethnicity, generate_ped)
            st.session_state['patient_metadata'] = patient_metadata
            st.success("Patient metadata generated!")

            st.subheader("Generated Patient Metadata:")
            for key, value in patient_metadata.items():
                st.write(f"**{key}:** {value}")

            if generate_ped:
                # Generate PED file content here as it depends on proband_id and sex
                disease_name_for_ped = st.session_state['selected_disease_info']['name'] if 'selected_disease_info' in st.session_state else "Unknown Disease"
                ped_content = generate_ped_file_content(patient_metadata["Proband ID"], patient_metadata["Sex"], disease_name_for_ped)
                st.session_state['ped_content'] = ped_content
                st.subheader("Generated PED File Content (Preview):")
                st.code(ped_content, language="text")

else:
    st.info("Please complete 'Step 2: Variant Assembly' first.")


# --- Step 4: VCF Generator ---
st.header("Step 4: VCF File Generation")
st.markdown("Generate the final VCF file containing the synthetic patient's variants.")

if 'assembled_variants' in st.session_state and 'patient_metadata' in st.session_state:
    if st.button("Generate VCF File", key="generate_vcf_btn"):
        with st.spinner("Generating VCF file..."):
            vcf_content = generate_vcf_file_content(
                st.session_state['assembled_variants'],
                st.session_state['patient_metadata']
            )
            st.session_state['vcf_content'] = vcf_content
            st.success("VCF file generated!")

            st.subheader("Generated VCF Content:")
            st.code(vcf_content, language="vcf")

            # Provide download button for VCF
            st.download_button(
                label="Download VCF File",
                data=vcf_content,
                file_name=f"{st.session_state['patient_metadata']['Proband ID']}.vcf",
                mime="text/vcf"
            )

            # Provide download button for PED if generated
            if st.session_state.get('ped_content'):
                st.download_button(
                    label="Download PED File",
                    data=st.session_state['ped_content'],
                    file_name=f"{st.session_state['patient_metadata']['Proband ID']}.ped",
                    mime="text/plain"
                )
else:
    st.info("Please complete 'Step 3: Patient Metadata Generation' first.")

st.markdown("---")
st.info("Developed with Streamlit for synthetic patient case generation.")
