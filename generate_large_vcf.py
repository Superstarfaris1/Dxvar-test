"""Command line tool to generate a synthetic patient case with a large VCF file."""

import argparse
from typing import Optional

import case_lib
from case_lib import (
    select_genetic_disorder_auto,
    assemble_variants,
    generate_patient_metadata,
    generate_ped_file_content,
    generate_vcf_file_content,
)


def main(
    disorder_name: Optional[str],
    pathogenic: int,
    vus: int,
    likely_benign: int,
    benign: int,
    output_prefix: str,
    sex: str,
    age: int,
    ethnicity: str,
    generation_mode: str,
):
    if disorder_name:
        disorder = next(
            (
                d
                for d in case_lib.MOCK_GENETIC_DISORDERS
                if d["name"].lower() == disorder_name.lower()
            ),
            None,
        )
        if disorder is None:
            print(f"Disease '{disorder_name}' not found. Using random disorder.")
            disorder = select_genetic_disorder_auto()
    else:
        disorder = select_genetic_disorder_auto()

    variants = assemble_variants(
        disorder,
        disorder["inheritance"],
        pathogenic,
        vus,
        likely_benign,
        benign,
        generation_mode,
        sex,
    )

    metadata = generate_patient_metadata(sex, age, output_prefix, ethnicity, True)
    ped_content = generate_ped_file_content(metadata["Proband ID"], sex, disorder["name"])
    vcf_content = generate_vcf_file_content(variants, metadata)

    vcf_file = f"{output_prefix}.vcf"
    with open(vcf_file, "w") as f:
        f.write(vcf_content)

    ped_file = f"{output_prefix}.ped"
    with open(ped_file, "w") as f:
        f.write(ped_content)

    print(f"Generated VCF: {vcf_file}")
    print(f"Generated PED: {ped_file}")
    print("Associated HPO terms:")
    print("\n".join(disorder["hpo"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a synthetic case")
    parser.add_argument("--disease", help="Name of disease to simulate")
    parser.add_argument("--pathogenic", type=int, default=2, help="Number of pathogenic variants")
    parser.add_argument("--vus", type=int, default=10, help="Number of VUS variants")
    parser.add_argument("--likely-benign", type=int, default=20, help="Number of likely benign variants")
    parser.add_argument("--benign", type=int, default=50, help="Number of benign variants")
    parser.add_argument("--output-prefix", default="synthetic_case", help="Prefix for output files")
    parser.add_argument("--sex", default="Male", help="Sex of proband")
    parser.add_argument("--age", type=int, default=5, help="Age of proband")
    parser.add_argument("--ethnicity", default="Not Specified", help="Ethnicity")
    parser.add_argument(
        "--generation-mode",
        default="De Novo/Synthetic (randomized)",
        choices=["De Novo/Synthetic (randomized)", "Realistic (from mock DB)"],
    )

    args = parser.parse_args()

    main(
        args.disease,
        args.pathogenic,
        args.vus,
        args.likely_benign,
        args.benign,
        args.output_prefix,
        args.sex,
        args.age,
        args.ethnicity,
        args.generation_mode,
    )
