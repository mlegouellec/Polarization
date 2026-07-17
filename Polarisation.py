
import pysam
import argparse
import sys
from tqdm import tqdm

def load_sample_list(file_path):
    try:
        with open(file_path, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        sys.exit(1)

def get_group_data(record, sample_ids):
    alleles = set()
    called_count = 0
    for sid in sample_ids:
        if sid not in record.samples: continue
        gt = record.samples[sid]['GT']
        if gt is not None and any(a is not None for a in gt):
            called_count += 1
            alleles.update(a for a in gt if a is not None)
    return alleles, called_count

def polarize_logic(record, pop_ids, out1_ids, out2_ids):
    pop_all, pop_called = get_group_data(record, pop_ids)
    out1_all, out1_called = get_group_data(record, out1_ids)
    out2_all, out2_called = get_group_data(record, out2_ids)

    # Filter for Missing Individuals (1/3 Targeted Pop, 2/3 of each Outgroup)
    if pop_called < (len(pop_ids) / 3):
        return "MISSING_DATA_POP", "NOPOL", None
    if out1_called < (2 * len(out1_ids) / 3):
        return "MISSING_DATA_OUT1", "NOPOL", None
    if out2_called < (2 * len(out2_ids) / 3):
        return "MISSING_DATA_OUT2", "NOPOL", None

    if not pop_all or not out1_all or not out2_all:
        return "CASE_MISSING", "NOPOL", None

    combined = pop_all.union(out1_all)

    if len(out2_all) > 1:
        if len(pop_all) == 1 and pop_all == out1_all:
            shared = list(pop_all)[0]
            if shared in out2_all:
                return "CASE4", "POL", shared
        return "CASE5", "NOPOL", None
    else:
        anc_cand = list(out2_all)[0]
        if len(combined) == 2:
            return "CASE1", "POL", anc_cand
        else:
            if len(combined) == 1 and combined == out2_all:
                return "CASE2", "POL", anc_cand
            if len(pop_all) == 1 and len(out1_all) == 1 and combined != out2_all:
                return "CASE3", "NOPOL", None
    return "OTHER", "NOPOL", None

def run_pipeline(input_path, out_vcf_path, out_tsv_path, pop_file, out1_file, out2_file):
    pop_list  = load_sample_list(pop_file)
    out1_list = load_sample_list(out1_file)
    out2_list = load_sample_list(out2_file)

    vcf_in = pysam.VariantFile(input_path)
    vcf_in.header.info.add("POL_STATUS", 1, "String", "Polarization status")
    vcf_in.header.info.add("POL_CASE", 1, "String", "Logic case used")

    if not out_vcf_path.endswith(".gz"): out_vcf_path += ".gz"
    vcf_out = pysam.VariantFile(out_vcf_path, 'wz', header=vcf_in.header)

    tsv_out = open(out_tsv_path, 'w')
    tsv_out.write("CHROM\tPOS\tREF_ORIG\tALT_ORIG\tANC_ALLELE\tDER_ALLELE\tCASE\tSTATUS\n")

    print(f"Processing VCF with full field :")

    try:
        for record in tqdm(vcf_in, unit=" variants"):
            ref_orig = record.ref
            alt_orig = record.alts[0] if (record.alts and len(record.alts) > 0) else "."

            case, action, anc_idx = polarize_logic(record, pop_list, out1_list, out2_list)

            status = "NOPOL"
            anc_str, der_str = ".", "."

            # First : copy all the fields from the VCF (INFO, FORMAT, etc.)
            new_record = record.copy()

            if action == "POL":
                if anc_idx == 0:
                    status = "POL1"
                    anc_str, der_str = ref_orig, alt_orig

                elif anc_idx == 1 and alt_orig != ".":
                    status = "POL2"
                    anc_str, der_str = alt_orig, ref_orig

                    # 1. Save opposite genotypes before changing the REF/ALT columns (avoid that pysam loses the FORMAT info)
                    flipped_genotypes = {}
                    for sample in record.samples:
                        gt = record.samples[sample]['GT']
                        if gt is not None:
                            flipped_genotypes[sample] = tuple((1 - a if a is not None else None) for a in gt)

                    # 2. Change REF/ALT alleles
                    new_record.ref = anc_str
                    new_record.alts = (der_str,)

                    # 3. Add flipped genotypes
                    for sample, new_gt in flipped_genotypes.items():
                        new_record.samples[sample]['GT'] = new_gt

                elif anc_idx == 0 and alt_orig == ".":
                    status = "POL1"
                    anc_str = ref_orig

            # Add polarisation taggs
            new_record.info["POL_STATUS"] = status
            new_record.info["POL_CASE"] = case

            vcf_out.write(new_record)
            tsv_out.write(f"{record.chrom}\t{record.pos}\t{ref_orig}\t{alt_orig}\t{anc_str}\t{der_str}\t{case}\t{status}\n")

    finally:
        vcf_in.close()
        vcf_out.close()
        tsv_out.close()

    print(f"\nDone. Output VCF: {out_vcf_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output_vcf", required=True)
    parser.add_argument("-t", "--output_tsv", required=True)
    parser.add_argument("--pop", required=True)
    parser.add_argument("--out1", required=True)
    parser.add_argument("--out2", required=True)
    args = parser.parse_args()
    run_pipeline(args.input, args.output_vcf, args.output_tsv, args.pop, args.out1, args.out2)
