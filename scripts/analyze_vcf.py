import pandas as 
import gzip
import sys

def parse_vcf_detailed(vcf_path):
    print(f"Parsing {vcf_path}...")
    variants = []
    
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            cols = line.strip().split("\t")
            
            # basic variant info
            variant_info = {
                'chrom' : cols[0],
                'pos' : cols[1],
                'ref' : cols[3],
                'alt' : cols[4],
                'qual' : cols[5],
                'filter_status' : cols[6],
                'info_str' : cols[7],
                'format_str' : cols[8],
                'tumor_sample' : cols[9] # assumes tumor is first sample col
                'normal_sample' : cols[10] # assumes normal is second sample col
            }
            
            # parse INFO column 
            info_dict = {x.split("=")[0]: x.split("=")[1] for x in info_str.split(";") if "=" in x}
            gene = "Intergenic"
            impact = "MODIFIER"
            consequence = "unknown"
            protein_change = "."
            
            if "ANN" in info_dict:
                ann_list = info_dict["ANN"].split(",")
                first_ann = ann_list[0].split("|")
                if len(first_ann) > 10:
                    consequence = first_ann[1]
                    impact = first_ann[2]
                    gene = first_ann[3]
                    protein_change = first_ann[10] # HGVS.p

            # parse FORMAT column
            fmt_keys = format_str.split(":")
            tumor_vals = variant_info['tumor_sample'].split(":")
            normal_vals = variant_info['normal_sample'].split(":")
            
            fmt_dict_tumor = dict(zip(fmt_keys, tumor_vals))
            fmt_dict_normal = dict(zip(fmt_keys, normal_vals))

            # extract allele depth (AD) -> ref_reads, alt_reads
            try:
                t_ref_count, t_alt_count = map(int, fmt_dict_tumor.get("AD", "0,0").split(","))
                n_ref_count, n_alt_count = map(int, fmt_dict_normal.get("AD", "0,0").split(","))
                
                # calculate variant allele frequency (VAF)
                t_vaf = t_alt_count / (t_ref_count + t_alt_count) if (t_ref_count + t_alt_count) > 0 else 0
                n_vaf = n_alt_count / (n_ref_count + n_alt_count) if (n_ref_count + n_alt_count) > 0 else 0
            except ValueError:
                t_vaf = 0
                n_vaf = 0

            variants.append({
                "Chrom": variant_info['chrom'],
                "Pos": variant_info['pos'],
                "Ref": variant_info['ref'],
                "Alt": variant_info['alt'],
                "Gene": gene,
                "Consequence": consequence,
                "Impact": impact,
                "Protein_Change": protein_change,
                "Tumor_VAF": round(t_vaf, 3),
                "Normal_VAF": round(n_vaf, 3),
                "Depth_Tumor": t_ref_count + t_alt_count,
                "Filter": variant_info['filter_status']
            })
            
    return pd.DataFrame(variants)

# ----------
vcf_file = "results/vcf/somatic/annotated/Tumor_1000G.somatic.ann.vcf.gz"
df = parse_vcf_detailed(vcf_file)

pass_df = df[df["Filter"] == "PASS"]

print(f"Total Variants Detected: {len(df)}")
print(f"High Quality (PASS) Variants: {len(pass_df)}")

high_impact = pass_df[pass_df["Impact"].isin(["HIGH", "MODERATE"])]
print("\n--- Top High Impact Mutations ---")
print(high_impact[["Gene", "Consequence", "Protein_Change", "Tumor_VAF"]].head(10).to_string(index=False))

pass_df.to_csv("results/clinical_report.csv", index=False)
print("\nFull clinical report saved to: results/clinical_report.csv")
