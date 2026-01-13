import pandas as pd
import gzip
import sys

def parse_vcf(vcf_path):
    print(f"parsing {vcf_path}...")
    variants = []
    csq_fields = {}
    t_idx = -1
    n_idx = -1
    
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                # you can dynamically pass in sample names here if you want, otherwise skip the CHROM parsing and just pass in the correct column indexes for your tumor sample and normal sample
                try:
                    t_idx = header.index("Tumor_1000G")
                    n_idx = header.index("Normal_1000G")
                    print(f"detected sample columns: tumor at {t_idx}, normal at {n_idx}")
                except ValueError:
                    print("error: could not find sample IDs 'Tumor_1000G' or 'Normal_1000G' in header.")
                    print(f"header: {header}")
                    sys.exit(1)
                continue

            if line.startswith("##INFO=<ID=CSQ"):
                try:
                    format_str = line.split("Format:")[1].strip().strip('">').strip("'")
                    headers = format_str.split("|")
                    csq_fields = {name: i for i, name in enumerate(headers)}
                except:
                    pass
                continue

            if line.startswith("#"):
                continue
            
            cols = line.strip().split("\t")
            
            info_str = cols[7]
            info_dict = {}
            for x in info_str.split(";"):
                if "=" in x:
                    key, val = x.split("=", 1)
                    info_dict[key] = val
            
            gene = "."
            impact = "MODIFIER"
            consequence = "."
            protein_change = "."
            
            tag = "CSQ" if "CSQ" in info_dict else ("ANN" if "ANN" in info_dict else None)
            
            if tag:
                ann_list = info_dict[tag].split(",")
                best_rank = 0
                rank_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
                
                for ann in ann_list:
                    fields = ann.split("|")
                    try:
                        curr_impact = fields[csq_fields.get('IMPACT', 2)]
                        curr_cons = fields[csq_fields.get('Consequence', 1)]
                        curr_gene = fields[csq_fields.get('SYMBOL', 3)]
                        curr_prot = fields[csq_fields.get('HGVSp', 10)] if 'HGVSp' in csq_fields else "."
                    except IndexError:
                        continue 

                    rank = rank_map.get(curr_impact, 0)
                    if rank > best_rank:
                        best_rank = rank
                        impact = curr_impact
                        consequence = curr_cons
                        gene = curr_gene
                        protein_change = curr_prot

            format_str = cols[8]
            tumor_sample = cols[t_idx]
            normal_sample = cols[n_idx]
            
            fmt_keys = format_str.split(":")
            fmt_tumor = dict(zip(fmt_keys, tumor_sample.split(":")))
            fmt_normal = dict(zip(fmt_keys, normal_sample.split(":")))

            try:
                t_ref, t_alt = map(int, fmt_tumor.get("AD", "0,0").split(","))
                t_depth = t_ref + t_alt
                t_vaf = t_alt / t_depth if t_depth > 0 else 0
            except ValueError:
                t_vaf = 0

            variants.append({
                "Gene": gene,
                "Impact": impact,
                "Consequence": consequence,
                "Protein_Change": protein_change,
                "Tumor_VAF": round(t_vaf, 3),
                "Filter": cols[6]
            })
            
    return pd.DataFrame(variants)

# ----------
vcf_file = "results/vcf/somatic/annotated/Tumor_1000G.somatic.ann.vcf.gz"
df = parse_vcf(vcf_file)

pass_df = df[df["Filter"] == "PASS"]
high = pass_df[pass_df["Impact"].isin(["HIGH", "MODERATE"])]

print(f"\n--- Top High Impact Mutations (Count: {len(high)}) ---")
if not high.empty:
    print(high[["Gene", "Consequence", "Protein_Change", "Tumor_VAF"]]
          .sort_values("Tumor_VAF", ascending=False)
          .head(15)
          .to_string(index=False))

pass_df.to_csv("results/clinical_report.csv", index=False)
