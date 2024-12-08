import vcfpy
import csv
from sys import argv
s,f=argv
# Open the annotated VCF file
vcf_reader = vcfpy.Reader.from_path(f)

# Open CSV for writing
with open(f.split("_")[0]+"_annotated.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(["CHROM", "POS", "REF", "ALT","QUAL", "GENE", "IMPACT", "EFFECT", "HGVS_C", "HGVS_P"])
    
    # Parse each record
    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alts = ",".join(str(alt) for alt in record.ALT)
        qual = record.QUAL
        
        # Extract annotations from INFO/ANN
        annotations = record.INFO.get("ANN", [])
        for ann in annotations:
            fields = ann.split("|")
            gene = fields[3]
            impact = fields[2]
            effect = fields[1]
            hgvs_c = fields[9]
            hgvs_p = fields[10]
            
            # Write row
            writer.writerow([chrom, pos, ref, alts, gene, impact, effect, hgvs_c, hgvs_p])

