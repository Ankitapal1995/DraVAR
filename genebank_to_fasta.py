from Bio import SeqIO
from sys import argv
s,f=argv
gbk_filename=f
extension_to_remove=['gb','gbk']
faa_filename=[f[:f.find(ext)] for ext in extension_to_remove if f.endswith(ext)][0]+"fasta"
print(faa_filename)
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print ("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
