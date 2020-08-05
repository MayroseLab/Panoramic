in_gff = "GCA_000149365.1_ASM14936v1_genomic.gff"
in_prot = "GCA_000149365.1_ASM14936v1_protein.faa"

replace = {}
with open(in_gff) as f:
  for line in f:
    if line.startswith('#'):
      continue
    fields = line.strip().split('\t')
    if fields[2] == "CDS":
      attr = {k.split('=')[0]: k.split('=')[1] for k in line.split('\t')[8].split(';')}
      replace[attr["Name"]] = attr["Parent"]
with open(in_prot) as f:
  for line in f:
    line = line.strip()
    if line.startswith('>'):
      name = line.split()[0][1:]
      print(">%s" % replace[name])
    else:
      print(line)

