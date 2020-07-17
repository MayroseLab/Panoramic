"""
If for any reason the proteins fasta output
of MAKER contains duplicate names, change
names to prevent this issue.
"""

from __future__ import print_function
import sys

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

names = {}
with open(in_fasta) as f, open(out_fasta, 'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('>'):
      name_parts = line[1:].split(' ', 1)
      name = name_parts[0]
      rest = ''
      if len(name_parts) == 2:
        rest = name_parts[1]
      if name in names:
        names[name] += 1
        name = "%s_%s" %(name, names[name])
      else:
        names[name] = 1
      print(">%s %s" %(name, rest), file=fo)
    else:
      print(line, file=fo)
