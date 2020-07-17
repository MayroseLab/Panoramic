"""
Simplify the ID and Parent attributes
of a gff file using a function specified
by the user.
"""

from __future__ import print_function
import sys

in_gff = sys.argv[1]
out_gff = sys.argv[2]
function = sys.argv[3]	# use lambda syntax

simp_function = eval(function)
with open(in_gff) as f, open(out_gff,'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      print(line, file=fo)
      continue
    fields = line.split('\t')
    attr = {k.split('=')[0]: k.split('=')[1] for k in fields[8].split(';')}
    if 'ID' in attr:
      attr['ID'] = simp_function(attr['ID'])
    if 'Parent' in attr:
      attr['Parent'] = simp_function(attr['Parent'])
    simp_attr_lst = ["%s=%s" %(key, attr[key]) for key in attr]
    simp_attr_str = ';'.join(simp_attr_lst)
    fields[8] = simp_attr_str
    print('\t'.join(fields), file=fo)
