"""
Allows easy editing of MAKER configurations
"""
from __future__ import print_function
import argparse

def pair(arg):
  return arg.split('=')

parser = argparse.ArgumentParser(description='Edit MAKER config')
parser.add_argument('in_conf', help="Input MAKER configuration")
parser.add_argument('out_conf', help="Output MAKER configuration (after edits)")
parser.add_argument('--edits', type=pair, nargs='+')
args = parser.parse_args()

edits_dict = {p[0]: p[1] for p in args.edits}
with open(args.in_conf) as f, open(args.out_conf,'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('#') or line == '':
      print(line, file=fo)
      continue
    key, val = line.split('=', 1)
    if key in edits_dict:
      print("%s=%s" %(key, edits_dict[key]), file=fo)
    else:
      print(line, file=fo)
