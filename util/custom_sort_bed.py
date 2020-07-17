"""
Allows sorting a bed file
by a custom chromosomes order.
Within each chromosome, the
bed features are sorted by
start position
"""

from __future__ import print_function
import sys

in_bed = sys.argv[1]
order_file = sys.argv[2]

# read order file and create a dict with indices
order_dict = {}
with open(order_file) as f:
  i = 0
  for line in f:
    order_dict[line.strip()] = i
    i += 1

# read in bed file and save as list of lists
bed_list = []
with open(in_bed) as f:
  for line in f:
    line = line.strip()
    fields = line.split('\t')
    if fields[0] not in order_dict:
      order_dict[fields[0]] = i
      i += 1
    feature_list = [fields[0], int(fields[1]), line]
    bed_list.append(feature_list)

# sort list
def sort_bed_list(l):
  return order_dict[l[0]], l[1]

sorted_list = sorted(bed_list, key = sort_bed_list)
for b in sorted_list:
  print(b[2])
