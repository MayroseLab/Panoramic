"""
Given a Panoramic pipeline logs dir,
calculate the total time the run took.
This is useful for discontinuous runs.
Optionally, can take a comma-separated
list of rules to work on. Otherwise will
calculate on all rules.
"""

from datetime import datetime
import os
import calendar
month_trans = {v: k for k,v in enumerate(calendar.month_abbr)}
from intervaltree import Interval, IntervalTree
import sys


logs_dir = sys.argv[1]
if len(sys.argv) >= 3:
  rules = sys.argv[2].split(',')
else:
  rules = ['']

def convert_time(time_str):
  """
  Converts Sat Aug 15 22:56:14 IDT 2020
  to time since epoch
  """
  week_day, month_name, day, time, timezone, year = time_str.split()
  year = int(year)
  month = month_trans[month_name]
  hour, minute, second = [int(x) for x in time.split(':')]
  day = int(day)
  return datetime(year,month,day,hour,minute,second).timestamp()

t = IntervalTree()
for log_file in os.listdir(logs_dir):
  if any([log_file.endswith(rule + '.out') for rule in rules]):
    with open(os.path.join(logs_dir,log_file)) as f:
      for line in f:
        if line.startswith('Start time:'):
          start_time = convert_time(line[12:-1])
        elif line.startswith('End time:'):
          end_time = convert_time(line[10:-1])
      if start_time == end_time:
        end_time += 1
      t[start_time:end_time] = log_file

t.merge_overlaps()
tot = 0
for inter in t:
  tot += inter.end - inter.begin
tot_hours = tot/(60*60)
print("Total run time: %s hours" % tot_hours)
