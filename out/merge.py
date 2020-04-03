import fileinput
import ROOT
import re
import numpy
import multiprocessing
import sys
import os

def ReadOutFile(fn):
  match = re.match('(\d+)(\w+).out', os.path.basename(fn))
  if match:
    jobnumber = int(match.group(1))
    logtype = match.group(2)
    descriptor = ''
    data = []
    print('Reading ' + fn)
    for line in fileinput.input(fn):
      line = line.strip()
      if fileinput.isfirstline():
        descriptor = line.replace(' ', ':')
      else:
        vals = [float(v) for v in line.split()]
        if len(vals) == len(descriptor.split(':')):
          data.append(numpy.array(vals))
        else:
          print('Line is missing entries!')
#    os.remove(fn)
    return logtype, descriptor, data
  else:
    print('Invalid filename ' + fn)

out = ROOT.TFile('out.root', 'RECREATE')
trees = {}
workers = multiprocessing.Pool()
results = workers.imap_unordered(ReadOutFile, sys.argv[1:])

for data in results:
  logtype = data[0]
  if logtype not in trees:
    descriptor = data[1]
    trees[logtype] = ROOT.TNtupleD(logtype, logtype, descriptor)
    print('{0} has {1} columns'.format(logtype, len(descriptor.split(':'))))
  for d in data[2]:
    trees[logtype].Fill(d)

for t in trees:
  print('{0} has {1} entries'.format(t, trees[t].GetEntries()))
out.Write()
