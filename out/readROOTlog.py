# this is an example of how to read PENTrack output with the ROOTlog option enabled
# run with
# python readROOTlog.py *.root

import sys
import ROOT

chain = ROOT.TChain('neutronend') # TChain is used to read a ROOT tree spanning multiple files
for fn in sys.argv[1:]: # read filenames from command line parameters and add them to TChain
  chain.Add(fn)

print('{0} neutrons simulated'.format(chain.GetEntries()))
print('{0} neutrons hit outer boundaries'.format(chain.GetEntries('stopID == -2')))
print('{0} neutrons encountered a serious error!'.format(chain.GetEntries('stopID == -3 || stopID == -6 || stopID == -7')))


# examples of reading config variables:

f = ROOT.TFile(sys.argv[1]) # read config variables from first file

configstr = f.Get('config/GLOBAL/ROOTlog') # check if the ROOTlog option was enabled
if not configstr or configstr.String().Data().split()[0] != '1':
  print('File {0} was created without the ROOTlog option enabled and might be incompatible with this script!'.format(sys.argv[1]))

geometrystr = f.Get('config/GEOMETRY/2') # print material of first solid
if geometrystr:
  print('Material of solid 2 is {0}'.format(geometrystr.String().Data().split()[1]))

sourcestr = f.Get('config/SOURCE/ActiveTime') # print time source was active for
if sourcestr:
  print('Particle source was active for {0} seconds'.format(sourcestr.String().Data().split()[0]))