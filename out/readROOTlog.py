# This script is an example of how to read PENTrack output with the ROOTlog option enabled.
# It will print some information about the simulated particles (total number, number with errors),
# plot some particle distributions, and print some configuration parameters
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

c = ROOT.TCanvas('c', 'c') # create canvas to draw on
chain.Draw('Estart') # plot histogram of starting energies
c.Print('spectrum.pdf') # save canvas to file

chain.Draw('tend', 'Estart > 100e-9') # plot histogram of stopping times of all neutrons with energy above 100neV (will overwrite contents of canvas)
chain.Draw('tend', 'Estart < 100e-9', 'same') # plot histogram of stopping times of all neutrons with energy below 100neV on same canvas
c.Print('tend.pdf') # save canvas to file

chain.Draw('zend:xend', 'stopID > 0') # plot spatial distribution of neutrons absorbed in a material
c.Print('absorbed.pdf')

chain.Draw('zend:xend:solidend', '', 'COLZ') # plot spatial distribution of all stopped neutrons with markers colored by part ID
c.Print('stopped.pdf')

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


