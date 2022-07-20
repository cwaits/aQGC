#---------------------------------------------------------
#run after running binMerger.py to merge bins for SM and other backgrounds,
#then pass in those background root files and INT or QUAD hist file
#and rebin INT or QUAD with the background binning
#---------------------------------------------------------
from sys import argv
from ROOT import *
from array import array
from binMerger import getBinning
#---------------------------------------------------------
if len(argv)>1:
    SMName=(argv[1])
else:
    print "Please give input root file with background histograms as second argument."
    print "Please give input root file with signal histograms as third argument."
    print "Please give the name of the histogram to be merged as fourth argument"
    print "Example: python2.7 binMerger.py path/to/SM_hists path/to/signal_hists R_dijet_mass"
    exit()
if len(argv)>2:
    signalName=str(argv[2])
if len(argv)>3:
    targetName=str(argv[3])
else:
    exit()
#---------------------------------------------------------
#get background binning
f=TFile(SMName)
o=f.Get(targetName)
a=getBinning(o)
print "background hist: ", a
f.Close()

#open signal file
f=TFile(signalName, "UPDATE")
o=f.Get(targetName)
print "Before rebinning: ", getBinning(o)
#rebin
o=o.Rebin(len(a)-1, o.GetName(), a)
print "signal hist: ", getBinning(o)
#print TFile.kOverwrite
f.Write(o.GetName(),TFile.kOverwrite)
f.Close()

