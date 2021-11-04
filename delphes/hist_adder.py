from ROOT import *
from sys import argv
import os
import pdb
inputs = argv[1:]
f0 = TFile(inputs[0])
keys = f0.GetListOfKeys()

#cross-sections and No. of events for normalization
cs = [1, 1, 1]
events = [1, 1, 1]

#creates list to hold added histograms and output file
Hlist = []
output = TFile("Name", "RECREATE")

for k in keys:
    for f in inputs:
        i = inputs.index(f)
        f = TFile(f)
        #SetOwnership(f,False)
	output.cd()
        o = f.Get(k.GetName()).Clone("h%i"%i)
	o.Scale(cs[i]/events[i])
        SetOwnership(o,False)
        #pdb.set_trace()
    h0 = output.Get('h0')
    h1 = output.Get("h1")
    h2 = output.Get("h2")
    h1 += h0
    h0.Draw()
    h1.Draw("Same")
    h2.Draw("Same")
    
output.Write()
