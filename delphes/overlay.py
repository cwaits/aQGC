# run with "python overlay.py <delphes_output_dir>/* <delphes_output_hists_dir>"
# <delphes_output_dir> is a directory with root files containing the histograms for every sample you want to plot
# <delphes_output_hists_dir> is the directory where you wish to save all the histograms at
# if <delphes_output_hists_dir> is in a parent directory to the directory where this script is run
# then you will need to edit the path variable
from ROOT import *
from sys import argv
import os

norm=False

# enter the cross section, luminosity and number of events for each sample in the same order it is listed in <delphes_output_dir>
cs = [ 1,1,1]
lumi = [1,1,1]
events = [1, 1, 1]
# edit this to point to your delphes dir
path = '/raid01/users/cwaits/MCProd/delphes/'+str(argv[-1])

inputs=argv[1:-1]

f0=TFile(inputs[0])
keys=f0.GetListOfKeys()

for k in keys:
    c=TCanvas()
    m=0
    for f in inputs:
        print f
        i=inputs.index(f)
        f=TFile(f)
        SetOwnership(f,False)
        o=f.Get(k.GetName())

        if type(o)!=type(TH1F()): continue

        try:
            if norm: o.Scale(1./o.Integral(0,o.GetNbinsX()+2))
        except: continue
        o.SetLineColor(i+1)
        o.SetLineStyle(i+1)
	o.Scale((cs[i]*lumi[i])/events[i])
        if i==0:
            first=o
            o.Draw()
        else: o.Draw("SAME")

        m=max(m,o.GetMaximum())
    first.SetMaximum(1.25*m)
    c.Modified()
    filename = os.path.join(path, o.GetName()+'.pdf')
    c.SaveAs(filename)
