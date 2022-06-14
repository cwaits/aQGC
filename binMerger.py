#---------------------------------------------------------
# This script is meant primarily to be run through setLimits.sh,
# so that the SM sample is sure to be run first with low bins settings set to True
# and then the INT and QUAD samples are run after, using array from SM_binning.pkl
# and low bins setting set to False
#---------------------------------------------------------
from sys import argv
from ROOT import *
from array import array
#---------------------------------------------------------
if len(argv)>1:
    inputName=(argv[1])
else:
    print "Please give input root file with histograms as second argument."
    print "Please give the name of the histogram to be merged as third argument"
    print "Example: python2.7 binMerger.py path/to/dir/ R_dijet_mass"
    exit()
if len(argv)>2:
    targetName=str(argv[2])
else:
    exit()
#---------------------------------------------------------
#returns a list of all bins with bin content below threshold
def findLowBins(hist, bin_threshold):
    #exclude underflow/overflow bins
    n_bins=hist.GetSize()-2
    low_bins=[]
    for i in range(1, n_bins+1):
        if hist.GetBinContent(i)<bin_threshold:
             low_bins.append(i)
    return low_bins
#---------------------------------------------------------
#merges bins with their smallest neighbor bin
def mergeBins(hist, low_bins_list, binning):
    n_bins=hist.GetSize()-2
    #binning[i-1] is the lower edge of bin i
    #to avoid changing index issues, will make a list of 
    #bin edges to remove and then remove those values from binning array
    trash_bins=[]
    for b in low_bins_list:
        #merge 1st bin to right neighbor
        if b == 1:
           trash_bins.append(binning[1])
        #merge last bin to left neighbor
        elif b == n_bins: 
            if binning[-2] not in trash_bins:
                trash_bins.append(binning[-2])
        else:
            if hist.GetBinContent(b-1) >= hist.GetBinContent(b+1):
                #don't duplicate bins in trash_bins
                if binning[b] not in trash_bins:
                    #right neighbor is the smaller, remove bin edge between bin b and b+1
                    trash_bins.append(binning[b])
            else:
                #don't duplicate bins in trash_bins
                if binning[b-2] not in trash_bins:
                    trash_bins.append(binning[b-2])
    for b in trash_bins:
        binning.remove(b)
    return binning
#---------------------------------------------------------
#returns array containing the binning of a histogram
def getBinning(hist):
    n_bins=hist.GetSize()-2
    binning=[]
    for i in range(1, n_bins+1):
        binning.append(hist.GetXaxis().GetBinLowEdge(i))
    #include top edge of last bin
    binning.append(hist.GetXaxis().GetBinLowEdge(n_bins) + hist.GetXaxis().GetBinWidth(n_bins))
    a=array('d', binning)
    return a
#---------------------------------------------------------
if __name__=='__main__':
    #get SM file
    f=TFile(inputName, "UPDATE") 
    o=f.Get(targetName)

    #bin threshold, all bins with bin content less than bin_threshold will be merged
    bin_threshold=4

    low_bins_list=findLowBins(o, bin_threshold)
    #repeat until there are no low bins
    while len(low_bins_list) != 0:
        a=getBinning(o)
        #merge the bins in low_bins_list and get the new binning
        new_binning=mergeBins(o, low_bins_list, a)
        print "Empty bins: ", low_bins_list
        print "New binning: ", new_binning
        print "Length: ", len(new_binning)
        print
        #re-bin histogram with new_binning
        o=o.Rebin(len(new_binning)-1, o.GetName(), new_binning)
        low_bins_list=findLowBins(o, bin_threshold)

    f.Write("",TFile.kOverwrite)
    f.Close()
