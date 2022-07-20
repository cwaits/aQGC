maxEvents=9E9
DEBUG=False
#---------------------------------------------------------
#choose event selection

#two muon, fully hadronic decay channel
#requires exactly 2 muons, a dimuon mass > 106 GeV, and exactly 2 Wjets (50<j_M<100)
hadronicDimuonChannel=True
#hadronicDimuonChannel=False

#semi-leptonic decay channel
#requires one or more muons, exactly one electron, exactly one Wjet
#semileptonicChannel=True
semileptonicChannel=False

#testing channel
#testing=True
testing=False

#---------------------------------------------------------
#normalization
#use to scale histogram for limits settings
#norm=[cross-section, luminosity, # of events]

#backgrounds
SM6=[2.09, 1000000, 100000]
ZZ6=[0.0080156, 1000000, 100000]
Quadboson6=[0.002399273966, 1000000, 100000]
Diboson6=[0.00181387522, 1000000, 100000]

SM10=[2.97534340446, 10000000, 100000]
ZZ10=[0.0032819, 10000000, 100000]
Quadboson10=[0.001880485374, 10000000, 100000]
Diboson10=[0.0019859432681, 10000000, 100000]

#6 TeV signal
INT6_T1=[0.0961303843956, 1000000, 100000]
QUAD6_T1=[1.30920764263, 1000000, 100000]
INT6_T2=[0.118476607978, 1000000, 100000]
QUAD6_T2=[0.76393138231, 1000000, 100000]

#10 TeV signal
INT10_T0=[ 1.04188805092, 10000000, 100000]
QUAD10_T0=[ 186.99326448, 10000000, 100000]
INT10_T1=[ 0.373949736706, 10000000, 100000]
QUAD10_T1=[ 29.7886385348, 10000000, 100000]
INT10_T2=[ 0.446552786477, 10000000, 100000]
QUAD10_T2=[ 16.3462335145, 10000000, 100000]
INT10_T6=[ 0.73951283237, 10000000, 100000]
QUAD10_T6=[ 150.574452024, 10000000, 100000]
INT10_T7=[ 0.891624621914, 10000000, 100000]
QUAD10_T7=[ 99.1337741295, 10000000, 100000]

processes={'SM6':SM6, 'ZZ6':ZZ6, 'Quadboson6':Quadboson6, 'Diboson6':Diboson6, 'SM10':SM10, 'ZZ10':ZZ10, 'Quadboson10':Quadboson10, 'Diboson10':Diboson10, 
'INT6_T1':INT6_T1, 'QUAD6_T1':QUAD6_T1, 'INT6_T2':INT6_T2, 'QUAD6_T2':QUAD6_T2, 
'INT10_T0':INT10_T0, 'QUAD10_T0':QUAD10_T0, 'INT10_T1':INT10_T1, 'QUAD10_T1':QUAD10_T1, 'INT10_T2':INT10_T2, 'QUAD10_T2':QUAD10_T2, 'INT10_T6':INT10_T6, 
'QUAD10_T6':QUAD10_T6, 'INT10_T7':INT10_T7, 'QUAD10_T7':QUAD10_T7,
'None':False, 'False':False}
#---------------------------------------------------------
jetType='VLCjetR10_inclusive'

particle={
    11:"e",
    13:"mu",
    22:"gamma", #using this for all jets
    21:"j",
    23:"z",
    24:"w",
    12:"nu"
}

#OBJECT SELECTION
pTMin={
    11:1,
    13:1,       
    21:5,       #using this for all jets
    22:1
}

etaMax={11:2.5,
        13:2.5,
        21:10,
        22:1
        }
#---------------------------------------------------------
import pdb
import operator
from sys import argv
#---------------------------------------------------------
from ROOT import *
TH1.SetDefaultSumw2()
import math
from array import array
import time
from cut_combo import cut_dict

gSystem.Load("libDelphes")
gStyle.SetOptStat(0)
#---------------------------------------------------------
def printHist(h):
    for i in range(h.GetNbinsX()+2):
        print h.GetBinContent(i),

#---------------------------------------------------------
#return TLorentzVector corresponding to sum of inputs
def parentConstructor(a,b):
    return a.P4()+b.P4()
#---------------------------------------------------------
def getParents(p):
    result=[p]

    motherIndices=[]
    if p.M1!=-1 and event.Particle[p.M1].PID==p.PID:
        motherIndices.append(p.M1)
    if p.M2!=-1 and event.Particle[p.M2].PID==p.PID:
        motherIndices.append(p.M2)
    result+=[getParents(event.Particle[i]) for i in motherIndices]

    return result

def isBeamRemnant(p):
    parents=getParents(p)
    while type(parents)==type([]): parents=parents[-1]
    return parents.Status==4
#---------------------------------------------------------
def sortFunc(x):
    return (x.PT)*math.cosh(x.Eta)
#---------------------------------------------------------
unstable=[6,23,24,25]

def getDaughters(p):
    result=[]
    if event.Particle[p.D1].PID==p.PID or event.Particle[p.D1].PID in unstable:
        result+=getDaughters(event.Particle[p.D1])
    else:
        result+=[event.Particle[p.D1]]

    if p.D2!=p.D1:
        if event.Particle[p.D2].PID==p.PID or event.Particle[p.D2].PID in unstable:
            result+=getDaughters(event.Particle[p.D2])
        else:
            result+=[event.Particle[p.D2]]

    return result
#---------------------------------------------------------
def hadronicDecay(p):
    #print p.PID,p.Status,p.D1,p.D2                                                                                           
    if event.Particle[p.D1].PID==p.PID:
        return hadronicDecay(event.Particle[p.D1])
    elif event.Particle[p.D2].PID==p.PID:
        return hadronicDecay(event.Particle[p.D2])
    else:
        #print tree.Particle[p.D1].PID                                                                                        
        if abs(event.Particle[p.D1].PID)<7:
            #print "q"                                                                                                        
            return True
        else:
            #print "l"                                                                                                        
            return False

#---------------------------------------------------------
def truth_selector(truthElectrons, truthMuons, truthWs, truthZs, truthPhotons, truthNeutrinos, particles):
    for p in particles:
        if p.Status==1 and abs(p.PID)==11:
            truthElectrons.append(p)
        elif p.Status==1 and abs(p.PID)==13:
            truthMuons.append(p)
        elif abs(p.Status)==22 and abs(p.PID)==24:
            truthWs.append(p)
        elif abs(p.Status) in range(21,30) and abs(p.PID)==23:
            truthZs.append(p)
        elif abs(p.Status)==1 and abs(p.PID)==22:
            truthPhotons.append(p)
        elif p.Status==1 and (abs(p.PID)==12 or abs(p.PID)==14 or abs(p.PID)==16):
            truthNeutrinos.append(p)
    return truthElectrons, truthMuons, truthWs, truthZs, truthPhotons, truthNeutrinos
#---------------------------------------------------------
def reco_selector(ROOTArray, pTCut, etaCut):
    #this is a little hacky, come up with smoother function later
    rejects=[]
    particles=[p for p in ROOTArray]
    for i in range(len(particles)):
        if particles[i].PT<pTCut or abs(particles[i].Eta)>etaCut:
            rejects.append(particles[i])
    for i in range(len(rejects)):
        particles.remove(rejects[i])
    return particles
#---------------------------------------------------------

if __name__=='__main__':
    start=time.time()
    inputName=str(argv[1])
    if len(argv)>2:
        sample=processes[str(argv[2])]
    else:
        sample=False
    if len(argv)>3:
        if str(argv[3]) == gridScan:
            gridScan = True
        else:
            gridScan = False
    else:
        gridScan =False
    print "gridScan = ",gridScan
    #if len(argv)>3:
    #    tag=str(argv[3])
    #else:
    tag=''
    outputName=inputName.replace('.root',(bool(tag)*('.%s'%tag))+'.hist.root')

    print "inputName:",inputName
    print "sample:", sample
    #print "tag:",tag

    f=TFile(inputName)
    t=f.Delphes
    output=TFile(outputName,"RECREATE")

    t.GetEntry(0)
    sqrtS=0
    sqrtS+=t.Particle[0].E
    sqrtS+=t.Particle[1].E
    #sets upper limit for bin range
    bin_range = sqrtS/2

    #declares truth-level pT, p, eta, and multiplicity histograms for e,mu,W,Z,gamma
    h_truth={}
    for p in [11,13,21,22,23,24,12]:
        h_truth[p]={}
        h_truth[p]['pT'] =TH1F('T_%s_pT'%particle[p],';Truth %s pT [GeV];Events'%particle[p], 10, 0, bin_range)
        h_truth[p]['p']  =TH1F('T_%s_p'%particle[p],';Truth %s p [GeV];Events'%particle[p], 10, 0, bin_range)
        h_truth[p]['eta']=TH1F('T_%s_eta'%particle[p],';Truth %s #eta;Events'%particle[p], 20, -4, 4)
        h_truth[p]['cos']=TH1F('T_%s_cos'%particle[p],';Truth %s cos(#theta);Events'%particle[p], 10, -1, 1)
        h_truth[p]['mult']=TH1F('T_%s_mult'%particle[p],';Truth %s multiplicity;Events'%particle[p], 7, -0.5, 6.5)
    
    #reco-level histograms for e,mu,jets,and photons
    h_reco={}
    for p in [11,13,21,22]:
        h_reco[p]={}
        h_reco[p]['pT']={}
        h_reco[p]['p']={}
        h_reco[p]['eta']={}
        h_reco[p]['cos']={}
        h_reco[p]['mass']={}
        h_reco[p]['mult']=TH1F('R_%s_mult'%particle[p],';%s multiplicity;Events'%particle[p], 7, -0.5, 6.5)
        for i in ['I']+range(4):
            h_reco[p]['pT'][i]=TH1F('R_%s_Pt_%s'%(particle[p],i),';%s pT [GeV];Events'%particle[p], 200, 0, bin_range)
            h_reco[p]['p'][i]=TH1F('R_%s_P_%s'%(particle[p],i),';%s P [GeV];Events'%particle[p], 20, 0, bin_range)
            h_reco[p]['eta'][i]=TH1F('R_%s_eta_%s'%(particle[p],i),';%s  #eta;Events'%particle[p], 20, -4, 4)
            h_reco[p]['cos'][i]=TH1F('R_%s_cos_%s'%(particle[p],i),';%s cos(#theta);Events'%particle[p], 10, -1, 1)
            h_reco[p]['mass'][i]=TH1F('R_%s_mass_%s'%(particle[p],i),';%s mass [GeV];Events'%particle[p], 20, 0, 200)

    #truth-level missing energy histograms
    T_missingEt = TH1F('T_missingEt', ';Missing Transverse Energy [GeV];Events', 200, 0, bin_range)
    T_missingE = TH1F('T_missingE', ';Missing Energy [GeV];Events', 200, 0, bin_range)

    #reco-level histograms for missing energy
    R_missingET = TH1F('R_MissingET', ';Missing Transverse Energy [GeV];Events', 20, 0, sqrtS)
    R_missingE = TH1F('R_MissingE', ';Missing Energy [GeV];Events', 20, 0, sqrtS)
    R_missingMass = TH1F('R_MissingMass', ';Missing Mass [GeV];Events' , 20, 0, sqrtS)
    R_missingP = TH1F('R_MissingP', ';Missing Momentum [Gev];Events', 20, 0, sqrtS)

    #histograms for OS pairs
    R_ee_pT = TH1F('ee_pT', ';pT [GeV];Events', 200, 0, bin_range)
    R_ee_eta = TH1F('ee_Eta', ';Eta;Events', 20, -4, 4)
    R_ee_mass = TH1F('ee_mass', ';Mass [GeV];Events', 200, 0, 200)
    R_ee_deltaEta = TH1F('ee_Delta Eta', ';Delta Eta;Events', 20, -4, 4)
    R_ee_multiplicity = TH1F('multiplicity', ';Multiplicity;Events', 5, -.5, 4.5)

    R_mumu_pT = TH1F('mumu_pT', ';pT [GeV];Events', 200, 0, bin_range)
    R_mumu_eta = TH1F('mumu_Eta', ';Eta;Events', 20, -4, 4)
    R_mumu_mass = TH1F('mumu_mass', ';Mass [GeV];Events', 20, 0, sqrtS)
    R_mumu_deltaEta = TH1F('mumu_Delta Eta', ';Delta Eta;Events', 20, -4, 4)
    R_mumu_multiplicity = TH1F('mumu_multiplicity', ';Multiplicity;Events', 5, -.5, 4.5)
    R_mumu_p = TH1F('mumu_P', ';P [GeV];Events', 20, 0, sqrtS)
    R_mumu_cos = TH1F('mumu_cos', ';cos(theta);Events', 10, -1, 1)

    R_emu_pT = TH1F('emu_pT', ';pT [GeV];Events', 200, 0, bin_range)
    R_emu_eta = TH1F('emu_Eta', ';Eta;Events', 20, -4, 4)
    R_emu_mass = TH1F('emu_mass', ';Mass [GeV];Events', 200, 0, 200)
    R_emu_deltaEta = TH1F('emu_Delta Eta', ';Delta Eta;Events', 20, -4, 4)
    R_emu_multiplicity = TH1F('emu_multplicity', ';Multiplicity;Events', 5, -.5, 4.5)

    #histograms for beam remnants
    T_beamRemnants_pT = TH1F('T_beamRemnants_Pt', ';pT [GeV];Events', 20, 0, bin_range)
    T_beamRemnants_p = TH1F('T_beamRemnants_p', ';p [GeV];Events', 200, 0, bin_range)
    T_beamRemnants_eta = TH1F('T_beamRemnants_eta', ';Eta;Events', 40, -10, 10)
    T_beamRemnants_multiplicity = TH1F('T_beamRemnants_multiplicity', ';Multiplicity;Events', 7, -.5, 6.5)
    T_beamRemnants_cosh = TH1F('T_beamRemnants_cosh', ';Cosh(Eta)/pT;Events', 20, 0, math.cosh(10)/(10*bin_range))
    T_beamRemnants_mass = TH1F('T_beamRemnants_mass', ';mass [GeV];Events', 300, 0, bin_range)

    #histograms for non-beam remnants
    T_nonBeamRemnants_pT = TH1F('T_nonBeamRemnants_Pt', ';pT [GeV];Events', 20, 0, bin_range)
    T_nonBeamRemnants_p = TH1F('T_nonBeamRemnants_p', ';p [GeV];Events', 200, 0, bin_range)
    T_nonBeamRemnants_eta = TH1F('T_nonBeamRemnants_eta', ';Eta;Events', 40, -10, 10)
    T_nonBeamRemnants_multiplicity = TH1F('T_nonBeamRemnants_multiplicity', ';Multiplicity;Events', 7, -.5, 6.5)
    T_nonBeamRemnants_cosh = TH1F('T_nonBeamRemnants_cosh', ';Cosh(Eta)/pT;Events', 20, 0, math.cosh(10)/(10*bin_range))

    #histograms for reco beam remnants
    R_beamRemnants_pT = TH1F('R_beamRemnants_Pt', ';pT [GeV];Events', 20, 100, bin_range)
    R_beamRemnants_p = TH1F('R_beamRemnants_p', ';p [GeV];Events', 200, 100, bin_range)
    R_beamRemnants_eta = TH1F('R_beamRemnants_eta', ';Eta;Events', 40, -2.5, 2.5)
    R_beamRemnants_multiplicity = TH1F('R_beamRemnants_multiplicity', ';Multiplicity;Events', 7, -.5, 6.5)
    R_beamRemnants_cosh = TH1F('R_beamRemnants_cosh', ';Cosh(Eta)/pT;Events', 20, 0, math.cosh(10)/(10*bin_range))
    R_beamRemnants_mass = TH1F('R_beamRemnants_mass', ';mass [GeV];Events', 600, 0, sqrtS)

    #histograms for reco non-beam remnants
    R_nonBeamRemnants_pT = TH1F('R_nonBeamRemnants_Pt', ';pT [GeV];Events', 20, 0, 100)
    R_nonBeamRemnants_p = TH1F('R_nonBeamRemnants_p', ';p [GeV];Events', 200, 0, 100)
    R_nonBeamRemnants_eta = TH1F('R_nonBeamRemnants_eta', ';Eta;Events', 40, -2.5, 2.5)
    R_nonBeamRemnants_multiplicity = TH1F('R_nonBeamRemnants_multiplicity', ';Multiplicity;Events', 7, -.5, 6.5)
    R_nonBeamRemnants_cosh = TH1F('R_nonBeamRemnants_cosh', ';Cosh(Eta)/pT;Events', 20, 0, math.cosh(10)/(10*bin_range))

    #histograms for Wjets
    R_Wjets_pT = TH1F('R_Wjets_Pt', ';pT [GeV];Events', 20, 0, bin_range)
    R_Wjets_p = TH1F('R_Wjets_p', ';p [GeV];Events', 20, 0, bin_range)
    R_Wjets_eta = TH1F('R_Wjets_eta', ';Eta;Events', 20, -2.5, 2.5)
    R_Wjets_multiplicity = TH1F('R_Wjets_multiplicity', ';Multiplicity;Events', 7, -.5, 6.5)
    R_Wjets_mass = TH1F('R_Wjets_mass', ';mass [GeV];Events', 50, 0, 250)
    R_Wjets_mass_pairs = TH1F('R_Wjets_mass_pairs', ';mass [GeV];Events', 20, 0, sqrtS)
    R_Wjets_p_leading = TH1F('R_Wjets_p_leading', ';p [GeV];Events', 20, 0, bin_range)
    R_Wjets_p_sub = TH1F('R_Wjets_p_sub', ';p [GeV];Events', 20, 0, bin_range)
    R_Wjets_cos_leading = TH1F('R_Wjets_cos_leading', ';cos(theta);Events', 10, -1, 1)
    R_Wjets_cos_sub = TH1F('R_Wjets_cos_sub', ';cos(theta);Events', 10, -1, 1)
    R_WJets_p_pair = TH1F('R_Wjets_p_pairs', ';p [GeV];Events', 20, 0, sqrtS)
    R_Wjets_cos_pair = TH1F('R_Wjets_cos_pairs', ';cos(theta);Events',10, -1, 1)

    #histograms for reconstructed W's (definition depedent on event selection being used)
    R_W_pT = TH1F('R_W_pT', ';pT [GeV]; Events', 20, 0, bin_range)
    R_W_p = TH1F('R_W_p', ';p [GeV]; Events', 20, 0, bin_range)
    R_W_mass = TH1F('R_W_mass', ';mass [GeV];Events', 50, 0, 250)
    R_W_eta = TH1F('R_W_eta', ';Eta; Events', 20, -2.5, 2.5)
    R_W_multiplicity = TH1F('R_W_mulitplicity', ';Multiplicity;Events', 7, -.5, 6.5)

    #dijet histograms for fitting
    #a=array('f', [0]+list(range(1000, 3001, 500))+[6000])
    R_dijet_mass = TH1F('R_dijet_mass','mass [GeV];Events', 10, 0, sqrtS)
    R_dijet_mass_10GeVbinning = TH1F('R_dijet_mass_10GeVbinning', ';mass [GeV];Events', 600, 0, sqrtS)

    CutFlow = TH1F('CutFlow', ';Cut;Events', 10, .5, 10.5)
   
    #misc histograms
    T_WW_mass = TH1F('WW_mass', ';Mass [GeV];Events', 20, 0, sqrtS) 
    Test_hist = TH1F('Test_hist', ';Mass;Events', 1, 0, 1) 
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    m=0
    countA=0
    countB=0
    for event in f.Delphes:
        m+=1
        #truth level particles
        truthElectrons=[]
        truthMuons=[]
        truthWs=[]
        truthZs=[]
        truthPhotons=[]
        truthNeutrinos=[]
        truth_selector(truthElectrons, truthMuons, truthWs, truthZs, truthPhotons, truthNeutrinos, event.Particle)

        #special truth level particles
        beamRemnantMuons=[p for p in truthMuons if isBeamRemnant(p)]
        nonBeamRemnantMuons=[p for p in truthMuons if not isBeamRemnant(p)]

        #reco level particles
        electrons=reco_selector(event.Electron, pTMin[11], etaMax[11]) 
        muons=reco_selector(event.Muon, pTMin[13], etaMax[13])
        photons=reco_selector(event.Photon, pTMin[22], etaMax[22])
        jets=reco_selector(event.__getattr__(jetType), pTMin[21], etaMax[21])
        leptons=electrons+muons
        
        # sort object lists by momentum
        for collection in [truthElectrons, truthMuons, truthWs, truthZs, truthPhotons, truthNeutrinos, beamRemnantMuons, nonBeamRemnantMuons]: collection.sort(key=sortFunc, reverse=True)
        for collection in [electrons, muons, photons, jets]: collection.sort(key=sortFunc, reverse=True)
        
        #special reco level particles
        W_jets=[p for p in jets if 50<p.P4().M()<100] 

        #--------------------------------------
        #hadronic dimuon channel event selection
        if hadronicDimuonChannel:
            #cuts applied through a dictionary made from gridScan.py and imported from cut_combo.py
            if gridScan:
                 muons = [p for p in muons if p.P4().P()<=cut_dict['mu_P']]
                 muons = [p for p in muons if abs(p.Eta)<=cut_dict['mu_eta']]
                 W_jets = [p for p in W_jets if p.P4().P()>=cut_dict['jet_P']]
                 W_jets = [p for p in W_jets if abs(p.Eta)<=cut_dict['jet_eta']]
            #fill CutFlow 
            CutFlow.Fill(1)
            if len(muons)==2:
                muon_4vec=muons[0].P4()+muons[1].P4()
                CutFlow.Fill(2)
                if muon_4vec.M()>=106:
                    CutFlow.Fill(3)
                    if len(W_jets)==2:
                        CutFlow.Fill(4)
            #select events
            if len(muons)!=2:
                continue
            if muon_4vec.M()<106:
                continue
            if len(W_jets)!=2:
                continue

            #fill channel-specific histograms
            R_mumu_mass.Fill(muon_4vec.M())
            R_mumu_p.Fill(muon_4vec.P())
            R_mumu_pT.Fill(muon_4vec.Pt())
            R_mumu_eta.Fill(muon_4vec.Eta())
            
            dijet=W_jets[0].P4()+W_jets[1].P4()
            R_dijet_mass.Fill(dijet.M())
            R_dijet_mass_10GeVbinning.Fill(dijet.M())
        #--------------------------------------
        #semi-leptonic channel event selection
        if semileptonicChannel:
            #find reco-level missing momentum
            #final and initial-state 4-vectors
            P4_i=TLorentzVector()
            P4_i+=event.Particle[0].P4()
            P4_i+=event.Particle[1].P4()
            P4_f=TLorentzVector()
            for i in range(len(leptons)):
                P4_f+=leptons[i].P4()
            for i in range(len(photons)):
                P4_f+=photons[i].P4()
            for i in range(len(jets)):
                P4_f+=jets[i].P4()
            missingP=(P4_f-P4_i).P()
            missingP_cut=sqrtS
            #cuts applied through a dictionary made from gridScan.py and imported from cut_combo.py
            if gridScan:
                muons = [p for p in muons if abs(p.Eta)<=cut_dict['mu_eta']]
                muons = [p for p in muons if p.P4().P()<=cut_dict['mu_P']]
                electrons = [p for p in electrons if p.P4().P()>=cut_dict['e_P']]
                missingP_cut = cut_dict['missing_P']
            #fill CutFlow
            CutFlow.Fill(1)
            if len(muons)>=1:
                CutFlow.Fill(2)
                if len(electrons)==1:
                    CutFlow.Fill(3)
                    if len(W_jets)==1:
                        CutFlow.Fill(4)
                        if missingP < missingP_cut:
                            CutFlow.Fill(5)

            #select events
            if len(muons)<1:
                continue
            if len(electrons)!=1:
                continue
            if len(W_jets)!=1:
                continue
            if missingP > missingP_cut:
                continue

            #fill channel-specific histograms
            #calling the electron and W-jet 4 vector the "dijet"
            dijet = electrons[0].P4() + W_jets[0].P4()
            R_dijet_mass.Fill(dijet.M())
        #--------------------------------------
        #testing channel
        if testing:
            ptCut=800
            truthWs=[p for p in truthWs if abs(p.Eta)<2.5]
           
            #fill CutFlow
            CutFlow.Fill(1)
            if len(truthWs)==2:
                CutFlow.Fill(2)
                if truthWs[0].PT>ptCut and truthWs[1].PT>ptCut:
                    CutFlow.Fill(3)
            
            #select events
            if len(truthWs)!=2:
                continue
            if truthWs[0].PT<ptCut or truthWs[1].PT<ptCut:
                continue

            #fill channel-specific histograms
            T_WW_mass.Fill((truthWs[0].P4()+truthWs[1].P4()).M())
        #--------------------------------------
        #fill all other truth histograms
        h_truth[11]['mult'].Fill(len(truthElectrons))
        h_truth[13]['mult'].Fill(len(truthMuons))
        h_truth[24]['mult'].Fill(len(truthWs))
        h_truth[23]['mult'].Fill(len(truthZs))
        h_truth[22]['mult'].Fill(len(truthPhotons))
        h_truth[12]['mult'].Fill(len(truthNeutrinos))
        T_beamRemnants_multiplicity.Fill(len(beamRemnantMuons))
        T_nonBeamRemnants_multiplicity.Fill(len(nonBeamRemnantMuons))
 
        for i in range(len(truthElectrons)):
            h_truth[11]['pT'].Fill(truthElectrons[i].PT)
            h_truth[11]['p'].Fill(truthElectrons[i].P4().P())
            h_truth[11]['eta'].Fill(truthElectrons[i].Eta)
            h_truth[11]['cos'].Fill(truthElectrons[i].P4().CosTheta()) 

        for i in range(len(truthMuons)):
            h_truth[13]['pT'].Fill(truthMuons[i].PT)
            h_truth[13]['p'].Fill(truthMuons[i].P4().P())
            h_truth[13]['eta'].Fill(truthMuons[i].Eta)
            h_truth[13]['cos'].Fill(truthMuons[i].P4().CosTheta())

        for i in range(len(truthWs)):
            h_truth[24]['pT'].Fill(truthWs[i].PT)
            h_truth[24]['p'].Fill(truthWs[i].P4().P())
            h_truth[24]['eta'].Fill(truthWs[i].Eta)
            h_truth[24]['cos'].Fill(truthWs[i].P4().CosTheta())

        for i in range(len(truthZs)):
            h_truth[23]['pT'].Fill(truthZs[i].PT)
            h_truth[23]['p'].Fill(truthZs[i].P4().P())
            h_truth[23]['eta'].Fill(truthZs[i].Eta)
            h_truth[23]['cos'].Fill(truthZs[i].P4().CosTheta())

        for i in range(len(truthPhotons)):
            h_truth[22]['pT'].Fill(truthPhotons[i].PT)
            h_truth[22]['p'].Fill(truthPhotons[i].P4().P())
            h_truth[22]['eta'].Fill(truthPhotons[i].Eta)
            h_truth[22]['cos'].Fill(truthPhotons[i].P4().CosTheta())

        for i in range(len(truthNeutrinos)):
            h_truth[12]['pT'].Fill(truthNeutrinos[i].PT)
            h_truth[12]['p'].Fill(truthNeutrinos[i].P4().P())
            h_truth[12]['eta'].Fill(truthNeutrinos[i].Eta)
            h_truth[12]['cos'].Fill(truthNeutrinos[i].P4().CosTheta())

        for i in range(len(beamRemnantMuons)):
            T_beamRemnants_pT.Fill(beamRemnantMuons[i].PT)
            T_beamRemnants_p.Fill(beamRemnantMuons[i].P4().P())
            T_beamRemnants_eta.Fill(beamRemnantMuons[i].Eta)
            T_beamRemnants_cosh.Fill(math.cosh(beamRemnantMuons[i].Eta)/beamRemnantMuons[i].PT)

        for i in range(len(nonBeamRemnantMuons)):
            T_nonBeamRemnants_pT.Fill(nonBeamRemnantMuons[i].PT)
            T_nonBeamRemnants_p.Fill(nonBeamRemnantMuons[i].P4().P())
            T_nonBeamRemnants_eta.Fill(nonBeamRemnantMuons[i].Eta)
            T_nonBeamRemnants_cosh.Fill(math.cosh(nonBeamRemnantMuons[i].Eta)/nonBeamRemnantMuons[i].PT)

        truthMissingP=TLorentzVector()
        for nu in truthNeutrinos:
            truthMissingP+=nu.P4()

        T_missingEt.Fill(truthMissingP.Et())
        T_missingE.Fill(truthMissingP.E())
        #--------------------------------------
        #fill all other reco histograms
        R_missingET.Fill(event.MissingET[0].MET)
        R_missingE.Fill(event.MissingET[0].P4().P())

        #final and initial-state 4-vectors
        P4_i=TLorentzVector()
        P4_i+=event.Particle[0].P4()
        P4_i+=event.Particle[1].P4()

        P4_f = TLorentzVector()

        h_reco[11]['mult'].Fill(len(electrons))
        for i in range(len(electrons)):
            P4_f+=electrons[i].P4()
            h_reco[11]['pT']['I'].Fill(electrons[i].PT)
            h_reco[11]['p']['I'].Fill(electrons[i].P4().P())
            h_reco[11]['eta']['I'].Fill(electrons[i].Eta)
            h_reco[11]['cos']['I'].Fill(electrons[i].P4().CosTheta())
            if i<4:
                h_reco[11]['pT'][i].Fill(electrons[i].PT)
                h_reco[11]['p'][i].Fill(electrons[i].P4().P())
                h_reco[11]['eta'][i].Fill(electrons[i].Eta)
                h_reco[11]['cos'][i].Fill(electrons[i].P4().CosTheta())

        h_reco[13]['mult'].Fill(len(muons))
        for i in range(len(muons)):
            P4_f+=muons[i].P4()
            h_reco[13]['pT']['I'].Fill(muons[i].PT)
            h_reco[13]['p']['I'].Fill(muons[i].P4().P())
            h_reco[13]['eta']['I'].Fill(muons[i].Eta)
            h_reco[13]['mass']['I'].Fill(muons[i].P4().M())
            h_reco[13]['cos']['I'].Fill(muons[i].P4().CosTheta())
            if i<4:
                h_reco[13]['pT'][i].Fill(muons[i].PT)
                h_reco[13]['p'][i].Fill(muons[i].P4().P())
                h_reco[13]['eta'][i].Fill(muons[i].Eta)
                h_reco[13]['mass'][i].Fill(muons[i].P4().M())
                h_reco[13]['cos'][i].Fill(muons[i].P4().CosTheta())

        h_reco[22]['mult'].Fill(len(photons))
        for i in range(len(photons)):
            P4_f+=photons[i].P4()
            h_reco[22]['pT']['I'].Fill(photons[i].PT)
            h_reco[22]['p']['I'].Fill(photons[i].P4().P())
            h_reco[22]['eta']['I'].Fill(photons[i].Eta)
            h_reco[22]['cos']['I'].Fill(photons[i].P4().CosTheta())
            if i<4:
                h_reco[22]['pT'][i].Fill(photons[i].PT)
                h_reco[22]['p'][i].Fill(photons[i].P4().P())
                h_reco[22]['eta'][i].Fill(photons[i].Eta)
                h_reco[22]['cos'][i].Fill(photons[i].P4().CosTheta())

        h_reco[21]['mult'].Fill(len(jets))
        for i in range(len(jets)):
            P4_f+=jets[i].P4()
            h_reco[21]['pT']['I'].Fill(jets[i].PT)
            h_reco[21]['p']['I'].Fill(jets[i].P4().P())
            h_reco[21]['eta']['I'].Fill(jets[i].Eta)
            h_reco[21]['mass']['I'].Fill(jets[i].P4().M())
            h_reco[21]['cos']['I'].Fill(jets[i].P4().CosTheta())
            if i<4:
                h_reco[21]['pT'][i].Fill(jets[i].PT)
                h_reco[21]['p'][i].Fill(jets[i].P4().P())
                h_reco[21]['eta'][i].Fill(jets[i].Eta)
                h_reco[21]['mass'][i].Fill(jets[i].P4().M())
                h_reco[21]['cos'][i].Fill(jets[i].P4().CosTheta())

        R_Wjets_multiplicity.Fill(len(W_jets))
        for i in range(len(W_jets)):
            R_Wjets_pT.Fill(W_jets[i].PT)
            R_Wjets_p.Fill(W_jets[i].P4().P())
            R_Wjets_eta.Fill(W_jets[i].Eta)
            R_Wjets_mass.Fill(W_jets[i].P4().M())

        R_missingMass.Fill((P4_f-P4_i).M())
        R_missingP.Fill((P4_f-P4_i).P())
        # end of loop
        #-------------------------------------------------------------------------------------------------
    #normalize target hist by cross-section*lumi/#ofEvents
    norm=1
    print "unscaled: " ,R_dijet_mass.Integral()
    if sample != False:
        norm=sample[0]*sample[1]/sample[2]
    #only need to scale one hist per channel, but will scale them all here for simplicity
    R_dijet_mass.Scale(norm)
    T_WW_mass.Scale(norm)
    print "scaled: ", R_dijet_mass.Integral()
    print 'Number of events:', CutFlow.GetBinContent(1)
    print 'Number of events passing cut 1:', CutFlow.GetBinContent(2)
    print 'Number of events passing cut 2:', CutFlow.GetBinContent(3)
    print 'Number of events passing cut 3:', CutFlow.GetBinContent(4)
    print 'Number of events passing cut 4:', CutFlow.GetBinContent(5)
    print 'Number of events passing cut 5:', CutFlow.GetBinContent(6)
    Test_hist.SetBinContent(1, 10)
    output.Write()
    print("Execution time:",time.time()-start)

