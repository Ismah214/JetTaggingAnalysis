import sys
import os
import glob
import ROOT
import plotter
import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-s", "--sample", action="store",         dest="sample",    default="ZprimeToTTM" )
parser.add_option("-t", "--test",   action="store_true",    dest="test",      default= False )
parser.add_option("-b", "--batch",  action="store_true",    dest="batch",     default= False )
parser.add_option("-p", "--part",   action="store",         dest="part",      default= "0" )
parser.add_option("-v", "--verb",   action="store_true",    dest="verbose",   default= False)

(options, args) = parser.parse_args()
BATCH = options.batch
TEST = options.test
SAMPLE = options.sample
PART = options.part
VERB = options.verbose

from collections import OrderedDict 
from array import array

from ZprimeToTT_read_helper import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Event
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import InputTree

ROOT.gROOT.SetBatch()
timer = ROOT.TStopwatch()
EOSUSER     = "root://eosuser.cern.ch/"
treeName    = "TreeFatJet"
outDir      = "/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/"
bins = 5000 
if not BATCH:
    if "Zprime" in SAMPLE:
        arg1 = ["500"] if TEST else ["500","750","1000","1250","1500","2000","2500","3000","3500","4000","5000","6000","7000","8000"]
    elif "QCD" in SAMPLE:
        arg1 = ["_part1"] if TEST else ["_part"+str(x) for x in range(0,81)]
        
else:
    if SAMPLE == "ZprimeToTTM":
        if PART == "1":
            arg1 = ["500","750","1000","1250"]
        elif PART == "2":
            arg1 = ["1500","2000","2500","3000","3500"]
        elif PART == "3":
            arg1 = ["4000","5000","6000","7000","8000"]
    

    if SAMPLE == "QCDPt15To7000":
        if PART == "1":
            arg1 = ["_part"+str(x) for x in range(0,21)]
        elif PART == "2":
            arg1 = ["_part"+str(x) for x in range(21,41)]
        elif PART == "3":
            arg1 = ["_part"+str(x) for x in range(41,61)]
        elif PART == "4":
            arg1 = ["_part"+str(x) for x in range(61,81)]

if SAMPLE == "ZprimeToTTM":
    h_pt_deepAK8    = ROOT.TH2F("h_pt_deepAK8", "Top-tagging DeepAK8 vs p_{T}; p_{T}; DeepAK8", 50, 0., 5000., 2000, 0., 1.)
    h_pt_tau32      = ROOT.TH2F("h_pt_tau32", "Top-tagging #tau_{32} vs p_{T}; p_{T}; #tau_{32}", 50, 0., 5000., 2000, 0., 1.)
    histoList = [h_pt_deepAK8,h_pt_tau32] 
    
if SAMPLE == "QCDPt15To7000":
    h_pt_T_deepAK8    = ROOT.TH2F("h_pt_T_deepAK8", "Top-tagging DeepAK8 vs p_{T}; p_{T}; DeepAK8", 50, 0., 5000., 2000, 0., 1.)
    h_pt_T_tau32      = ROOT.TH2F("h_pt_T_tau32", "Top-tagging #tau_{32} vs p_{T}; p_{T}; #tau_{32}", 50, 0., 5000., 2000, 0., 1.)
    histoList = [h_pt_T_deepAK8,h_pt_T_tau32]

p = [arg1]

### make sure no same part in the list ###
def checkList(gp, gpfinal):
    for p in gpfinal:
        if getFinal(gp)._index==p._index:
            return False          
    return True

### check if gp is within AK8 jet ###
def inAK8(gp):
    for fj in ak8jet:
        if gp.DeltaR(fj)<0.6:
            return True
    return False

### get final decay of part using dauIdx ###
def getFinal(gp):
  #seeing if each particle has a daughter and indexing (idx) the value of the daughter id
    for idx in gp.dauIdx:
    #placing the daughters in array if exist 
        dau = particles[idx]
        #if a daughter exists, then go to the top of the stack 
        if dau.pdgId == gp.pdgId:
            #loop again until daughter with no daughter is obtained, i.e. final decay / simulation ends   
            return getFinal(dau) 
    #this line takes effect when the final decay part is found in the loop and allows the loop to exit
    return gp 

### check if gp -> W -> qq ###
def decayToW(gp):
    if len(gp.dauIdx) == 0:
        raise ValueError('Particle has no daughters!')
    for idx in gp.dauIdx:
        dp=particles[idx]
        dpart=getFinal(dp)
        if abs(dpart.pdgId) == 24 and decayTo2Q(dpart):
            return True
    return False

### check if gp -> b ###
def decayToB(gp):
    if len(gp.dauIdx) == 0:
        raise ValueError('Particle has no daughters!')
    for idx in gp.dauIdx:
        dpart=particles[idx]
        if abs(dpart.pdgId) == 5:
            return True
    return False

## check if gp decays to at least 2 quarks ###
def decayTo2Q(gp):
    if len(gp.dauIdx) == 0:
        raise ValueError('Particle has no daughters!')
    for idx in gp.dauIdx:
        dpart=particles[idx]
        if abs(dpart.pdgId) < 6:            
            gp.dauIdx.reverse()
            for idx in gp.dauIdx:
                dpart=particles[idx]
                if abs(dpart.pdgId)<6:
                    return True
    return False

### LOOP OVER MASS ###
timer.Start()
for y in p[0]:
    print ""
    print "****************************************"
    print ""
    print "START Processing sample:"
    print "%s%s"%(SAMPLE,y)
    print ""
    print "****************************************"
    INDIR   =  "/eos/user/n/nsubrama/AnaJetTagging/Ntuples/Ntuple_"
    INDIR   += SAMPLE
    INDIR   += y
    print INDIR
    ### GET LIST OF FILES ###
    inFileList = [EOSUSER+f for f in glob.glob(INDIR+".root")]
    
    tree = ROOT.TChain(treeName)
    
    for inFile in inFileList:
        tree.Add(inFile)
    
    if "Zprime" in SAMPLE:    
        INDIR2  =  INDIR
        INDIR2  += "_ext1"
        inFileListExt = [EOSUSER+f for f in glob.glob(INDIR2+".root")]
        
        for inFile in inFileListExt:
          tree.Add(inFile)  
        
    ### USE TCHAIN AND SETUP TTREEREADER ###
    inTree  = InputTree(tree)
    nEvents = inTree.GetEntries()
    if TEST:  nEvents = 1000 
    ### LOOP OVER EVENTS ###
    print "Total events: ", nEvents
    #xrange returns xrange objects so list operations cannot be performed on them. range returns list but is slower than xrange.
    for i in xrange(0,nEvents):     
        if i%10000==0: 
            print "Running %d out of %d" %(i, nEvents)
        
        evt       = Event(inTree,i)
        if "Zprime" in SAMPLE:
            particles = Collection(evt, "GenPart")
            fatjet    = Collection(evt, "FatJet")
            jetak8    = Collection(evt, "GenJetAK8")
            
            ###############################################
            ### efficiency with respect to top quark pt ### 
            ###############################################
            
            AssignDauIdx(particles)
            
            ak8jet  =   [x for x in fatjet if (x.pt>200. and abs(x.eta)<2.4)]
            nT      =   [x for x in particles if (abs(x.eta)<2.4 and abs(x.pdgId)==6)]
            
            if len(nT)<1: 
                continue 
                #continue to next loop if no top quarks in this set
            
            Tfinal = []
            
            for t in nT:
                if checkList(t,Tfinal) and decayToW(t) and decayToB(t): 
                #if t is a new entry into Tfinal, and it decays to a W boson and a b quark
                    Tfinal.append(getFinal(t)) 
                    #then add this final decaying top into Tfinal
            
            if len(Tfinal)<1:
            #continue to next loop if no candidate top quarks in this set 
                continue 
            
            for t in Tfinal:
                if inAK8(t):
                    for fj in ak8jet:
                    #fill the histos with the data #discriminator value of deepAK8 of top from QCD vs pT of the top 
                    #n-subjettiness value is the discriminator value of the top from QCD; vs pT of the top
                        h_pt_deepAK8.Fill(fj.pt, fj.deepTag_TvsQCD)                                    
                        h_pt_tau32.Fill(fj.pt, fj.tau32)

        if "QCD" in SAMPLE: 
            
            fatjet  = Collection(evt, "FatJet")
            ak8jet  = [x for x in fatjet if (x.pt>200. and abs(x.eta)<2.4)]
            
            for fj in ak8jet:
                    #fill the histos with the data #discriminator value of deepAK8 of top from QCD vs pT of the top 
                    h_pt_T_deepAK8.Fill(fj.pt, fj.deepTag_TvsQCD)
                    #n-subjettiness value is the discriminator value of the top from QCD; vs pT of the top            
                    h_pt_T_tau32.Fill(fj.pt, fj.tau32)    
    
color = [46,42,30,36,38]
lgd = [ROOT.kFullSquare, ROOT.kFullCircle, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kOpenCircle, ROOT.kOpenSquare]

#################
### ROOT FILE ###
#################
filename = SAMPLE + "_"+ PART +"_tagger_pt"
outFile = ROOT.TFile("%s%s.root" %(outDir, filename), "RECREATE") #/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/ZprimeToTTM_T_tagger_pt.root
for a in histoList:
    a.SetStats(0)
    a.SetOption("COLZ")
    a.Write()
outFile.Close()

#Start plotting the histos
## CANVAS ##

#legname = ["#bf{#DeepAK8}", "#bf{max#DeltaR(AK8,q)<0.8}"]
#leg = BeginLeg(len(legname))
i=0
c = ROOT.TCanvas("","",877, 620)
for a in histoList:
    outputName = SAMPLE + "_"+ PART
    if "deep" in a.GetName():
        outputName += "_deepAK8"
    else:
        outputName += "_tau32"
    a.SetStats(0)
    #a.SetStatisticOption(0)
    a.SetLineWidth(1)
    a.SetLineColor(color[i])
    a.SetLineColor(color[i])
    a.SetMarkerColor(color[i])
    a.SetMarkerStyle(lgd[i])
    a.SetFillColor(color[i])
    a.Draw()
    #leg.AddEntry(a,legname[i],"p")
    #leg.Draw()
    c.SetTicks(1,1)
    c.Print("%s%s.png" %(outDir+"pngs/", outputName))
    c.Print("%s%s.pdf" %(outDir, outputName))
    c.Clear()

#c.Update() 
#graph = a.GetPaintedGraph()
#graph.GetXaxis().SetLimits(0, 2000)
#graph.SetMinimum(0)
#graph.SetMaximum(1.4)
#c.Update()

################################

timer.Stop()
timer.Print()