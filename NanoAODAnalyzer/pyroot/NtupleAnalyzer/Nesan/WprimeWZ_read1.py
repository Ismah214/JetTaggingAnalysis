import sys
import os
import glob
import ROOT
import plotter
import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-s", "--sample",      action="store",   dest="sample",      default="ZprimeToTTM" )
parser.add_option("-wb", "--test", action="store_true",   dest="test",      default= False )

(options, args) = parser.parse_args()

TEST = options.test
SAMPLE = options.sample

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
if SAMPLE == "WprimeToWZM":
  arg1       = ["600"] if TEST else ["600","800","1000","1200","1400","1600","2000","2500","3000","3500","4000","4500","5000","5500","6000","6500","7000","7500","8000"]
if SAMPLE == "QCDPt15To7000":
  arg1 = ["_part1"] if TEST else ["_part"+str(x) for x in range(0,81)]
p           = [arg1]

### Define histograms ###
bins = 5000 

### pt vs taggers histos ###
h_pt_deepAK8    = ROOT.TH2F("h_pt_deepAK8", "Top-tagging DeepAK8 vs p_{T}; p_{T}; DeepAK8", 50, 0., 5000., 2000, 0., 1.)
h_pt_tau21      = ROOT.TH2F("h_pt_tau21", "Top-tagging #tau_{21} vs p_{T}; p_{T}; #tau_{21}", 50, 0., 5000., 2000, 0., 1.)

### pt vs taggers histos QCD ###
h_pt_W_deepAK8    = ROOT.TH2F("h_pt_W_deepAK8", "Top-tagging DeepAK8 vs p_{T}; p_{T}; DeepAK8", 50, 0., 5000., 2000, 0., 1.)
h_pt_W_tau21      = ROOT.TH2F("h_pt_W_tau21", "Top-tagging #tau_{21} vs p_{T}; p_{T}; #tau_{21}", 50, 0., 5000., 2000, 0., 1.)

### make sure no same part in the list ###
def checkList(gp, gpfinal):
    for p in gpfinal:
        if getFinal(gp)._index == p._index:
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
    
    if "Wprime" in SAMPLE:
        
        INDIR2  =  INDIR
        INDIR2  += "_ext1"
        inFileListExt = [EOSUSER+f for f in glob.glob(INDIR2+".root")]
        
        for inFile in inFileListExt:
            tree.Add(inFile)  
        
    ### USE TCHAIN AND SETUP TTREEREADER ###
    inTree  = InputTree(tree)
    nEvents = inTree.GetEntries()
    
    ### LOOP OVER EVENTS ###
    print "Total events: ", nEvents
    #xrange returns xrange objects so list operations cannot be performed on them. range returns list but is slower than xrange.
    for i in xrange(0,nEvents):     
        if i%10000==0: 
          print "Running %d out of %d" %(i, nEvents)
        
        evt       = Event(inTree,i)
        if "Wprime" in SAMPLE:
            particles = Collection(evt, "GenPart")
            fatjet    = Collection(evt, "FatJet")
            jetak8    = Collection(evt, "GenJetAK8")
            
            ###############################################
            ### efficiency with respect to top quark pt ### 
            ###############################################
            
            AssignDauIdx(particles)
            
            ak8jet  = [x for x in fatjet if (x.pt>200. and abs(x.eta)<2.4)]
            nW      = [x for x in particles if (abs(x.pdgId)==24 and abs(x.eta)<2.4)]
            
            if len(nW)<1: 
                continue 
                #continue to next loop if no W's in this set
            
            Wfinal = []
            
            for wb in nW:
                if checkList(wb, Wfinal) and decayTo2Q(wb): 
                #if wb is a new entry into Wfinal, and it decays to 2 quarks
                    Wfinal.append(getFinal(wb)) 
                    #then add this final decaying W into Wfinal
            
            if len(Wfinal)<1:
            #continue to next loop if no candidate W's in this set 
                continue 
            
            for wb in Wfinal:
                if inAK8(wb):
                    for fj in ak8jet:
                    #fill the histos with the data #discriminator value of deepAK8 of top from QCD vs pT of the top 
                        h_pt_deepAK8.Fill(fj.pt, fj.deepTag_WvsQCD)
                        #n-subjettiness value is the discriminator value of the top from QCD; vs pT of the top            
                        h_pt_tau21.Fill(fj.pt, fj.tau21)
            
        if "QCD" in SAMPLE: 
        
            fatjet  = Collection(evt, "FatJet")
            ak8jet  =   [x for x in fatjet if (x.pt>200. and abs(x.eta)<2.4)]
        
            for fj in ak8jet:
                #fill the histos with the data #discriminator value of deepAK8 of top from QCD vs pT of the top 
                h_pt_W_deepAK8.Fill(fj.pt, fj.deepTag_WvsQCD)
                #n-subjettiness value is the discriminator value of the top from QCD; vs pT of the top            
                h_pt_W_tau21.Fill(fj.pt, fj.tau21)
        
        
color = [46,42,30,36,38]
lgd = [ROOT.kFullSquare, ROOT.kFullCircle, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kOpenCircle, ROOT.kOpenSquare]

#################
### ROOT FILE ###
#################
filename = SAMPLE
filename += "_W_tagger_pt"
outFile = ROOT.TFile("%s%s.root" %(outDir, filename), "RECREATE") #/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/WprimeToWZM_W_tagger_pt.root

if "QCD" in SAMPLE:
    for a in [h_pt_W_deepAK8,h_pt_W_tau21]:
        a.SetStats(0)
        a.SetOption("COLZ")
        a.Write()
    outFile.Close()

if "Wprime" in SAMPLE:
    for a in [h_pt_deepAK8,h_pt_tau21]:
        a.SetStats(0)
        a.SetOption("COLZ")
        a.Write()
    outFile.Close()

################
### PDF FILE ###
################
pdfname = SAMPLE
pdfname += "_eff_pt"

## CANVAS ##
legname=["#bf{#DeltaR(AK8,wb)<0.8}", "#bf{max#DeltaR(AK8,q)<0.8}"]
leg=BeginLeg(len(legname))
if "Wprime" in SAMPLE:
    i=0
    c = ROOT.TCanvas("","",877, 620)
    for a in [h_pt_deepAK8,h_pt_tau21]:
        a.SetStats(0)
        #a.SetStatisticOption(0)
        a.SetLineWidth(1)
        a.SetLineColor(color[i])
        a.SetLineColor(color[i])
        a.SetMarkerColor(color[i])
        a.SetMarkerStyle(lgd[i])
        a.SetFillColor(color[i])
  
    if i==0:
        a.Draw()
        leg.AddEntry(a,legname[i],"p")
    else:
        a.Draw("SAME")
        leg.AddEntry(a, legname[i], "p")
    i+=1

    leg.Draw()
    c.SetTicks(1,1)
    c.Print("%s%s.pdf" %(outDir, pdfname))
    
if "QCD" in SAMPLE:
    j=0
    c = ROOT.TCanvas("","",877, 620)
    for a in [h_pt_W_deepAK8,h_pt_W_tau21]:
        a.SetStats(0)
        #a.SetStatisticOption(0)
        a.SetLineWidth(1)
        a.SetLineColor(color[j])
        a.SetLineColor(color[j])
        a.SetMarkerColor(color[j])
        a.SetMarkerStyle(lgd[j])
        a.SetFillColor(color[j])

    if j==0:
        a.Draw()
        leg.AddEntry(a,legname[j],"p")
    else:
        a.Draw("SAME")
        leg.AddEntry(a, legname[j], "p")
    j+=1

    leg.Draw()
    c.SetTicks(1,1)
    c.Print("%s%s.pdf" %(outDir, pdfname))
  #c.Update() 
  #graph = a.GetPaintedGraph()
  #graph.GetXaxis().SetLimits(0, 2000)
  #graph.SetMinimum(0)
  #graph.SetMaximum(1.4)
  #c.Update()

################################

timer.Stop()
timer.Print()