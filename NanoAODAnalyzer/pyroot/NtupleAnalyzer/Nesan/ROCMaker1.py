import os
import ROOT
import collections
import math
import sys
from ROOT import TGraph
import array as array
  
ROOT.gROOT.SetBatch()

color = [46,42,41,30,36,38]
linestyle = [1,3,5,7,9,11]
EOSUSER = "root://eosuser.cern.ch/"
outDir = "/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/ROC/"
outName = [
("W_ROC_DeepAK8",
"W_ROC_tau21"),
("T_ROC_DeepAK8",
"T_ROC_tau32"),
]

def main():
  inFileName_sig = [
  EOSUSER+"/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/WprimeWZM_tagger_pt.root",
  EOSUSER+"/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/ZprimeToTTM_tagger_pt.root",
  ]
  
  inFileName_bg  = EOSUSER+"/eos/user/n/nsubrama/AnaJetTagging/ANALYSIS/tagger_pt/QCDPt15To7000_tagger_pt.root"

  part = ["T", "W"]
  taggers = [
  ["deepAK8"],
  ["tau21"], 
  ["tau32"],
  ]

  legname = [
  "300 < p_{T}^{truth} < 500 GeV", 
  "500 < p_{T}^{truth} < 1000 GeV",
  "1000 < p_{T}^{truth} < 1500 GeV",
  "1500 < p_{T}^{truth} < 2000 GeV",
  "2000 < p_{T}^{truth} < 2500 GeV",
  ]
  
  ProjCut=[(1,1000),(1001,2000),(2001,3000),(3001,4000),(4001,5000)]

  g = [None]*len(ProjCut)
  h_sig=[None]*len(ProjCut)
  h_bg=[None]*len(ProjCut)

  yaxisname = "background mistag"
  xaxisname = "signal efficiency"
  titlename = [
  ("ROC Curve for W-tagging (DeepAK8)",
  "ROC Curve for W-tagging (#tau_{21})"),
  ("ROC Curve for T-tagging (DeepAK8)",
  "ROC Curve for T-tagging (#tau_{32})"),
  ]

  #user chooses which ROC curve to make
  a1=input("0.T 1.W :: ")

  # Get root file with histograms
  inputFile_sig = ROOT.TFile.Open(inFileName_sig[a1])
  inputFile_bg  = ROOT.TFile.Open(inFileName_bg)

  # Make TCanvas
  c = ROOT.TCanvas("","", 877, 620)
  c.SetFillStyle(4000)
  c.SetLeftMargin(0.15)
  c.SetRightMargin(0.08)
  c.SetTopMargin(0.1)
  c.SetBottomMargin(0.15)
  # c.SetGrid()
  c.SetLogy()

  ### LEGENDS ###
  def BeginLeg(ll):
    ## LEGEND POSITION ##
    legpos = "Right"
    if legpos == "Left":
      xLat = 0.13
    elif legpos == "Right":
      xLat = 0.65 # move left and right
    else:
      xLat = 0.2
    yLat = 0.4 # move legend up and down, the larger( about 1 already at the edge of canvas) the higher it goes
    xLeg = xLat
    yLeg = yLat

    leg_g =  0.04 * ll # num of legend entries
    leg = ROOT.TLegend( xLeg+0.05, yLeg - leg_g, xLeg+0.15, yLeg )
    leg.SetNColumns(1)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(14)
    leg.SetBorderSize(0)
    return leg

  ########################################################
  #################################
  ### EACH PROJCUT IN SAME HIST ###
  ### WITH DIFFERENT TAGGER     ###
  #################################
  leg = BeginLeg(len(legname))

  ### prompt input to choose tag and taggers ###
  a3 = input("0.deeptag, 1.tau :: ")
  if a3 == 0:
    a2=0 #deepTag
  elif a1 == 0 and a3 == 1:
    a2=1 #tau21
  elif a1 == 1 and a3 == 1:
    a2=2 #tau32

  for j, tag in enumerate(taggers[a2]):
    print "running %s" %(tag)
    
    # get 2D histograms from root file
    if tag == "deepAK8": ### deeptag
      h2d_sig = inputFile_sig.Get("h_pt_%s" %(tag)) #ZprimeToTTM_tagger_pt.root get h_pt_deepAK8
      h2d_bg  = inputFile_bg.Get("h_pt_%s_%s" %(part[a1], tag)) #QCDPt15To7000_tagger_pt.root get h_pt_T_deepAK8
  
    else: ### taus
      h2d_sig = inputFile_sig.Get("h_pt_%s" %(tag)) #ZprimeToTTM_tagger_pt.root get h_pt_tau32
      h2d_bg  = inputFile_bg.Get("h_pt_%s_%s" %(part[a1], tag)) #QCDPt15To7000_tagger_pt.root get h_pt_T_tau32
  
    # 1D projections
    for i, p in enumerate(ProjCut):
      a = p[0] #4
      b = p[1] #6
      h_sig[i] = h2d_sig.ProjectionY("py_sig_%dTo%dBins_"%(a,b) +tag , a, b) #Projects 2D histo to 1D histo along Y
      h_bg[i]  = h2d_bg.ProjectionY("py_bg_%dTo%dBins_"%(a,b) +tag, a, b)
      
      bins = 4
      h_sig[i].Rebin(bins)
      h_bg[i].Rebin(bins)

      fbin = h_sig[i].FindFirstBinAbove()
      nbin = h_sig[i].FindLastBinAbove()
    
      # define arrays #
      x= array.array('d')  
      y= array.array('d')

      for ibin in xrange(fbin, nbin+1):
        if tag == "deepAK8":
          # signal efficiency
          total_sig = h_sig[i].Integral(fbin, nbin+1)
          cut_sig   = h_sig[i].Integral(ibin, nbin+1) #increasing tagger value
          eff = (cut_sig/total_sig)
          x.append(eff)
        
          # background mistag/acceptance
          total_bg = h_bg[i].Integral(fbin, nbin+1)
          cut_bg   = h_bg[i].Integral(ibin, nbin+1)
          misid = (cut_bg/total_bg)
          y.append(misid)

        else: #taus
          # signal efficiency
          total_sig = h_sig[i].Integral(fbin, nbin+1)
          cut_sig   = h_sig[i].Integral(fbin, ibin)
          eff = (cut_sig/total_sig)
          x.append(eff)
        
          # background mistag
          total_bg = h_bg[i].Integral(fbin, nbin+1)
          cut_bg   = h_bg[i].Integral(fbin, ibin)
          misid = (cut_bg/total_bg)
          y.append(misid)

      n = nbin +1 - fbin
      g[i] = TGraph(n,x,y)
      g[i].SetTitle("%s; %s; %s" %(titlename[a1][a3] ,xaxisname, yaxisname))
      g[i].GetXaxis().SetRangeUser(0,1.1)
      g[i].GetYaxis().SetRangeUser(10e-5,5)
      g[i].SetLineStyle(linestyle[i])
      g[i].SetLineWidth(2)
      g[i].SetLineColor(color[i])
      g[i].SetMarkerColor(color[i])
      leg.AddEntry(g[i], legname[i])
      if i==0:
        g[i].Draw("AC")
      else:
        g[i].Draw("C SAME")

  leg.Draw()
  c.SetTicks(1,1)
  c.Print("%s%s_ROC_%s_test.pdf" %(outDir, part[a1], taggers[a2][0]))
  h2d1.RebinY(20)
  h2d2.RebinY(20)
  h2d1.Write("pt_%s" %(taggers[0]))
  h2d2.Write("pt_%s" %(taggers[1]))
  outFile.Close()

  ########################################################
if __name__ == "__main__":
  main()