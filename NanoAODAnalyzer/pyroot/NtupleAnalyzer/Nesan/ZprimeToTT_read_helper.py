import sys
import os
import glob
import ROOT
import plotter

from collections import OrderedDict 
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Event
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import InputTree

ROOT.gROOT.SetBatch()

#### check if gp -> W -> qq in AK8 ###
#def decayToWAK8(gp):
#        if len(gp.dauIdx) == 0:
#            raise ValueError('Particle has no daughters!')
#        for idx in gp.dauIdx:
#          dp=particles[idx]
#          dpart=getFinal(dp)
#          if abs(dpart.pdgId) == 24 and inAK8(dpart) and decayTo2QAK8(dpart):
#            return True
#        return False

#### check if gp -> b in AK8 ###
#def decayToBAK8(gp):
#  if len(gp.dauIdx) == 0:
#      raise ValueError('Particle has no daughters!')
#  for idx in gp.dauIdx:
#    dpart=particles[idx]
#    if abs(dpart.pdgId) == 5 and inAK8(dpart):
#      return True
#  return False

### check if gp decays to at least 2 quarks in AK8 ###
#def decayTo2QAK8(gp):
#        if len(gp.dauIdx) == 0:
#            raise ValueError('Particle has no daughters!')
#        for idx in gp.dauIdx:
#          dpart=particles[idx]
#          if abs(dpart.pdgId) < 6 and inAK8(dpart):            
#            gp.dauIdx.reverse()
#            for idx in gp.dauIdx:
#              dpart=particles[idx]
#              if abs(dpart.pdgId)<6 and inAK8(dpart):
#                return True
#        return False

#looping over all particles and storing how many quark daughters each particle produces
def isHadronic(gp):
  dautemp=[x for x in gp.dauIdx if abs(particles[x].pdgId)<6 and abs(particles[x].eta)<2.4] 
  if len(dautemp)==2: #fullhadron
    return 0
  elif len(dautemp)==1: #semilept
    return 1
  elif len(dautemp)==0: #fulllept
    return 2

def BeginLeg(ll):
    ## LEGEND POSITION ##
    legpos = "Right"
    if legpos == "Left":
      xLat = 0.13
    elif legpos == "Right":
      xLat = 0.65 # move left and right
    else:
      xLat = 0.2
    yLat = 0.85 # move legend up and down, the larger( about 1 already at the edge of canvas) the higher it goes
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

## define dauIdx for the mother particles ###
def AssignDauIdx(gp):   
    for idx, dp in enumerate(gp):
      if not hasattr(dp, 'dauIdx'): 
        dp.dauIdx = []

      if dp.genPartIdxMother >= 0: #access particles with a mother Id i.e. all daughter particles of some mother
        mom = gp[dp.genPartIdxMother] #the mother particle is accessed using the genpartIdxMother and stored in mom
        if not hasattr(mom, 'dauIdx'): 
            mom.dauIdx = [idx] #feed the index of the daughter to the mother particle

        else: 
            mom.dauIdx.append(idx)