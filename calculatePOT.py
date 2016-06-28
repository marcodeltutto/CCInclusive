#!/usr/bin/env python
#
# read a root anatree using Pyroot
#
# run with python calculatePOT.py

import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain
#from ROOT import *
from array import array
from glob import glob

# Opening root file
#fname = glob("/pnfs/uboone/persistent/users/aschu/MC_BNB_Cosmic/prodgenie_bnb_nu_cosmic_uboone_v05_08_00_anatree.root")
#fname = glob("/pnfs/uboone/persistent/users/aschu/MEC/MECmerge.root")
#fname = glob("/pnfs/uboone/persistent/users/aschu/TEM/TEMmerge.root")
#fname = glob("/pnfs/uboone/scratch/users/mdeltutt/v05_11_00/anatree_bnb_cosmic_eventWeight_MACCQE_ccincl_wflash_2/*/ana_hist.root") # Ma CCQE Selection I
fname = glob("/pnfs/uboone/scratch/users/mdeltutt/v05_11_00/anatree_prodgenie_bnb_nu_cosmic_uboone_mcc7_reco2_reweightMACCQE/*/ana_hist.root") # Ma CCQE Selection II
#print "Input file: ", fname

# Creating TChain
chain    = TChain("analysistree/anatree")
chainPOT = TChain("analysistree/pottree")
for f in fname: 
  chain.Add(f)
  chainPOT.Add(f)
           
# Printing the number of entries
entries    = chain.GetEntries()
entriesPOT = chainPOT.GetEntries()
print "Number of entries in the ANA tree: ", entries
print "Number of entries in the POT tree: ", entriesPOT

# Getting POT number
#NominalPOT = 6.0e20
TotalPOT = 0
for jentry in xrange( entriesPOT ):
  entry = chainPOT.GetEntry( jentry )
  #print '{:e}'.format(float(chainPOT.pot))
  TotalPOT += chainPOT.pot 
print "Accumulated POT: ", '{:e}'.format(float(TotalPOT))
#print "Nominal POT:     ", '{:e}'.format(float(NominalPOT))



raw_input("Please press enter to exit.")



