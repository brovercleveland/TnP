#! /usr/bin/env python

import ROOT
import numpy as np

ROOT.gROOT.ProcessLine('.L ./tdrstyle.C')
ROOT.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

can= ROOT.TCanvas('can','canvas',800,600)

oldPtBins = [12.5,17.5,25.0,35.0,45.0,70.0]
newPtBins = [15.0,30.0,50.0,80.0,150.0]

oldMCEta1 = [0.623,0.78,0.90, 0.95, 0.97, 0.98]
oldMCEta2 = [0.623,0.78,0.92, 0.96, 0.98, 0.98]
oldMCEta3 = [0.57,0.75,0.79, 0.89, 0.95, 0.95]
oldMCEta4 = [0.57,0.75,0.90, 0.95, 0.97, 0.97]
oldMCEta5 = [0.57,0.75,0.87, 0.92, 0.94, 0.95]

newMCEta1 = [0.76,0.92,0.94, 0.96, 0.97]
newMCEta2 = [0.74,0.92,0.96, 0.97, 0.97]
newMCEta3 = [0.63,0.82,0.92, 0.89, 0.94]
newMCEta4 = [0.64,0.88,0.95, 0.96, 0.97]
newMCEta5 = [0.64,0.85,0.89, 0.91, 0.93]

oldMCGraph1 = ROOT.TGraph(len(oldMCEta1),np.array(oldPtBins),np.array(oldMCEta1))
newMCGraph1 = ROOT.TGraph(len(newMCEta1),np.array(newPtBins),np.array(newMCEta1))

can.cd()
tmpHist = ROOT.TH2F('tmp','',10,9,200,10,0.5,1)

tmpHist.Draw()
newMCGraph1.Draw('LP')
oldMCGraph1.Draw('LP')
can.SaveAs('test.pdf')

