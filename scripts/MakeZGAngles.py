#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import *
import numpy as np


def ComparisonSuite():
  FileP6= TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/simpleHistograms_ggM125_8TeV_v2.root")
  FileP8_153= TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/simpleHistograms_ggM125_8TeV_pythia8_v2.root")

  #selection = ['el','mu']
  plotList = []
  fileList = [FileP6,FileP8_153]
  labelList = ['pythia6','pythia8_153']
  colorList = [kMagenta+1,kRed]
  physList = ['costheta_lm','costheta_lp','phi','cosTheta']
  can= TCanvas('can','canvas',800,600)
  can.SetGrid(2,2)
  gStyle.SetOptStat(0)
  for phys in physList:
    can.cd()
    leg = TLegend(0.83,0.75,0.99,0.92,'',"brNDC")
    leg.SetBorderSize(1)
    leg.SetTextSize(0.03)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)
    for i,thisFile in enumerate(fileList):
      plotList.append(thisFile.Get(phys))
      plotList[-1].SetLineColor(colorList[i])
      plotList[-1].SetLineWidth(2)
      if i == 1: plotList[-1].SetLineStyle(2)
      #if thisLeaf.Draw(phys)
      plotList[-1].Scale(1/plotList[-1].Integral())
      leg.AddEntry(plotList[-1],labelList[i],'l')
    ymax = max(map(lambda x:x.GetMaximum(),plotList[-2:]))*1.1
    #ymin = max(map(lambda x:x.GetMinimum(),plotList[-2:]))*0.9
    plotList[-2].Draw()
    plotList[-2].SetTitle(phys)
    plotList[-2].GetXaxis().SetTitle('')
    plotList[-2].GetYaxis().SetTitle('')
    plotList[-2].SetMaximum(ymax)
    plotList[-1].Draw('same')
    leg.Draw()
    can.Update()
    #raw_input()
    can.SaveAs('Pythia6vs8/'+phys+'_CMSSW.pdf')

if __name__=="__main__":
  ComparisonSuite()
