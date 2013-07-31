#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import *
import numpy as np

gROOT.ProcessLine('.L ./tdrstyle.C')
setTDRStyle()

def ComparisonSuite():
  File1= TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/localHistos/higgsHistograms_ggM125_8TeV_local.root")
  File2= TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/localHistos/higgsHistograms_ggM125_8TeV_pythia8_local.root")

  #selection = ['el','mu']
  plotList = []
  fileList = [File1,File2]
  labelList = ['pythia6','pythia8']
  colorList = [kRed,kBlue]
  selection = ['mu']
  catList = ['CAT1','CAT2','CAT3','CAT4','','PostGen','ZGAngles']
  distList = ['Pt','Eta','Phi','Mass']
  physList1 = ['Photon','LeadingLepton','TrailingLepton','DiLep','ThreeBody']
  physList2 = ['costhetaLP','costhetaLM','phi','costhetaZG']
  can= TCanvas('trailingCan','canvas',800,600)
  can.cd()
  can.SetGrid(2,2)
  gStyle.SetOptStat(0)
  for sel in selection:
    for cat in catList:
      folder = cat
      if cat == '': folder = 'pT-Eta-Phi'
      if cat == 'ZGAngles': physListTmp = physList2
      else: physListTmp = physList1
      for phys in physListTmp:
        if phys in physList2: distListTmp = ['']
        elif phys is 'DiLep' or phys is 'ThreeBody': distListTmp = distList
        else: distListTmp = distList[0:3]

        for dist in distListTmp:
          for i,thisFile in enumerate(fileList):
            if i == 0: suffix = ''
            elif i ==1: suffix = '_pythia8'
            if sel == 'mu':
              if phys == 'Photon':
                thisPlot = thisFile.GetDirectory(folder).Get('h1_photon'+dist+cat+'_Signal2012ggM125'+suffix)
              if phys == 'LeadingLepton':
                thisPlot = thisFile.GetDirectory(folder).Get('h1_leadLeptonStd'+dist+cat+'_Signal2012ggM125'+suffix)
              if phys == 'TrailingLepton':
                thisPlot = thisFile.GetDirectory(folder).Get('h1_trailingLeptonStd'+dist+cat+'_Signal2012ggM125'+suffix)
              if phys == 'DiLep':
                thisPlot = thisFile.GetDirectory(folder).Get('h1_diLep'+dist+cat+'_Signal2012ggM125'+suffix)
              if phys == 'ThreeBody':
                thisPlot = thisFile.GetDirectory(folder).Get('h1_threeBody'+dist+cat+'_Signal2012ggM125'+suffix)
              elif phys not in physList1:
                thisPlot = thisFile.GetDirectory(folder).Get('h1_'+phys+'_Signal2012ggM125'+suffix)



            print sel, cat, phys, dist, thisFile
            thisPlot.Scale(1/thisPlot.Integral())
            thisPlot.SetLineColor(colorList[i])
            thisPlot.SetLineWidth(2)
            thisPlot.SetMarkerColor(colorList[i])
            thisPlot.SetMarkerStyle(20+i)
            if i == 0:
              thisPlot.GetYaxis().SetTitleSize(0.05)
              thisPlot.GetYaxis().CenterTitle()
              thisPlot.GetXaxis().SetTitleSize(0.05)
              thisPlot.GetYaxis().SetLabelSize(0.05)
              thisPlot.GetXaxis().SetLabelSize(0.05)
              thisPlot.GetYaxis().SetTitle("Events")
              #thisPlot.GetXaxis().SetTitle(dist)
              thisPlot.SetTitle(sel+sel+' '+phys+' '+cat+' '+dist)
              thisPlot.SetTitleOffset(1.2,"Y")

              leg = TLegend(0.85,0.75,0.95,0.92,'',"brNDC")
              leg.SetBorderSize(1)
              leg.SetTextSize(0.03)
              leg.SetFillColor(0)
              leg.SetFillStyle(1001)
            plotList.append(thisPlot)

          ymax = max(map(lambda x:x.GetMaximum(),plotList))*1.1
          ymin = max(map(lambda x:x.GetMinimum(),plotList))*0.9
          for i,plot in enumerate(plotList):
            if i == 0:
              plot.SetMaximum(ymax)
              plot.SetMinimum(ymin)
              plot.Draw('hist')
            else:
              plot.Draw('histsame')
            leg.AddEntry(plot,labelList[i],'l')
          leg.Draw()
          can.SaveAs('Pythia6vs8/'+sel+'_'+dist+'_'+cat+'_'+phys+'.pdf')
          plotList = []

if __name__=="__main__":
  ComparisonSuite()
