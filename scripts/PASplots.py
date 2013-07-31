#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import *
import numpy as np

gROOT.ProcessLine('.L ./tdrstyle.C')
setTDRStyle()

def PtRatios():
  el2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_MoriondV5_2-10-13.root')
  mu2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_MoriondV6_2-20-13.root')
  leptons = ['mu','el']
  objects = ['photon','leadLeptonStd','trailingLeptonStd', 'energy']

  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(2.7)
  gStyle.SetTitleH(0.08) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(0.99)    # Set the position of the title box
  gStyle.SetLineWidth(2)

  canvas1 = TCanvas('canvas1','canvas1',900,700)
  canvas2 = TCanvas('canvas2','canvas2',900,700)
  pad1 = TPad('pad1', '', 0.01, 0.34, 0.95, 0.98, 0)
  pad1.SetBottomMargin(0.07)
#  pad1.SetGrid(1,1)
  pad2 = TPad('pad2', '', 0.01, 0.02, 0.95, 0.36, 0)
  pad2.SetTopMargin(0.)
  pad2.SetBottomMargin(0.2)
  pad2.SetGrid(2,2)

  dict2012_Moriond = {'mu':mu2012_Moriond,'el':el2012_Moriond}
  legend = TLegend(0.6,0.7,0.85,0.85,'',"brNDC")
  legend.SetBorderSize(0)
  legend.SetTextSize(0.025)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  for lep in leptons:
    for obj in objects:
      if obj == 'energy':
        pt125 = dict2012_Moriond[lep].GetDirectory('ZGamma').Get('h1_GammaEnergy_Signal2012ggM125')
        pt120 = dict2012_Moriond[lep].GetDirectory('ZGamma').Get('h1_GammaEnergy_Signal2012ggM120')
      else:
        pt125 = dict2012_Moriond[lep].GetDirectory('pT-Eta-Phi').Get('h1_'+obj+'Pt_Signal2012ggM125')
        pt120 = dict2012_Moriond[lep].GetDirectory('pT-Eta-Phi').Get('h1_'+obj+'Pt_Signal2012ggM120')


      if obj == 'energy':
        pt125.SetName('photon energy, M125')
      else:
        pt125.SetName(obj+' pt, M125')
      pt125.SetLineColor(kBlack)
      pt125.GetYaxis().CenterTitle()
      pt125.GetYaxis().SetTitle('A.U.')
      pt125.Scale(1/pt125.Integral())
      if obj == 'energy':
        pt120.SetName('photon energy, M120')
      else:
        pt120.SetName(obj+' pt, M120')
      pt120.SetLineColor(kRed)
      pt120.GetYaxis().CenterTitle()
      pt120.GetYaxis().SetTitle('A.U.')
      pt120.Scale(1/pt120.Integral())

      if (lep == 'el'):
        pt120.SetTitle(pt120.GetTitle()+' Electron')
        pt125.SetTitle(pt125.GetTitle()+' Electron')
      else:
        pt120.SetTitle(pt120.GetTitle()+' Muon')
        pt125.SetTitle(pt125.GetTitle()+' Muon')

      canvas1.cd()
      pad1.Draw()
      pad2.Draw()
      pad1.cd()
      ratio = pt120.Clone()
      ratio.Divide(pt125)
      ratio.SetTitle('')
      ratio.SetLineColor(kBlack)
      ratio.GetYaxis().SetTitle('M120/M125')
      ratio.GetXaxis().SetTitleSize(0.055)
      ratio.GetXaxis().SetLabelSize(0.05)
      ratio.GetYaxis().SetTitleSize(0.055)
      ratio.GetYaxis().SetLabelSize(0.05)
      #ratio.SetMinimum(0)
      #ratio.SetMaximum(15)
      pt120.Draw('hist')
      pt125.Draw('histsame')
      if obj == 'energy':
        legend.AddEntry(pt125, 'photon energy M125', 'l')
        legend.AddEntry(pt120, 'photon energy M120', 'l')
      else:
        legend.AddEntry(pt125, obj+' pt M125', 'l')
        legend.AddEntry(pt120, obj+' pt M120', 'l')
      legend.Draw()
      pad2.cd()
      ratio.Draw()
      if obj == 'energy':
        canvas1.SaveAs('PhotonEnergyRatio_'+lep+'.pdf')
      else:
        canvas1.SaveAs(obj+'PtRatio_'+lep+'.pdf')
      legend.Clear()
      pad1.Clear()
      pad2.Clear()


def RMSandFWHMandYieldTable():
  el2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_MoriondV5_2-10-13.root')
  mu2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_MoriondV6_2-20-13.root')
  el2011 = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2011_Medium_ICHEPmzmzg_1-19-13.root')
  mu2011 = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_ICHEPmzmzg_1-19-13.root')
  catList = ['CAT1','CAT2','CAT3','CAT4']
  leptons = ['mu','el']
  dict2012_Moriond = {'mu':mu2012_Moriond,'el':el2012_Moriond}
  dict2011 = {'mu':mu2011,'el':el2011}
  yearDict = {2011:dict2011,2012:dict2012_Moriond}

  canvas1 = TCanvas('canvas1','canvas1',900,700)
  canvas1.cd()
  sigByCatMu2011 = []
  sigByCatMu2012 = []
  sigByCatEl2011 = []
  sigByCatEl2012 = []
  significanceByCatMu = {2011:sigByCatMu2011,2012:sigByCatMu2012}
  significanceByCatEl = {2011:sigByCatEl2011,2012:sigByCatEl2012}
  sigYieldTot = 0
  dataYieldTot = 0
  for year in yearDict.keys():
    for lep in leptons:
      for cat in catList:
        m125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'ggM125')
        m125.Fit('gaus','Q','',121,135)
        #m125.Draw()
        #canvas1.SaveAs('signal_'+lep+'_'+cat+'.pdf')
        sigma = m125.GetFunction('gaus').GetParameter(2)
        fwhm = sigma*2.35482
        print year, lep, cat, 'sigma', '%.2f'%sigma, 'FWHM', '%.2f'%fwhm

        gg125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'ggM125')
        vbf125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'vbfM125')
        wh125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'whM125')
        zh125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'zhM125')
        tth125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal'+str(year)+'tthM125')
        if year == 2012:
          if lep == 'el':
            lumi = 19.6195
          else:
            lumi = 19.6175
        else:
          if lep == 'el':
            lumi = 4.98
          else:
            lumi = 5.05
        gg125.Scale(lumi/(99991/(19.52*0.00154*0.100974*1000)))
        vbf125.Scale(lumi/(99885/(1.578*0.00154*0.10098*1000)))
        wh125.Scale(lumi/(656101/(0.6966*0.00154*1000)))
        zh125.Scale(lumi/(344143/(0.3943*0.00154*1000)))
        tth125.Scale(lumi/(100048/(0.1302*0.00154*0.10098*1000)))
        sigYield = gg125.Integral() + vbf125.Integral() + wh125.Integral() + zh125.Integral() + tth125.Integral()

        data125 = yearDict[year][lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_DATA')
        dataYield = data125.Integral(data125.FindBin(120), data125.FindBin(130))
        dataYieldTot += dataYield
        signif = sigYield/(sigYield+dataYield)
        if lep == 'mu':
          significanceByCatMu[year].append(signif)
        else:
          significanceByCatEl[year].append(signif)

      nEvents_100_180 = yearDict[year][lep].GetDirectory('Misc').Get('h1_acceptanceByCut_DATA').Integral(50,50)
      nEvents_120_150 = yearDict[year][lep].GetDirectory('Misc').Get('h1_acceptanceByCut_DATA').Integral(70,70)

      gg125 = yearDict[year][lep].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal'+str(year)+'ggM125')
      vbf125 = yearDict[year][lep].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal'+str(year)+'vbfM125')
      wh125 = yearDict[year][lep].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal'+str(year)+'whM125')
      zh125 = yearDict[year][lep].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal'+str(year)+'zhM125')
      tth125 = yearDict[year][lep].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal'+str(year)+'tthM125')

      if year == 2012:
        if lep == 'el':
          lumi = 19.6195
        else:
          lumi = 19.6175
      else:
        if lep == 'el':
          lumi = 4.98
        else:
          lumi = 5.05

      gg125.Scale(lumi/(99991/(19.52*0.00154*0.100974*1000)))
      vbf125.Scale(lumi/(99885/(1.578*0.00154*0.10098*1000)))
      wh125.Scale(lumi/(656101/(0.6966*0.00154*1000)))
      zh125.Scale(lumi/(344143/(0.3943*0.00154*1000)))
      tth125.Scale(lumi/(100048/(0.1302*0.00154*0.10098*1000)))

      sigYield = gg125.Integral() + vbf125.Integral() + wh125.Integral() + zh125.Integral() + tth125.Integral()
      sigYieldTot += sigYield

      print year
      print 'gg'.ljust(10), 'vbf'.ljust(10), 'wh'.ljust(10), 'zh'.ljust(10), 'tth'.ljust(10), 'Total'.ljust(10)
      print '{0:<{1}.5}'.format(gg125.Integral(),10), '{0:<{1}.5}'.format(vbf125.Integral(),10), '{0:<{1}.5}'.format(wh125.Integral(),10), '{0:<{1}.5}'.format(zh125.Integral(),10),'{0:<{1}.5}'.format(tth125.Integral(),10), '{0:<{1}.5}'.format(sigYield,10)
      print lep, 'lumi', lumi, '100-180', '%i'%nEvents_100_180, '120-150', '%i'%nEvents_120_150, 'M125', '%.1f'%sigYield
    print
    print year, 'muon signif', zip(catList,significanceByCatMu[year])
    print year, 'elec signif', zip(catList,significanceByCatEl[year])
  print 'total signal (125):', sigYieldTot
  print 'total data (120-130):', dataYieldTot



def BGComp():
  el2011 = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2011_Medium_ICHEPmzmzg_1-19-13.root')
  mu2011 = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_ICHEPmzmzg_1-19-13.root')
  el2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_MoriondV4_2-7-13.root')
  mu2012_Moriond = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_MoriondV4_2-7-13.root')
  el2012_ICHEP = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_ICHEPmzmzg_1-28-13.root')
  mu2012_ICHEP = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_ICHEPmzmzg_1-28-13.root')

  catList = ['CAT1','CAT2','CAT3','CAT4']
  leptons = ['mu','el']
  dict2011 = {'mu':mu2011,'el':el2011}
  dict2012_Moriond = {'mu':mu2012_Moriond,'el':el2012_Moriond}
  dict2012_ICHEP = {'mu':mu2012_ICHEP,'el':el2012_ICHEP}

  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(2.7)
  gStyle.SetTitleH(0.08) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(0.99)    # Set the position of the title box
  gStyle.SetLineWidth(2)

  canvas1 = TCanvas('canvas1','canvas1',900,700)
  canvas2 = TCanvas('canvas2','canvas2',900,700)
  pad1 = TPad('pad1', '', 0.01, 0.34, 0.95, 0.98, 0)
  pad1.SetBottomMargin(0.07)
#  pad1.SetGrid(1,1)
  pad2 = TPad('pad2', '', 0.01, 0.02, 0.95, 0.36, 0)
  pad2.SetTopMargin(0.)
  pad2.SetBottomMargin(0.2)
  pad2.SetGrid(2,2)

  legend = TLegend(0.6,0.7,0.85,0.85,'',"brNDC")
  legend.SetBorderSize(0)
  legend.SetTextSize(0.025)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)

  for lep in leptons:
    for cat in catList:
      mzg_2011 = dict2011[lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_DATA')
      mzg_2012_Moriond = dict2012_Moriond[lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_DATA')
      mzg_2012_ICHEP = dict2012_ICHEP[lep].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_DATA')

      mzg_2011.SetName('2011 data')
      mzg_2011.SetLineColor(kBlue)
      mzg_2012_Moriond.SetName('2012 data, Moriond Selection')
      mzg_2012_Moriond.SetLineColor(kRed)
      mzg_2012_Moriond.GetYaxis().CenterTitle()
      mzg_2012_Moriond.GetYaxis().SetTitle('Events/GeV')
      mzg_2012_ICHEP.SetName('2012 data, ICHEP Selection')
      mzg_2012_ICHEP.SetLineColor(kGreen+1)
      mzg_2012_ICHEP.GetYaxis().CenterTitle()
      mzg_2012_ICHEP.GetYaxis().SetTitle('Events/GeV')

      if (lep == 'el'):
        mzg_2012_ICHEP.SetTitle(mzg_2012_ICHEP.GetTitle()+' Electron, '+cat)
        mzg_2012_Moriond.SetTitle(mzg_2012_Moriond.GetTitle()+' Electron, '+cat)
      else:
        mzg_2012_ICHEP.SetTitle(mzg_2012_ICHEP.GetTitle()+' Muon, '+cat)
        mzg_2012_Moriond.SetTitle(mzg_2012_Moriond.GetTitle()+' Muon, '+cat)



      canvas1.cd()
      pad1.Draw()
      pad2.Draw()
      pad1.cd()
      ratio = mzg_2012_ICHEP.Clone()
      ratio.Divide(mzg_2011)
      ratio.SetTitle('')
      ratio.SetLineColor(kBlack)
      ratio.GetYaxis().SetTitle('ICHEP/2011')
      ratio.GetXaxis().SetTitleSize(0.055)
      ratio.GetXaxis().SetLabelSize(0.05)
      ratio.GetYaxis().SetTitleSize(0.055)
      ratio.GetYaxis().SetLabelSize(0.05)
      ratio.SetMinimum(0)
      ratio.SetMaximum(15)
      mzg_2012_ICHEP.Draw('hist')
      mzg_2011.Draw('histsame')
      legend.AddEntry(mzg_2012_ICHEP, '2012 ICHEP Selection', 'l')
      legend.AddEntry(mzg_2011, '2011 Selection', 'l')
      legend.Draw()
      pad2.cd()
      ratio.Draw()
      canvas1.SaveAs('mzg_2011vs2012ICHEP_'+lep+'_'+cat+'.pdf')
      legend.Clear()
      pad1.Clear()
      pad2.Clear()

      canvas2.cd()
      pad1.Draw()
      pad2.Draw()
      pad1.cd()
      ratio = mzg_2012_Moriond.Clone()
      ratio.Divide(mzg_2012_ICHEP)
      ratio.SetTitle('')
      ratio.SetLineColor(kBlack)
      ratio.GetYaxis().SetTitle('Moriond/ICHEP')
      ratio.GetXaxis().SetTitleSize(0.055)
      ratio.GetXaxis().SetLabelSize(0.05)
      ratio.GetYaxis().SetTitleSize(0.055)
      ratio.GetYaxis().SetLabelSize(0.05)
      ratio.SetMinimum(0)
      ratio.SetMaximum(3)
      mzg_2012_Moriond.Draw('hist')
      mzg_2012_ICHEP.Draw('histsame')
      legend.AddEntry(mzg_2012_Moriond, '2012 Moriond Selection', 'l')
      legend.AddEntry(mzg_2012_ICHEP, '2012 ICHEP Selection', 'l')
      legend.Draw()
      pad2.cd()
      ratio.Draw()
      canvas2.SaveAs('mzg_2012Moriondvs2012ICHEP_'+lep+'_'+cat+'.pdf')
      legend.Clear()
      pad1.Clear()
      pad2.Clear()


def RMSComp():
  mu2012_52X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_Nominal_12-8-12.root')
  mu2012_53X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012CD_signalStudy_1-15-13.root')
  el2012_52X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012_Medium_Nominal_12-8-12.root')
  el2012_53X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012CD_signalStudy_1-15-13.root')

  #mu2012_52X = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_ICHEPnoSF_1-20-13.root')
  #mu2012_53X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_Medium_ICHEP_1-17-13.root')
  #el2012_52X = TFile('/uscms_data/d2/bpollack/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012_Medium_ICHEPnoSF_1-20-13.root')
  #el2012_53X = TFile('/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_Medium_ICHEP_1-17-13.root')

  massList = ['M120','M125','M130','M135','M140','M145','M150']
  catList = ['CAT1','CAT2','CAT3','CAT4']
  leptons = ['mu','el']
  fileDict52X = {'mu':mu2012_52X, 'el':el2012_52X}
  fileDict53X = {'mu':mu2012_53X, 'el':el2012_53X}
  massArray = np.array([120.0,125.0,130.0,135.0,140.0,145.0,150.0])

  gStyle.SetOptStat(1)
  gStyle.SetTitleFontSize(2.7)
  gStyle.SetTitleH(0.08) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(0.99)    # Set the position of the title box
  gStyle.SetLineWidth(2)

  canvas = TCanvas('canvas','canvas',800,600)
  canvas.SetGrid()

  for lepton in leptons:
    for cat in catList:
      rms_52X = []
      rms_53X = []
      for mass in massList:
        if cat == 'CAT1':
          canvas.cd()

          signal_52X = fileDict52X[lepton].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal2012gg'+mass)
          signal_52X.SetName('52X')
          signal_53X = fileDict53X[lepton].GetDirectory('ZGamma').Get('h1_InvariantMass_Signal2012gg'+mass)
          signal_53X.SetName('53X')
          signal_52X.SetLineColor(kBlue)
          signal_53X.SetLineColor(kRed)
          signal_52X.Scale(1/signal_52X.Integral())
          signal_53X.Scale(1/signal_53X.Integral())
          signal_52X.SetTitle('M_{ll#gamma}, '+lepton+' 2012, gg'+mass)
          signal_52X.Draw('hist')
          signal_52X.SetTitleOffset(1.2,"Y")
          signal_52X.GetYaxis().CenterTitle()
          signal_52X.GetYaxis().SetTitle('Entries (Normalized)')

          gPad.Update()
          stats52 = signal_52X.GetListOfFunctions().FindObject("stats")
          stats52.SetLineColor(kBlue)
          stats52.SetTextColor(kBlue)
          stats52.SetLineWidth(2)
          stats52.SetY1NDC(stats52.GetY1NDC()-0.1)
          stats52.SetY2NDC(stats52.GetY2NDC()-0.1)
          stats52.SetX1NDC(stats52.GetX1NDC()-0.1)
          stats52.SetX2NDC(stats52.GetX2NDC()-0.1)
          stats52.AddText('52X')
          stats52.Draw()

          signal_53X.Draw('histsames')
          gPad.Update()

          stats53 = signal_53X.GetListOfFunctions().FindObject("stats")
          stats53.SetLineColor(kRed)
          stats53.SetTextColor(kRed)
          stats53.SetLineWidth(2)
          stats53.SetY1NDC(stats53.GetY1NDC()-0.27)
          stats53.SetY2NDC(stats53.GetY2NDC()-0.27)
          stats53.SetX1NDC(stats53.GetX1NDC()-0.1)
          stats53.SetX2NDC(stats53.GetX2NDC()-0.1)
          stats52.AddText('53X')
          stats53.Draw()
          canvas.SaveAs(lepton+'_52Xvs53X_'+mass+'.pdf')

        signalCat_52X = fileDict52X[lepton].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal2012gg'+mass)
        signalCat_53X = fileDict53X[lepton].Get('h1_InvariantMassReco1Gev'+cat+'FULLRANGE_Signal2012gg'+mass)
        rms_52X.append(signalCat_52X.GetRMS())
        rms_53X.append(signalCat_53X.GetRMS())

      canvas.cd()

      rmsGraph_52X = TGraph(7,massArray,np.asarray(rms_52X))
      rmsGraph_53X = TGraph(7,massArray,np.asarray(rms_53X))

      rmsGraph_52X.SetTitle(lepton+' '+cat+' #sigma per Mass')
      rmsGraph_52X.GetXaxis().SetTitle('m_{H} (GeV)')
      rmsGraph_52X.GetYaxis().SetTitle('#sigma (GeV)')
      rmsGraph_52X.GetYaxis().CenterTitle()
      rmsGraph_52X.SetMinimum(0.9*min(rms_52X+rms_53X))
      rmsGraph_52X.SetMaximum(1.1*max(rms_52X+rms_53X))

      rmsGraph_52X.SetLineColor(kBlue)
      rmsGraph_53X.SetLineColor(kRed)

      rmsGraph_52X.Draw('ALP')
      rmsGraph_53X.Draw('LPsame')

      legend = TLegend(0.7,0.2,0.85,0.4,'',"brNDC")
      legend.SetBorderSize(1)
      legend.SetTextSize(0.04)
      legend.SetFillColor(0)
      legend.AddEntry(rmsGraph_52X, '52X', 'l')
      legend.AddEntry(rmsGraph_53X, '53X', 'l')
      legend.Draw()

      canvas.SaveAs(lepton+' '+cat+'_52Xvs53X_sigma.pdf')



def AccTables():
  mu2012File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_NewPhoCutV2_3-29-13.root")
  el2012File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_NewPhoCutV2_3-29-13.root")
  FileDict = {'mu2012':mu2012File,'el2012':el2012File}
  for selection in FileDict.keys():
    accHisto = FileDict[selection].GetDirectory('Misc').Get('h1_acceptanceByCut_Signal2012ggM125')
    initialEvents = accHisto.Integral(1,1)/3
    accXeffList =[]
    accXeffList.append(accHisto.Integral(50,50)/initialEvents)
    accXeffList.append(accHisto.Integral(60,60)/initialEvents)
    accXeffList.append(accHisto.Integral(61,61)/initialEvents)
    accXeffList.append(accHisto.Integral(62,62)/initialEvents)
    accXeffList.append(accHisto.Integral(63,63)/initialEvents)

    print selection
    print 'total', 'cat1', 'cat2', 'cat3', 'cat4'
    print(' '.join('{0:.2f}'.format(k) for k in accXeffList))



def AccPlots():

  mu2012FileWithR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_Nominal_12-5-12.root")
  mu2011FileWithR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_Nominal_12-8-12.root")
  el2012FileWithR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012_Medium_Nominal_12-5-12.root")
  el2011FileWithR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2011_Medium_Nominal_12-5-12.root")

  mu2012FileWithOutR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_noR9_12-5-12.root")
  mu2011FileWithOutR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_noR9_12-5-12.root")
  el2012FileWithOutR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012_Medium_noR9_12-5-12.root")
  el2011FileWithOutR9 = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2011_Medium_noR9_12-5-12.root")

  mu2012FileTest = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/localHistos/higgsHistograms_ggM125_SYNC_V11.root")
  mu2011FileTest = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/localHistos/higgsHistograms_ggM125_7TeV.root")

  FileDictWithR9 = {'mu2012':mu2012FileWithR9,'el2012':el2012FileWithR9,'mu2011':mu2011FileWithR9,'el2011':el2011FileWithR9}
  FileDictWithOutR9 = {'mu2012':mu2012FileWithOutR9,'el2012':el2012FileWithOutR9,'mu2011':mu2011FileWithOutR9,'el2011':el2011FileWithOutR9}
  massList = ['M120','M125','M130','M135','M140']
  massArray = np.array([120.0,125.0,130.0,135.0,140.0])


  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(3)
  gStyle.SetTitleH(0.09) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(1.01)    # Set the position of the title box
  gStyle.SetLineWidth(2)

  canvas = TCanvas('canvas','canvas',800,600)
  canvas.SetGrid()
  for selection in FileDictWithR9.keys():

    eventByCutWithR9 = []
    accByMassCat1WithR9 = []
    accByMassCat2WithR9 = []
    accByMassCat3WithR9 = []
    accByMassCat4WithR9 = []
    accByMassWithR9 = []
    eventByCutWithOutR9 = []
    accByMassCat1WithOutR9 = []
    accByMassCat4WithOutR9 = []
    for i,mass in enumerate(massList):
      if '2012' in selection:
        eventByCutWithR9.append(FileDictWithR9[selection].GetDirectory('Misc').Get('h1_acceptanceByCut_Signal2012gg'+mass))
        eventByCutWithOutR9.append(FileDictWithOutR9[selection].GetDirectory('Misc').Get('h1_acceptanceByCut_Signal2012gg'+mass))
      else:
        eventByCutWithR9.append(FileDictWithR9[selection].GetDirectory('Misc').Get('h1_acceptanceByCut_Signal2011gg'+mass))
        eventByCutWithOutR9.append(FileDictWithOutR9[selection].GetDirectory('Misc').Get('h1_acceptanceByCut_Signal2011gg'+mass))
      accByMassCat1WithR9.append(eventByCutWithR9[i].Integral(60,60)/eventByCutWithR9[i].Integral(1,1)*300)
      accByMassCat2WithR9.append(eventByCutWithR9[i].Integral(61,61)/eventByCutWithR9[i].Integral(1,1)*300)
      accByMassCat3WithR9.append(eventByCutWithR9[i].Integral(62,62)/eventByCutWithR9[i].Integral(1,1)*300)
      accByMassCat4WithR9.append(eventByCutWithR9[i].Integral(63,63)/eventByCutWithR9[i].Integral(1,1)*300)

      accByMassCat1WithOutR9.append(eventByCutWithOutR9[i].Integral(60,60)/eventByCutWithOutR9[i].Integral(1,1)*300)
      accByMassCat4WithOutR9.append(eventByCutWithOutR9[i].Integral(63,63)/eventByCutWithOutR9[i].Integral(1,1)*300)

      accByMassWithR9.append(eventByCutWithR9[i].Integral(50,50)/eventByCutWithR9[i].Integral(1,1)*300)

    cat1GraphWithR9 = TGraph(5,massArray,np.asarray(accByMassCat1WithR9))
    cat2GraphWithR9 = TGraph(5,massArray,np.asarray(accByMassCat2WithR9))
    cat3GraphWithR9 = TGraph(5,massArray,np.asarray(accByMassCat3WithR9))
    cat4GraphWithR9 = TGraph(5,massArray,np.asarray(accByMassCat4WithR9))

    cat1GraphWithOutR9 = TGraph(5,massArray,np.asarray(accByMassCat1WithOutR9))
    cat4GraphWithOutR9 = TGraph(5,massArray,np.asarray(accByMassCat4WithOutR9))

    GraphWithR9 = TGraph(5,massArray,np.asarray(accByMassWithR9))

    cat1GraphWithR9.SetMaximum(30)
    cat1GraphWithR9.SetMinimum(0)
    cat1GraphWithR9.SetTitle(selection+' Acc*Eff for R9 corrections')
    cat1GraphWithR9.GetXaxis().SetTitle('m_{H} (GeV)')
    cat1GraphWithR9.GetYaxis().SetTitle('Acc*Eff (%)')
    cat1GraphWithR9.GetYaxis().CenterTitle()
    cat1GraphWithR9.SetLineColor(kBlue)
    cat1GraphWithR9.SetMarkerColor(kBlue)
    cat1GraphWithR9.SetMarkerStyle(20)
    cat1GraphWithR9.Draw('ALP')

    cat1GraphWithOutR9.SetLineStyle(2)
    cat1GraphWithOutR9.SetLineColor(kBlue)
    cat1GraphWithOutR9.SetMarkerColor(kBlue)
    cat1GraphWithOutR9.SetMarkerStyle(21)
    cat1GraphWithOutR9.Draw('LPsame')

    cat4GraphWithR9.SetLineColor(kGreen)
    cat4GraphWithR9.SetMarkerColor(kGreen)
    cat4GraphWithR9.SetMarkerStyle(20)
    cat4GraphWithR9.Draw('LPsame')

    cat4GraphWithOutR9.SetLineColor(kGreen)
    cat4GraphWithOutR9.SetLineStyle(2)
    cat4GraphWithOutR9.SetMarkerColor(kGreen)
    cat4GraphWithOutR9.SetMarkerStyle(21)
    cat4GraphWithOutR9.Draw('LPsame')

    '''
    cat2GraphWithR9.SetLineColor(kGreen)
    cat2GraphWithR9.SetMarkerColor(kGreen)
    cat2GraphWithR9.SetMarkerStyle(20)
    cat2GraphWithR9.Draw('LPsame')

    cat3GraphWithR9.SetLineColor(kBlack)
    cat3GraphWithR9.SetMarkerColor(kBlack)
    cat3GraphWithR9.SetMarkerStyle(20)
    cat3GraphWithR9.Draw('LPsame')

    GraphWithR9.SetLineColor(kGreen)
    GraphWithR9.SetMarkerColor(kGreen)
    GraphWithR9.SetMarkerStyle(20)
    GraphWithR9.Draw('LPsame')
    '''

    legend = TLegend(0.30,0.55,0.63,0.89,'',"brNDC")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.05)
    legend.SetFillColor(0)
    legend.AddEntry(cat1GraphWithR9, 'CAT1 with R9 correction','pl')
    legend.AddEntry(cat1GraphWithOutR9, 'CAT1 without R9 correction','pl')
    legend.AddEntry(cat4GraphWithR9, 'CAT4 with R9 correction','pl')
    legend.AddEntry(cat4GraphWithOutR9, 'CAT4 without R9 correction','pl')
    legend.Draw()

    canvas.SaveAs(selection+'accR9.pdf')

    if '2011' in selection:
      photonPtPreHistWithR9 = FileDictWithR9[selection].GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2011ggM125')
      photonPtPostHistWithR9 = FileDictWithR9[selection].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2011ggM125')
      photonPtPreHistWithOutR9 = FileDictWithOutR9[selection].GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2011ggM125')
      photonPtPostHistWithOutR9 = FileDictWithOutR9[selection].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2011ggM125')
    else:
      photonPtPreHistWithR9 = FileDictWithR9[selection].GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2012ggM125')
      photonPtPostHistWithR9 = FileDictWithR9[selection].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2012ggM125')
      photonPtPreHistWithOutR9 = FileDictWithOutR9[selection].GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2012ggM125')
      photonPtPostHistWithOutR9 = FileDictWithOutR9[selection].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2012ggM125')

    photonPtPostHistWithR9.Divide(photonPtPostHistWithR9, photonPtPreHistWithR9,1,1, "B")
    photonPtPostHistWithR9.SetFillColor(kBlue)
    photonPtPostHistWithR9.SetFillStyle(3001)
    photonPtPostHistWithR9.SetLineColor(kBlack)
    photonPtPostHistWithR9.SetLineWidth(2)
    photonPtPostHistWithR9.SetTitle(selection+' Acc*Eff for Photon pT')
    photonPtPostHistWithR9.GetXaxis().SetTitle('pT (GeV)')
    photonPtPostHistWithR9.GetYaxis().SetTitle('Acc*Eff')
    photonPtPostHistWithR9.GetYaxis().CenterTitle()
    photonPtPostHistWithR9.Draw('e5')
    photonPtPostHistWithR9.SetMarkerStyle(8)
    gStyle.SetErrorX(0)
    photonPtPostHistWithR9.Draw('pX0same')

    photonPtPostHistWithOutR9.Divide(photonPtPostHistWithOutR9, photonPtPreHistWithOutR9,1,1, "B")
    photonPtPostHistWithOutR9.SetFillColor(kGreen)
    photonPtPostHistWithOutR9.SetFillStyle(3001)
    photonPtPostHistWithOutR9.SetLineColor(kBlack)
    photonPtPostHistWithOutR9.SetLineWidth(2)
    photonPtPostHistWithOutR9.Draw('samee5')
    photonPtPostHistWithOutR9.SetMarkerStyle(8)
    gStyle.SetErrorX(0)
    photonPtPostHistWithOutR9.Draw('pX0same')

    legend = TLegend(0.18,0.2,0.3,0.4,'',"brNDC")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.SetFillColor(0)
    legend.AddEntry(photonPtPostHistWithR9, 'Photon Acc*Eff, w/ R9 corr', 'F')
    legend.AddEntry(photonPtPostHistWithOutR9, 'Photon Acc*Eff, w/o R9 corr', 'F')
    legend.Draw()

    canvas.SaveAs('photonPtAcc_'+selection+'.pdf')

  photonPtPreHistTest2012 = mu2012FileTest.GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2012ggM125')
  photonPtPostHistTest2012 = mu2012FileTest.GetDirectory('PostGen').Get('h1_photonPtPostGen_Signal2012ggM125')
  photonPtPreHistTest2011 = mu2011FileTest.GetDirectory('PreGen').Get('h1_photonPtPreGen_Signal2011ggM125')
  photonPtPostHistTest2011 = mu2011FileTest.GetDirectory('PostGen').Get('h1_photonPtPostGen_Signal2011ggM125')

  photonPtPostHistTest2012.Divide(photonPtPostHistTest2012, photonPtPreHistTest2012,1,1, "B")
  photonPtPostHistTest2011.Divide(photonPtPostHistTest2011, photonPtPreHistTest2011,1,1, "B")
  photonPtPostHistTest2012.SetFillColor(kBlue)
  photonPtPostHistTest2012.SetFillStyle(3001)
  photonPtPostHistTest2012.SetLineColor(kBlack)
  photonPtPostHistTest2012.SetLineWidth(2)
  photonPtPostHistTest2011.SetFillColor(kGreen)
  photonPtPostHistTest2011.SetFillStyle(3001)
  photonPtPostHistTest2011.SetLineColor(kBlack)
  photonPtPostHistTest2011.SetLineWidth(2)
  photonPtPostHistTest2012.SetTitle('mu2011 vs mu2012 Acc for Photon pT')
  photonPtPostHistTest2012.GetXaxis().SetTitle('pT (GeV)')
  photonPtPostHistTest2012.GetYaxis().SetTitle('Acc')
  photonPtPostHistTest2012.GetYaxis().CenterTitle()
  photonPtPostHistTest2012.Draw('e5')
  photonPtPostHistTest2012.SetMarkerStyle(8)
  gStyle.SetErrorX(0)
  photonPtPostHistTest2012.Draw('pX0same')
  photonPtPostHistTest2011.Draw('e5same')
  photonPtPostHistTest2011.SetMarkerStyle(8)
  photonPtPostHistTest2011.Draw('pX0same')

  legend = TLegend(0.18,0.2,0.3,0.4,'',"brNDC")
  legend.SetBorderSize(0)
  legend.SetTextSize(0.04)
  legend.SetFillColor(0)
  legend.AddEntry(photonPtPostHistTest2012, 'Photon Acc 2012', 'F')
  legend.AddEntry(photonPtPostHistTest2011, 'Photon Acc 2011', 'F')
  legend.Draw()

  canvas.SaveAs('photonPtAcc_TEST.pdf')

  relChange = np.array([-4.4,-4.9,-5.4,-5.2,-5.4])
  relChangeGraph = TGraph(5,massArray,relChange)
  relChangeGraph.SetTitle('Relative Change in Limit Sensitivity Due to R9 Migration')
  relChangeGraph.GetXaxis().SetTitle('m_{H} GeV')
  relChangeGraph.GetYaxis().SetTitle('Relative Change in Limit (%)')
  relChangeGraph.GetYaxis().SetNoExponent()
  relChangeGraph.GetYaxis().CenterTitle()
  relChangeGraph.Draw('ALP')

  canvas.SaveAs('relChangeR9.pdf')








def BGPlots():
  mu2012File = TFile("mu2012_higgsOutput_Medium_11-30-12.root")
  mu2011File = TFile("mu2011_higgsOutput_Medium_11-30-12.root")
  el2012File = TFile("el2012_higgsOutput_Medium_11-30-12.root")
  el2011File = TFile("el2011_higgsOutput_Medium_11-30-12.root")

  FileList = {'mu2012':mu2012File,'mu2011':mu2011File,'el2012':el2012File,'el2011':el2011File}

  gStyle.SetOptStat(0)
  gStyle.SetOptFit(0)
  gStyle.SetTitleFontSize(3)
  gStyle.SetTitleH(0.09) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(1)    # Set the position of the title box

  canvas = TCanvas('canvas','canvas',800,600)
  for selection in FileList.keys():

    blindData = FileList[selection].Get('data_obs')
    if '2012' in selection:
      blindData.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 5.2 fb^{-1}, H#rightarrowZ#gamma")
    else:
      blindData.SetTitle("CMS Preliminary #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}, H#rightarrowZ#gamma")
    fit = TF1("expo_data_fit","expo")
    fit.SetLineColor(kBlack)
    fit.SetLineWidth(2)
    blindData.Fit(fit)
    DYfit = TF1("expo_DY_fit","expo",115,180)
    if selection == 'mu2012':
      DYfit.SetParameters(8.8470, -0.045299)
    elif selection == 'mu2011':
      DYfit.SetParameters(7.6875, -0.038857)
    elif selection == 'el2012':
      DYfit.SetParameters(8.6192, -0.044034)
    elif selection == 'el2011':
      DYfit.SetParameters(7.33350, -0.036891)
    DYfit.SetLineWidth(2)
    DYfit.SetLineStyle(2)
    DYfit.SetLineColor(kBlue)
    blindData.GetXaxis().SetTitleSize(0.055)
    blindData.GetXaxis().SetLabelSize(0.05)
    blindData.GetYaxis().SetTitleSize(0.055)
    blindData.GetYaxis().SetLabelSize(0.05)
    blindData.GetYaxis().CenterTitle()
    if 'el' in selection:
      blindData.GetXaxis().SetTitle("m_{ee#gamma} (GeV)")
    else:
      blindData.GetXaxis().SetTitle("m_{#mu#mu#gamma} (GeV)")
    blindData.GetYaxis().SetTitle("Events per GeV")
    blindData.Draw('ep')
    fit.Draw('same')
    DYfit.Draw('same')


    SignalM125 = FileList[selection].Get(selection+'_sigMC_M125')
    SignalM125.Scale(100)
    SignalM125.SetLineWidth(2)
    SignalM125.SetLineColor(kRed)
    SignalM125.SetLineStyle(2)
    SignalM125.Draw('histsame')

    legend = TLegend(0.44,0.55,0.63,0.89,'',"brNDC")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.06)
    legend.SetFillColor(0)
    legend.AddEntry(blindData, 'Data','pe')
    legend.AddEntry(fit, 'Fit to Data', 'l')
    legend.AddEntry(DYfit, 'Drell-Yan Estimate', 'l')
    legend.AddEntry(SignalM125, 'm_{H} = 125 GeV (x100)','l')
    legend.Draw()

    txtBox = TPaveText(0.22,0.85,0.42,0.90,"brNDC")
    txtBox.SetBorderSize(0)
    txtBox.SetTextSize(0.055)
    txtBox.SetFillColor(0)

    if selection == 'mu2012':
      txtBox.AddText("2012 Muon")
    elif selection == 'mu2011':
      txtBox.AddText("2011 Muon")
    elif selection == 'el2012':
      txtBox.AddText("2012 Electron")
    elif selection == 'el2011':
      txtBox.AddText("2011 Electron")
    txtBox.Draw()


    canvas.SaveAs(selection+'_bg_excl_PAS.pdf')
    raw_input('Hit any key to continue')
    legend.Clear()
    txtBox.Clear()


def RochCorPlots():
  mu2012FileCor = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_MuonCor_11-30-12.root")
  mu2012FileNoCor = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012_Medium_MuonNoCor_11-30-12.root")
  mu2011FileCor = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_MuonCor_11-30-12.root")
  mu2011FileNoCor = TFile("/uscms/home/bpollack/work/CMSSW_5_3_2/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_MuonNoCor_11-30-12.root")


  FileDict = {'mu2012cor':mu2012FileCor,'mu2012no':mu2012FileNoCor,'mu2011cor':mu2011FileCor,'mu2011no':mu2011FileNoCor}
  FileList2011 = [mu2011FileCor,mu2011FileNoCor]
  FileList2012 = [mu2012FileCor,mu2012FileNoCor]

  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(3)
  gStyle.SetTitleH(0.09) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(1)    # Set the position of the title box

  for muonkey in FileDict.keys():
    diLepMassDATA = FileDict[muonkey].GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_DATA')
    if '2012' in muonkey:
      diLepMassSignal = FileDict[muonkey].GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_Signal2012ggM125')
    else:
      diLepMassSignal = FileDict[muonkey].GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_Signal2011ggM125')

    diLepCan = TCanvas('diLepCan','canvas',800,600)
    diLepCan.cd()
    diLepMassDATA.SetLineColor(kBlack)
    diLepMassDATA.SetLineWidth(2)
    diLepMassSignal.SetLineColor(kRed)
    diLepMassSignal.SetLineWidth(2)
    diLepMassSignal.GetYaxis().SetTitleSize(0.055)
    diLepMassSignal.GetYaxis().CenterTitle()
    diLepMassSignal.GetXaxis().SetTitleSize(0.055)
    diLepMassSignal.GetYaxis().SetLabelSize(0.05)
    diLepMassSignal.GetXaxis().SetLabelSize(0.05)
    diLepMassSignal.GetYaxis().SetTitle("Arbitrary Units")
    if 'cor' in muonkey:
      diLepMassSignal.GetXaxis().SetTitle("m_{#mu#mu} (GeV), corrected")
    else:
      diLepMassSignal.GetXaxis().SetTitle("m_{#mu#mu} (GeV), uncorrected")
    if '2012' in muonkey:
      diLepMassSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 5.2 fb^{-1}")
    else:
      diLepMassSignal.SetTitle("CMS Preliminary #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}")
    diLepMassSignal.SetTitleOffset(1.2,"Y")
    diLepMassSignal.DrawNormalized('hist')
    diLepMassDATA.DrawNormalized('histsame')
    diLepLeg = TLegend(0.56,0.65,0.75,0.89,'',"brNDC")
    diLepLeg.SetBorderSize(0)
    diLepLeg.SetTextSize(0.06)
    diLepLeg.SetFillColor(0)
    diLepLeg.AddEntry(diLepMassSignal,'m_{H} = 125 GeV', 'l')
    diLepLeg.AddEntry(diLepMassDATA, 'Data', 'l')
    diLepLeg.Draw()
    diLepCan.SaveAs('diLeptonM_'+muonkey+'_PAS.pdf')

  for i,muonFile in enumerate(FileList2011):
    diLepMassDATA = muonFile.GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_DATA')
    diLepMassSignal =  muonFile.GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_Signal2011ggM125')

    diLepMassDATA.SetLineColor(kBlack)
    diLepMassDATA.SetLineWidth(2)
    diLepMassSignal.SetLineColor(kRed)
    diLepMassSignal.SetLineWidth(2)

    if i == 0:
      diLepCanDATA = TCanvas('diLepCanDATA','canvas',800,600)
      diLepCanDATA.cd()
      diLepMassDATA.GetYaxis().SetTitleSize(0.055)
      diLepMassDATA.GetYaxis().CenterTitle()
      diLepMassDATA.GetXaxis().SetTitleSize(0.055)
      diLepMassDATA.GetYaxis().SetLabelSize(0.05)
      diLepMassDATA.GetXaxis().SetLabelSize(0.05)
      diLepMassDATA.GetYaxis().SetTitle("Arbitrary Units")
      diLepMassDATA.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
      diLepMassDATA.SetTitle("CMS Preliminary #sqrt{s} = 7 TeV, L = 5.0 fb^{-1}")
      diLepMassDATA.SetTitleOffset(1.2,"Y")
      diLepMassDATA.DrawNormalized('hist')
      diLepLegDATA = TLegend(0.56,0.65,0.75,0.89,'',"brNDC")
      diLepLegDATA.SetBorderSize(0)
      diLepLegDATA.SetTextSize(0.05)
      diLepLegDATA.SetFillColor(0)
      diLepLegDATA.AddEntry(diLepMassDATA, 'Data, corrected', 'l')
    if i == 1:
      diLepCanDATA.cd()
      diLepMassDATA.SetLineStyle(2)
      diLepMassDATA.DrawNormalized('histsame')
      diLepLegDATA.AddEntry(diLepMassDATA, 'Data, uncorrected', 'l')
      diLepLegDATA.Draw()

    if i == 0:
      diLepCanSignal = TCanvas('diLepCanSignal','canvas',800,600)
      diLepCanSignal.cd()
      diLepMassSignal.GetYaxis().SetTitleSize(0.055)
      diLepMassSignal.GetYaxis().CenterTitle()
      diLepMassSignal.GetXaxis().SetTitleSize(0.055)
      diLepMassSignal.GetYaxis().SetLabelSize(0.05)
      diLepMassSignal.GetXaxis().SetLabelSize(0.05)
      diLepMassSignal.GetYaxis().SetTitle("Arbitrary Units")
      diLepMassSignal.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
      diLepMassSignal.SetTitle("2011 Signal, m_{H} = 125 (GeV)")
      diLepMassSignal.SetTitleOffset(1.2,"Y")
      diLepMassSignal.SetMaximum(1670)
      diLepMassSignal.DrawNormalized('hist')
      diLepLegSignal = TLegend(0.56,0.65,0.75,0.89,'',"brNDC")
      diLepLegSignal.SetBorderSize(0)
      diLepLegSignal.SetTextSize(0.05)
      diLepLegSignal.SetFillColor(0)
      diLepLegSignal.AddEntry(diLepMassSignal, 'MC, corrected', 'l')
    if i == 1:
      diLepCanSignal.cd()
      diLepMassSignal.SetLineStyle(2)
      diLepMassSignal.DrawNormalized('histsameA')
      diLepLegSignal.AddEntry(diLepMassSignal, 'MC, uncorrected', 'l')
      diLepLegSignal.Draw()

  diLepCanDATA.SaveAs('diLeptonM_2011_DATA_PAS.pdf')
  diLepCanSignal.SaveAs('diLeptonM_2011_Signal_PAS.pdf')

  for i,muonFile in enumerate(FileList2012):
    diLepMassDATA = muonFile.GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_DATA')
    diLepMassSignal =  muonFile.GetDirectory('pT-Eta-Phi').Get('h1_diLepMass_Signal2012ggM125')

    diLepMassDATA.SetLineColor(kBlack)
    diLepMassDATA.SetLineWidth(2)
    diLepMassSignal.SetLineColor(kRed)
    diLepMassSignal.SetLineWidth(2)

    if i == 0:
      diLepCanDATA = TCanvas('diLepCanDATA','canvas',800,600)
      diLepCanDATA.cd()
      diLepMassDATA.GetYaxis().SetTitleSize(0.055)
      diLepMassDATA.GetYaxis().CenterTitle()
      diLepMassDATA.GetXaxis().SetTitleSize(0.055)
      diLepMassDATA.GetYaxis().SetLabelSize(0.05)
      diLepMassDATA.GetXaxis().SetLabelSize(0.05)
      diLepMassDATA.GetYaxis().SetTitle("Arbitrary Units")
      diLepMassDATA.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
      diLepMassDATA.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 5.2 fb^{-1}")
      diLepMassDATA.SetTitleOffset(1.2,"Y")
      diLepMassDATA.DrawNormalized('hist')
      diLepLegDATA = TLegend(0.56,0.65,0.75,0.89,'',"brNDC")
      diLepLegDATA.SetBorderSize(0)
      diLepLegDATA.SetTextSize(0.05)
      diLepLegDATA.SetFillColor(0)
      diLepLegDATA.AddEntry(diLepMassDATA, 'Data, corrected', 'l')
    if i == 1:
      diLepCanDATA.cd()
      diLepMassDATA.SetLineStyle(2)
      diLepMassDATA.DrawNormalized('histsame')
      diLepLegDATA.AddEntry(diLepMassDATA, 'Data, uncorrected', 'l')
      diLepLegDATA.Draw()

    if i == 0:
      diLepCanSignal = TCanvas('diLepCanSignal','canvas',800,600)
      diLepCanSignal.cd()
      diLepMassSignal.GetYaxis().SetTitleSize(0.055)
      diLepMassSignal.GetYaxis().CenterTitle()
      diLepMassSignal.GetXaxis().SetTitleSize(0.055)
      diLepMassSignal.GetYaxis().SetLabelSize(0.05)
      diLepMassSignal.GetXaxis().SetLabelSize(0.05)
      diLepMassSignal.GetYaxis().SetTitle("Arbitrary Units")
      diLepMassSignal.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
      diLepMassSignal.SetTitle("2012 Signal, m_{H} = 125 (GeV)")
      diLepMassSignal.SetTitleOffset(1.2,"Y")
      diLepMassSignal.SetMaximum(1500)
      diLepMassSignal.DrawNormalized('hist')
      diLepLegSignal = TLegend(0.56,0.65,0.75,0.89,'',"brNDC")
      diLepLegSignal.SetBorderSize(0)
      diLepLegSignal.SetTextSize(0.05)
      diLepLegSignal.SetFillColor(0)
      diLepLegSignal.AddEntry(diLepMassSignal, 'MC, corrected', 'l')
    if i == 1:
      diLepCanSignal.cd()
      diLepMassSignal.SetLineStyle(2)
      diLepMassSignal.DrawNormalized('histsameA')
      diLepLegSignal.AddEntry(diLepMassSignal, 'MC, uncorrected', 'l')
      diLepLegSignal.Draw()

  diLepCanDATA.SaveAs('diLeptonM_2012_DATA_PAS.pdf')
  diLepCanSignal.SaveAs('diLeptonM_2012_Signal_PAS.pdf')



def PTPlots():
  mu2012File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_Medium_MoriondV3_1-28-13.root")
  #mu2011File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2011_Medium_DataSyncAllCor_11-30-12.root")
  el2012File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_Medium_MoriondV3_1-28-13.root")
  #el2011File = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2011_Medium_MoriondV1_1-28-13.root")


  FileList = [mu2012File,el2012File]
  #FileList = [mu2012File,el2012File,mu2011File,el2011File]

  photonPtDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_DATA')
  photonPtSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2012ggM125')
  photonEtaDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonEta_DATA')
  photonEtaSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonEta_Signal2012ggM125')
  photonPhiDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonPhi_DATA')
  photonPhiSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_photonPhi_Signal2012ggM125')

  leadingPtDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPt_DATA')
  leadingPtSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPt_Signal2012ggM125')
  leadingEtaDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdEta_DATA')
  leadingEtaSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdEta_Signal2012ggM125')
  leadingPhiDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPhi_DATA')
  leadingPhiSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPhi_Signal2012ggM125')

  trailingPtDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPt_DATA')
  trailingPtSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPt_Signal2012ggM125')
  trailingEtaDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdEta_DATA')
  trailingEtaSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdEta_Signal2012ggM125')
  trailingPhiDATA = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPhi_DATA')
  trailingPhiSignal = FileList[0].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPhi_Signal2012ggM125')

  r9DATA = FileList[0].Get('h1_R9CorFull_DATA')
  r9Signal = FileList[0].Get('h1_R9CorFull_Signal2012ggM125')

  gStyle.SetOptStat(0)
  gStyle.SetTitleFontSize(3)
  gStyle.SetTitleH(0.09) # Set the height of the title box
  gStyle.SetTitleW(1)    # Set the width of the title box
  gStyle.SetTitleX(0)    # Set the position of the title box
  gStyle.SetTitleY(1)    # Set the position of the title box

  for i in range(1,2):
    photonPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_DATA'))
    leadingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPt_DATA'))
    trailingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPt_DATA'))
    photonEtaDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonEta_DATA'))
    leadingEtaDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdEta_DATA'))
    trailingEtaDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdEta_DATA'))
    photonPhiDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonPhi_DATA'))
    leadingPhiDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPhi_DATA'))
    trailingPhiDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPhi_DATA'))
    r9DATA.Add(FileList[i].Get('h1_R9CorFull_DATA'))
    if i < 2:
      photonPtSignal.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2012ggM125'))
      leadingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPt_Signal2012ggM125'))
      trailingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPt_Signal2012ggM125'))
      photonEtaSignal.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonEta_Signal2012ggM125'))
      leadingEtaDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdEta_Signal2012ggM125'))
      trailingEtaDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdEta_Signal2012ggM125'))
      photonPhiSignal.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonPhi_Signal2012ggM125'))
      leadingPhiDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPhi_Signal2012ggM125'))
      trailingPhiDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPhi_Signal2012ggM125'))
      r9DATA.Add(FileList[i].Get('h1_R9CorFull_Signal2012ggM125'))
    else:
      photonPtSignal.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_photonPt_Signal2011ggM125'))
      leadingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_leadLeptonStdPt_Signal2011ggM125'))
      trailingPtDATA.Add(FileList[i].GetDirectory('pT-Eta-Phi').Get('h1_trailingLeptonStdPt_Signal2011ggM125'))
      r9DATA.Add(FileList[i].Get('h1_R9CorFull_Signal2011ggM125'))

  phoCan = TCanvas('phoCan','canvas',800,600)
  phoCan.cd()
  photonPtDATA.SetLineColor(kBlack)
  photonPtDATA.SetLineWidth(2)
  photonPtDATA.SetMarkerStyle(8)
  photonPtSignal.SetLineColor(kRed)
  photonPtSignal.SetLineWidth(2)
  photonPtDATA.GetYaxis().SetTitleSize(0.055)
  photonPtDATA.GetYaxis().CenterTitle()
  photonPtDATA.GetXaxis().SetTitleSize(0.055)
  photonPtDATA.GetYaxis().SetLabelSize(0.05)
  photonPtDATA.GetXaxis().SetLabelSize(0.05)
  photonPtDATA.GetYaxis().SetTitle("Arbitrary Units")
  photonPtDATA.GetXaxis().SetTitle("#gamma p_{T} (GeV)")
  photonPtDATA.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  photonPtDATA.SetTitleOffset(1.2,"Y")
  photonPtDATA.DrawNormalized('ep')
  photonPtSignal.DrawNormalized('histsame')
  phoLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  phoLeg.SetBorderSize(0)
  phoLeg.SetTextSize(0.06)
  phoLeg.SetFillColor(0)
  phoLeg.AddEntry(photonPtSignal,'m_{H} = 125 GeV', 'l')
  phoLeg.AddEntry(photonPtDATA, 'Data', 'ep')
  phoLeg.Draw()
  phoCan.SaveAs('photonPT_PAS.pdf')

  photonEtaDATA.SetLineColor(kBlack)
  photonEtaDATA.SetLineWidth(2)
  photonEtaDATA.SetMarkerStyle(8)
  photonEtaSignal.SetLineColor(kRed)
  photonEtaSignal.SetLineWidth(2)
  photonEtaDATA.GetYaxis().SetTitleSize(0.055)
  photonEtaDATA.GetYaxis().CenterTitle()
  photonEtaDATA.GetXaxis().SetTitleSize(0.055)
  photonEtaDATA.GetYaxis().SetLabelSize(0.05)
  photonEtaDATA.GetXaxis().SetLabelSize(0.05)
  photonEtaDATA.GetYaxis().SetTitle("Arbitrary Units")
  photonEtaDATA.GetXaxis().SetTitle("#gamma #eta")
  photonEtaDATA.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  photonEtaDATA.SetTitleOffset(1.2,"Y")
  photonEtaSignal.DrawNormalized('hist')
  photonEtaDATA.DrawNormalized('epsame')
  phoLeg = TLegend(0.34,0.25,0.55,0.49,'',"brNDC")
  phoLeg.SetBorderSize(0)
  phoLeg.SetTextSize(0.06)
  phoLeg.SetFillColor(0)
  phoLeg.AddEntry(photonEtaSignal,'m_{H} = 125 GeV', 'l')
  phoLeg.AddEntry(photonEtaDATA, 'Data', 'ep')
  phoLeg.Draw()
  phoCan.SaveAs('photonEta_PAS.pdf')

  photonPhiDATA.SetLineColor(kBlack)
  photonPhiDATA.SetLineWidth(2)
  photonPhiDATA.SetMarkerStyle(8)
  photonPhiSignal.SetLineColor(kRed)
  photonPhiSignal.SetLineWidth(2)
  photonPhiDATA.GetYaxis().SetTitleSize(0.055)
  photonPhiDATA.GetYaxis().CenterTitle()
  photonPhiDATA.GetXaxis().SetTitleSize(0.055)
  photonPhiDATA.GetYaxis().SetLabelSize(0.05)
  photonPhiDATA.GetXaxis().SetLabelSize(0.05)
  photonPhiDATA.GetYaxis().SetTitle("Arbitrary Units")
  photonPhiDATA.GetXaxis().SetTitle("#gamma #phi")
  photonPhiDATA.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  photonPhiDATA.SetTitleOffset(1.2,"Y")
  photonPhiDATA.DrawNormalized('ep')
  photonPhiSignal.DrawNormalized('histsame')
  phoLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  phoLeg.SetBorderSize(0)
  phoLeg.SetTextSize(0.06)
  phoLeg.SetFillColor(0)
  phoLeg.AddEntry(photonPhiSignal,'m_{H} = 125 GeV', 'l')
  phoLeg.AddEntry(photonPhiDATA, 'Data', 'ep')
  phoLeg.Draw()
  phoCan.SaveAs('photonPhi_PAS.pdf')

  leadingCan = TCanvas('leadingCan','canvas',800,600)
  leadingCan.cd()
  leadingPtDATA.SetLineColor(kBlack)
  leadingPtDATA.SetLineWidth(2)
  leadingPtDATA.SetMarkerStyle(8)
  leadingPtSignal.SetLineColor(kRed)
  leadingPtSignal.SetLineWidth(2)
  leadingPtSignal.GetYaxis().SetTitleSize(0.055)
  leadingPtSignal.GetYaxis().CenterTitle()
  leadingPtSignal.GetXaxis().SetTitleSize(0.055)
  leadingPtSignal.GetYaxis().SetLabelSize(0.05)
  leadingPtSignal.GetXaxis().SetLabelSize(0.05)
  leadingPtSignal.GetYaxis().SetTitle("Arbitrary Units")
  leadingPtSignal.GetXaxis().SetTitle("Leading Lepton p_{T} (GeV)")
  leadingPtSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  leadingPtSignal.SetTitleOffset(1.2,"Y")
  leadingPtSignal.DrawNormalized('hist')
  leadingPtDATA.DrawNormalized('epsame')
  leadingLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  leadingLeg.SetBorderSize(0)
  leadingLeg.SetTextSize(0.06)
  leadingLeg.SetFillColor(0)
  leadingLeg.AddEntry(leadingPtSignal,'m_{H} = 125 GeV', 'l')
  leadingLeg.AddEntry(leadingPtDATA, 'Data', 'ep')
  leadingLeg.Draw()
  leadingCan.SaveAs('leadingPT_PAS.pdf')

  leadingEtaDATA.SetLineColor(kBlack)
  leadingEtaDATA.SetLineWidth(2)
  leadingEtaDATA.SetMarkerStyle(8)
  leadingEtaSignal.SetLineColor(kRed)
  leadingEtaSignal.SetLineWidth(2)
  leadingEtaSignal.GetYaxis().SetTitleSize(0.055)
  leadingEtaSignal.GetYaxis().CenterTitle()
  leadingEtaSignal.GetXaxis().SetTitleSize(0.055)
  leadingEtaSignal.GetYaxis().SetLabelSize(0.05)
  leadingEtaSignal.GetXaxis().SetLabelSize(0.05)
  leadingEtaSignal.GetYaxis().SetTitle("Arbitrary Units")
  leadingEtaSignal.GetXaxis().SetTitle("Leading Lepton #eta")
  leadingEtaSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  leadingEtaSignal.SetTitleOffset(1.2,"Y")
  leadingEtaSignal.DrawNormalized('hist')
  leadingEtaDATA.DrawNormalized('epsame')
  leadingLeg = TLegend(0.34,0.25,0.55,0.49,'',"brNDC")
  leadingLeg.SetBorderSize(0)
  leadingLeg.SetTextSize(0.06)
  leadingLeg.SetFillColor(0)
  leadingLeg.AddEntry(leadingEtaSignal,'m_{H} = 125 GeV', 'l')
  leadingLeg.AddEntry(leadingEtaDATA, 'Data', 'ep')
  leadingLeg.Draw()
  leadingCan.SaveAs('leadingEta_PAS.pdf')

  leadingPhiDATA.SetLineColor(kBlack)
  leadingPhiDATA.SetLineWidth(2)
  leadingPhiDATA.SetMarkerStyle(8)
  leadingPhiSignal.SetLineColor(kRed)
  leadingPhiSignal.SetLineWidth(2)
  leadingPhiSignal.GetYaxis().SetTitleSize(0.055)
  leadingPhiSignal.GetYaxis().CenterTitle()
  leadingPhiSignal.GetXaxis().SetTitleSize(0.055)
  leadingPhiSignal.GetYaxis().SetLabelSize(0.05)
  leadingPhiSignal.GetXaxis().SetLabelSize(0.05)
  leadingPhiSignal.GetYaxis().SetTitle("Arbitrary Units")
  leadingPhiSignal.GetXaxis().SetTitle("Leading Lepton #phi")
  leadingPhiSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  leadingPhiSignal.SetTitleOffset(1.2,"Y")
  leadingPhiSignal.DrawNormalized('hist')
  leadingPhiDATA.DrawNormalized('epsame')
  leadingLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  leadingLeg.SetBorderSize(0)
  leadingLeg.SetTextSize(0.06)
  leadingLeg.SetFillColor(0)
  leadingLeg.AddEntry(leadingPhiSignal,'m_{H} = 125 GeV', 'l')
  leadingLeg.AddEntry(leadingPhiDATA, 'Data', 'ep')
  leadingLeg.Draw()
  leadingCan.SaveAs('leadingPhi_PAS.pdf')

  trailingCan = TCanvas('trailingCan','canvas',800,600)
  trailingCan.cd()
  trailingPtDATA.SetLineColor(kBlack)
  trailingPtDATA.SetLineWidth(2)
  trailingPtDATA.SetMarkerStyle(8)
  trailingPtSignal.SetLineColor(kRed)
  trailingPtSignal.SetLineWidth(2)
  trailingPtSignal.GetYaxis().SetTitleSize(0.055)
  trailingPtSignal.GetYaxis().CenterTitle()
  trailingPtSignal.GetXaxis().SetTitleSize(0.055)
  trailingPtSignal.GetYaxis().SetLabelSize(0.05)
  trailingPtSignal.GetXaxis().SetLabelSize(0.05)
  trailingPtSignal.GetYaxis().SetTitle("Arbitrary Units")
  trailingPtSignal.GetXaxis().SetTitle("Trailing Lepton p_{T} (GeV)")
  trailingPtSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  trailingPtSignal.SetTitleOffset(1.2,"Y")
  trailingPtSignal.DrawNormalized('hist')
  trailingPtDATA.DrawNormalized('epsame')
  trailingLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  trailingLeg.SetBorderSize(0)
  trailingLeg.SetTextSize(0.06)
  trailingLeg.SetFillColor(0)
  trailingLeg.AddEntry(trailingPtSignal,'m_{H} = 125 GeV', 'l')
  trailingLeg.AddEntry(trailingPtDATA, 'Data', 'ep')
  trailingLeg.Draw()
  trailingCan.SaveAs('trailingPT_PAS.pdf')

  trailingEtaDATA.SetLineColor(kBlack)
  trailingEtaDATA.SetLineWidth(2)
  trailingEtaDATA.SetMarkerStyle(8)
  trailingEtaSignal.SetLineColor(kRed)
  trailingEtaSignal.SetLineWidth(2)
  trailingEtaSignal.GetYaxis().SetTitleSize(0.055)
  trailingEtaSignal.GetYaxis().CenterTitle()
  trailingEtaSignal.GetXaxis().SetTitleSize(0.055)
  trailingEtaSignal.GetYaxis().SetLabelSize(0.05)
  trailingEtaSignal.GetXaxis().SetLabelSize(0.05)
  trailingEtaSignal.GetYaxis().SetTitle("Arbitrary Units")
  trailingEtaSignal.GetXaxis().SetTitle("Trailing Lepton #eta")
  trailingEtaSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  trailingEtaSignal.SetTitleOffset(1.2,"Y")
  trailingEtaSignal.DrawNormalized('hist')
  trailingEtaDATA.DrawNormalized('epsame')
  trailingLeg = TLegend(0.34,0.25,0.55,0.49,'',"brNDC")
  trailingLeg.SetBorderSize(0)
  trailingLeg.SetTextSize(0.06)
  trailingLeg.SetFillColor(0)
  trailingLeg.AddEntry(trailingEtaSignal,'m_{H} = 125 GeV', 'l')
  trailingLeg.AddEntry(trailingEtaDATA, 'Data', 'ep')
  trailingLeg.Draw()
  trailingCan.SaveAs('trailingEta_PAS.pdf')

  trailingPhiDATA.SetLineColor(kBlack)
  trailingPhiDATA.SetLineWidth(2)
  trailingPhiDATA.SetMarkerStyle(8)
  trailingPhiSignal.SetLineColor(kRed)
  trailingPhiSignal.SetLineWidth(2)
  trailingPhiSignal.GetYaxis().SetTitleSize(0.055)
  trailingPhiSignal.GetYaxis().CenterTitle()
  trailingPhiSignal.GetXaxis().SetTitleSize(0.055)
  trailingPhiSignal.GetYaxis().SetLabelSize(0.05)
  trailingPhiSignal.GetXaxis().SetLabelSize(0.05)
  trailingPhiSignal.GetYaxis().SetTitle("Arbitrary Units")
  trailingPhiSignal.GetXaxis().SetTitle("Trailing Lepton #phi")
  trailingPhiSignal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  trailingPhiSignal.SetTitleOffset(1.2,"Y")
  trailingPhiSignal.DrawNormalized('hist')
  trailingPhiDATA.DrawNormalized('epsame')
  trailingLeg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  trailingLeg.SetBorderSize(0)
  trailingLeg.SetTextSize(0.06)
  trailingLeg.SetFillColor(0)
  trailingLeg.AddEntry(trailingPhiSignal,'m_{H} = 125 GeV', 'l')
  trailingLeg.AddEntry(trailingPhiDATA, 'Data', 'ep')
  trailingLeg.Draw()
  trailingCan.SaveAs('trailingPhi_PAS.pdf')

  r9Can = TCanvas('r9Can','canvas',800,600)
  r9Can.cd()
  r9DATA.SetLineColor(kBlack)
  r9DATA.SetLineWidth(2)
  r9Signal.SetLineColor(kRed)
  r9Signal.SetLineWidth(2)
  r9Signal.GetYaxis().SetTitleSize(0.055)
  r9Signal.GetYaxis().CenterTitle()
  r9Signal.GetXaxis().SetTitleSize(0.055)
  r9Signal.GetYaxis().SetLabelSize(0.05)
  r9Signal.GetXaxis().SetLabelSize(0.05)
  r9Signal.GetYaxis().SetTitle("Arbitrary Units")
  r9Signal.GetXaxis().SetTitle("#gamma R9")
  r9Signal.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
  r9Signal.SetTitleOffset(1.2,"Y")
  r9Signal.DrawNormalized('hist')
  r9DATA.DrawNormalized('histsame')
  r9Leg = TLegend(0.54,0.65,0.75,0.89,'',"brNDC")
  r9Leg.SetBorderSize(0)
  r9Leg.SetTextSize(0.06)
  r9Leg.SetFillColor(0)
  r9Leg.AddEntry(r9Signal,'m_{H} = 125 GeV', 'l')
  r9Leg.AddEntry(r9DATA, 'Data', 'l')
  r9Leg.Draw()
  r9Can.SaveAs('r9_PAS.pdf')

  raw_input('Hit any key to continue')


def PTPlots2():
  mu2012FileMoriond = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_MoriondV4_2-7-13.root")
  mu2012FileICHEP = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_MuMu2012ABCD_ICHEPmzmzg_1-28-13.root")
  el2012FileMoriond = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_MoriondV4_2-7-13.root")
  el2012FileICHEP = TFile("/uscms/home/bpollack/work/CMSSW_5_3_6/src/PollackPrograms/HiggsZGAnalyzer/batchHistos/higgsHistograms_EE2012ABCD_ICHEPmzmzg_1-28-13.root")

  selection = ['el','mu']
  #selection = ['mu']
  catList = ['CAT1','CAT2','CAT3','CAT4']
  distList = ['Pt','Eta','Phi']
  physList = ['Photon','LeadingLepton','TrailingLepton']
  can= TCanvas('trailingCan','canvas',800,600)
  can.cd()
  gStyle.SetOptStat(0)
  for sel in selection:
    for cat in catList:
      for dist in distList:
        for phys in physList:
          if sel == 'mu':
            if phys == 'Photon':
              moriondPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_photon'+dist+cat+'_DATA')
              signalPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_photon'+dist+cat+'_Signal2012ggM120')
            if phys == 'LeadingLepton':
              moriondPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_leadLeptonStd'+dist+cat+'_DATA')
              signalPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_leadLeptonStd'+dist+cat+'_Signal2012ggM120')
            if phys == 'TrailingLepton':
              moriondPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_trailingLeptonStd'+dist+cat+'_DATA')
              signalPlot = mu2012FileMoriond.GetDirectory(cat).Get('h1_trailingLeptonStd'+dist+cat+'_Signal2012ggM120')
          else:
            if phys == 'Photon':
              moriondPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_photon'+dist+cat+'_DATA')
              signalPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_photon'+dist+cat+'_Signal2012ggM120')
            if phys == 'LeadingLepton':
              moriondPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_leadLeptonStd'+dist+cat+'_DATA')
              signalPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_leadLeptonStd'+dist+cat+'_Signal2012ggM120')
            if phys == 'TrailingLepton':
              moriondPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_trailingLeptonStd'+dist+cat+'_DATA')
              signalPlot = el2012FileMoriond.GetDirectory(cat).Get('h1_trailingLeptonStd'+dist+cat+'_Signal2012ggM120')


          moriondPlot.Scale(1/moriondPlot.Integral())
          signalPlot.Scale(1/signalPlot.Integral())
          moriondPlot.SetLineColor(kBlack)
          moriondPlot.SetLineWidth(2)
          moriondPlot.SetMarkerStyle(8)
          moriondPlot.SetMarkerColor(kBlack)
          signalPlot.SetLineColor(kRed)
          signalPlot.SetMarkerColor(kRed)
          signalPlot.SetLineWidth(2)
          signalPlot.SetMarkerStyle(9)
          moriondPlot.GetYaxis().SetTitleSize(0.05)
          moriondPlot.GetYaxis().CenterTitle()
          moriondPlot.GetXaxis().SetTitleSize(0.05)
          moriondPlot.GetYaxis().SetLabelSize(0.05)
          moriondPlot.GetXaxis().SetLabelSize(0.05)
          moriondPlot.GetYaxis().SetTitle("Events")
          moriondPlot.GetXaxis().SetTitle(sel+sel+' '+phys+' '+cat+' '+dist)
          moriondPlot.SetTitle("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}")
          moriondPlot.SetTitleOffset(1.2,"Y")
          moriondPlot.SetMinimum(signalPlot.GetMinimum()*0.9)
          moriondPlot.SetMaximum(max(signalPlot.GetMaximum()*1.1, moriondPlot.GetMaximum()*1.1))
          moriondPlot.Draw('hist')
          signalPlot.Draw('histsame')
          leg = TLegend(0.64,0.71,0.85,0.89,'',"brNDC")
          leg.SetBorderSize(0)
          leg.SetTextSize(0.03)
          leg.SetFillColor(0)
          leg.SetFillStyle(0)
          leg.AddEntry(signalPlot,'M120 Signal', 'l')
          leg.AddEntry(moriondPlot, 'DATA', 'l')
          leg.Draw()
          can.SaveAs(sel+'_'+dist+'_'+cat+'_'+phys+'.pdf')




def CATPlots():
  mu2012FileCAT1High = TFile("mu2012_CAT1_higgsOutput_Nominal_HighR9.root")
  mu2011FileCAT1High = TFile("mu2011_CAT1_higgsOutput_Nominal_HighR9.root")
  el2012FileCAT1High = TFile("el2012_CAT1_higgsOutput_Nominal_HighR9.root")
  el2011FileCAT1High = TFile("el2011_CAT1_higgsOutput_Nominal_HighR9.root")

  mu2012FileCAT2High = TFile("mu2012_CAT2_higgsOutput_Nominal_HighR9.root")
  mu2011FileCAT2High = TFile("mu2011_CAT2_higgsOutput_Nominal_HighR9.root")
  el2012FileCAT2High = TFile("el2012_CAT2_higgsOutput_Nominal_HighR9.root")
  el2011FileCAT2High = TFile("el2011_CAT2_higgsOutput_Nominal_HighR9.root")

  mu2012FileCAT3High = TFile("mu2012_CAT3_higgsOutput_Nominal_HighR9.root")
  mu2011FileCAT3High = TFile("mu2011_CAT3_higgsOutput_Nominal_HighR9.root")
  el2012FileCAT3High = TFile("el2012_CAT3_higgsOutput_Nominal_HighR9.root")
  el2011FileCAT3High = TFile("el2011_CAT3_higgsOutput_Nominal_HighR9.root")


  mu2012FileCAT1Low = TFile("mu2012_CAT1_higgsOutput_Nominal_LowR9.root")
  mu2011FileCAT1Low = TFile("mu2011_CAT1_higgsOutput_Nominal_LowR9.root")
  el2012FileCAT1Low = TFile("el2012_CAT1_higgsOutput_Nominal_LowR9.root")
  el2011FileCAT1Low = TFile("el2011_CAT1_higgsOutput_Nominal_LowR9.root")

  mu2012FileCAT2Low = TFile("mu2012_CAT2_higgsOutput_Nominal_LowR9.root")
  mu2011FileCAT2Low = TFile("mu2011_CAT2_higgsOutput_Nominal_LowR9.root")
  el2012FileCAT2Low = TFile("el2012_CAT2_higgsOutput_Nominal_LowR9.root")
  el2011FileCAT2Low = TFile("el2011_CAT2_higgsOutput_Nominal_LowR9.root")

  mu2012FileCAT3Low = TFile("mu2012_CAT3_higgsOutput_Nominal_LowR9.root")
  mu2011FileCAT3Low = TFile("mu2011_CAT3_higgsOutput_Nominal_LowR9.root")
  el2012FileCAT3Low = TFile("el2012_CAT3_higgsOutput_Nominal_LowR9.root")
  el2011FileCAT3Low = TFile("el2011_CAT3_higgsOutput_Nominal_LowR9.root")

  #mu2012List = [mu2012FileCAT1Low,mu2012FileCAT1High,mu2012FileCAT2Low,mu2012FileCAT2High,mu2012FileCAT3Low,mu2012FileCAT3High]
  #mu2011List = [mu2011FileCAT1Low,mu2011FileCAT1High,mu2011FileCAT2Low,mu2011FileCAT2High,mu2011FileCAT3Low,mu2011FileCAT3High]
  #el2012List = [el2012FileCAT1Low,el2012FileCAT1High,el2012FileCAT2Low,el2012FileCAT2High,el2012FileCAT3Low,el2012FileCAT3High]
  #el2011List = [el2011FileCAT1Low,el2011FileCAT1High,el2011FileCAT2Low,el2011FileCAT2High,el2011FileCAT3Low,el2011FileCAT3High]

  mu2012ListHigh = [mu2012FileCAT1High,mu2012FileCAT2High,mu2012FileCAT3High]
  mu2011ListHigh = [mu2011FileCAT1High,mu2011FileCAT2High,mu2011FileCAT3High]
  el2012ListHigh = [el2012FileCAT1High,el2012FileCAT2High,el2012FileCAT3High]
  el2011ListHigh = [el2011FileCAT1High,el2011FileCAT2High,el2011FileCAT3High]

  mu2012ListLow = [mu2012FileCAT1Low,mu2012FileCAT2Low,mu2012FileCAT3Low]
  mu2011ListLow = [mu2011FileCAT1Low,mu2011FileCAT2Low,mu2011FileCAT3Low]
  el2012ListLow = [el2012FileCAT1Low,el2012FileCAT2Low,el2012FileCAT3Low]
  el2011ListLow = [el2011FileCAT1Low,el2011FileCAT2Low,el2011FileCAT3Low]

  #catNames = ['CAT1Low','CAT1High','CAT2Low','CAT2High','CAT3Low','CAT3High']
  catNames = ['CAT1','CAT2','CAT3']

  #TotalList = {'mu2012':mu2012List,'mu2011':mu2011List,'el2012':el2012List,'el2011':el2011List}
  TotalList = {'mu2012High':mu2012ListHigh,'mu2011High':mu2011ListHigh,'el2012High':el2012ListHigh,'el2011High':el2011ListHigh,
      'mu2012Low':mu2012ListLow,'mu2011Low':mu2011ListLow,'el2012Low':el2012ListLow,'el2011Low':el2011ListLow}

  i = 1
  canList = []
  legList = []
  for selection in TotalList.keys():
    canList.append(TCanvas('c'+str(i),'canvas',800,600))
    canList[-1].cd()
    legList.append(TLegend(0.30,0.55,0.77,0.89,'',"brNDC"))
    legList[-1].SetBorderSize(0)
    legList[-1].SetTextSize(0.028)
    legList[-1].SetFillColor(0)
    bgList = []
    upList = []
    downList = []
    fresh = True
    for j,cat in enumerate(catNames):
      if '3' in cat:
        #bgList.append(TotalList[selection][j].Get(selection+'_expoFull_'+cat.strip('High').strip('Low')))
        bgList.append(TotalList[selection][j].Get(selection.strip('HighLow')+'_expoFull_'+cat))
      else:
        #bgList.append(TotalList[selection][j].Get(selection+'_pol3Full_'+cat.strip('High').strip('Low')))
        bgList.append(TotalList[selection][j].Get(selection.strip('HighLow')+'_pol3Full_'+cat))
      upList.append(TotalList[selection][j].Get(selection.strip('HighLow')+'_background_bgShapeUp'))
      downList.append(TotalList[selection][j].Get(selection.strip('HighLow')+'_background_bgShapeDown'))
      print bgList[-1].GetName()
      if 'CAT1' in cat:
        bgList[-1].SetLineColor(kBlue)
        upList[-1].SetLineColor(kBlue+2)
        upList[-1].SetLineStyle(2)
        downList[-1].SetLineColor(kBlue-8)
        downList[-1].SetLineStyle(2)
      elif 'CAT2' in cat:
        bgList[-1].SetLineColor(kRed)
        upList[-1].SetLineColor(kRed+1)
        upList[-1].SetLineStyle(2)
        downList[-1].SetLineColor(kRed-4)
        downList[-1].SetLineStyle(2)
      elif 'CAT3' in cat:
        bgList[-1].SetLineColor(kGreen+1)
        upList[-1].SetLineColor(kGreen+3)
        upList[-1].SetLineStyle(2)
        downList[-1].SetLineColor(kGreen-4)
        downList[-1].SetLineStyle(2)
      #if 'High' in cat: bgList[-1].SetLineStyle(1)
      #elif 'Low' in cat: bgList[-1].SetLineStyle(2)
      if fresh:
        bgList[-1].SetTitle(selection+' Backgrounds per Category')
        bgList[-1].SetMinimum(0)
        bgList[-1].SetMaximum(50)
        bgList[-1].Draw()
      else: bgList[-1].Draw("same")
      upList[-1].Draw("sameC")
      downList[-1].Draw("sameC")
      fresh = False
      if cat == 'CAT1' and ('Low' in selection):
        legTitle = 'Both Leptons in Barrel, #gamma in Barrel, R9<0.94'
      elif cat == 'CAT1' and ('High' in selection):
        legTitle = 'Both Leptons in Barrel, #gamma in Barrel, R9>0.94'
      elif cat == 'CAT2' and ('Low' in selection):
        legTitle = 'One or More Leptons in Endcap, #gamma in Barrel, R9<0.94'
      elif cat == 'CAT2' and ('High' in selection):
        legTitle = 'One or More Leptons in Endcap, #gamma in Barrel, R9>0.94'
      elif cat == 'CAT3' and ('Low' in selection):
        legTitle = '#gamma in Endcap, R9<0.94'
      elif cat == 'CAT3' and ('High' in selection):
        legTitle = '#gamma in Endcap, R9>0.94'
      legList[-1].AddEntry(bgList[-1],legTitle,'l')
    legList[-1].Draw()
    canList[-1].SaveAs(selection+'_catsPAS.pdf')
    i+=1



  raw_input('Hit any key to continue')

def CATSigPlots():
  mu2012FileCAT1High = TFile("mu2012_CAT1_higgsOutput_TEST.root")
  mu2011FileCAT1High = TFile("mu2011_CAT1_higgsOutput_TEST.root")
  el2012FileCAT1High = TFile("el2012_CAT1_higgsOutput_TEST.root")
  el2011FileCAT1High = TFile("el2011_CAT1_higgsOutput_TEST.root")

  mu2012FileCAT2High = TFile("mu2012_CAT2_higgsOutput_TEST.root")
  mu2011FileCAT2High = TFile("mu2011_CAT2_higgsOutput_TEST.root")
  el2012FileCAT2High = TFile("el2012_CAT2_higgsOutput_TEST.root")
  el2011FileCAT2High = TFile("el2011_CAT2_higgsOutput_TEST.root")

  mu2012FileCAT3High = TFile("mu2012_CAT3_higgsOutput_TEST.root")
  mu2011FileCAT3High = TFile("mu2011_CAT3_higgsOutput_TEST.root")
  el2012FileCAT3High = TFile("el2012_CAT3_higgsOutput_TEST.root")
  el2011FileCAT3High = TFile("el2011_CAT3_higgsOutput_TEST.root")


  mu2012FileCAT1Low = TFile("mu2012_CAT4_higgsOutput_TEST.root")
  mu2011FileCAT1Low = TFile("mu2011_CAT4_higgsOutput_TEST.root")
  el2012FileCAT1Low = TFile("el2012_CAT4_higgsOutput_TEST.root")
  el2011FileCAT1Low = TFile("el2011_CAT4_higgsOutput_TEST.root")

  mu2012FileCAT2Low = TFile("mu2012_CAT5_higgsOutput_TEST.root")
  mu2011FileCAT2Low = TFile("mu2011_CAT5_higgsOutput_TEST.root")
  el2012FileCAT2Low = TFile("el2012_CAT5_higgsOutput_TEST.root")
  el2011FileCAT2Low = TFile("el2011_CAT5_higgsOutput_TEST.root")

  mu2012FileCAT3Low = TFile("mu2012_CAT6_higgsOutput_TEST.root")
  mu2011FileCAT3Low = TFile("mu2011_CAT6_higgsOutput_TEST.root")
  el2012FileCAT3Low = TFile("el2012_CAT6_higgsOutput_TEST.root")
  el2011FileCAT3Low = TFile("el2011_CAT6_higgsOutput_TEST.root")

  mu2012List = [mu2012FileCAT1Low,mu2012FileCAT1High,mu2012FileCAT2Low,mu2012FileCAT2High,mu2012FileCAT3Low,mu2012FileCAT3High]
  mu2011List = [mu2011FileCAT1Low,mu2011FileCAT1High,mu2011FileCAT2Low,mu2011FileCAT2High,mu2011FileCAT3Low,mu2011FileCAT3High]
  el2012List = [el2012FileCAT1Low,el2012FileCAT1High,el2012FileCAT2Low,el2012FileCAT2High,el2012FileCAT3Low,el2012FileCAT3High]
  el2011List = [el2011FileCAT1Low,el2011FileCAT1High,el2011FileCAT2Low,el2011FileCAT2High,el2011FileCAT3Low,el2011FileCAT3High]

  catNames = ['CAT1Low','CAT1High','CAT2Low','CAT2High','CAT3Low','CAT3High']
  TotalList = {'mu2012':mu2012List,'mu2011':mu2011List,'el2012':el2012List,'el2011':el2011List}

  i = 1
  canList = []
  legList = []
  for selection in TotalList.keys():
    for j,cat in enumerate(catNames):
      corList = []
      uncorList = []
      canList.append(TCanvas('c'+str(i),'canvas',800,600))
      canList[-1].cd()
      gStyle.SetOptStat(1001100)
      legList.append(TLegend(0.20,0.75,0.47,0.90,'',"brNDC"))
      legList[-1].SetBorderSize(0)
      legList[-1].SetTextSize(0.028)
      legList[-1].SetFillStyle(0)
      uncorList.append(TotalList[selection][j].Get('uncorrected_M125'))
      uncorList[-1].SetTitle(cat+' '+selection+' '+'M125;M_{ll#gamma} (GeV);Entries')
      uncorList[-1].SetLineStyle(2)
      uncorList[-1].SetLineWidth(2)
      uncorList[-1].SetLineColor(kRed)
      corList.append(TotalList[selection][j].Get(selection+'_sigMC_M125'))
      if 'el' in selection:
        uncorList[-1].GetXaxis().SetRangeUser(115,145)
        uncorList[-1].GetYaxis().SetRangeUser(0,uncorList[-1].GetBinContent(uncorList[-1].GetMaximumBin())*1.2)
      else:
        uncorList[-1].GetXaxis().SetRangeUser(115,135)
      uncorList[-1].Draw("hist")
      gPad.Update()
      uncorStat = uncorList[-1].GetListOfFunctions().FindObject("stats")
      uncorStat.SetLineColor(kRed)
      uncorStat.SetLineWidth(2)
      uncorStat.SetLineStyle(2)
      uncorStat.SetY1NDC(uncorStat.GetY1NDC()-0.1)
      uncorStat.SetY2NDC(uncorStat.GetY2NDC()-0.1)
      uncorStat.SetX1NDC(uncorStat.GetX1NDC()-0.1)
      uncorStat.SetX2NDC(uncorStat.GetX2NDC()-0.1)
      corList[-1].Draw("histsames")
      gPad.Update()
      corStat = corList[-1].GetListOfFunctions().FindObject("stats")
      #print corStat.GetY1NDC(), corStat.GetY2NDC()
      corStat.SetLineWidth(2)
      corStat.SetY1NDC(corStat.GetY1NDC()-0.25)
      corStat.SetY2NDC(corStat.GetY2NDC()-0.25)
      corStat.SetX1NDC(corStat.GetX1NDC()-0.1)
      corStat.SetX2NDC(corStat.GetX2NDC()-0.1)
      corStat.Draw()
      gPad.Update()
      legList[-1].AddEntry(uncorList[-1],'uncorrected','l')
      legList[-1].AddEntry(corList[-1],'corrected','l')
      legList[-1].Draw()
      gPad.Update()
      canList[-1].SaveAs(selection+'_'+cat+'_M125.pdf')
      i+=1

  raw_input('Hit any key to continue')








if __name__=="__main__":
  if len(sys.argv) < 2:
    print 'choose: bgs, pts, cats, sigs'
  elif sys.argv[1] == 'bgs':
    BGPlots()
  elif sys.argv[1] == 'pts':
    PTPlots()
  elif sys.argv[1] == 'pts2':
    PTPlots2()
  elif sys.argv[1] == 'cats':
    CATPlots()
  elif sys.argv[1] == 'sigs':
    CATSigPlots()
  elif sys.argv[1] == 'roch':
    RochCorPlots()
  elif sys.argv[1] == 'acc':
    AccPlots()
  elif sys.argv[1] == 'rms':
    RMSComp()
  elif sys.argv[1] == 'bgc':
    BGComp()
  elif sys.argv[1] == 'table':
    RMSandFWHMandYieldTable()
  elif sys.argv[1] == 'ptrats':
    PtRatios()
  elif sys.argv[1] == 'accTable':
    AccTables()
  else:
    print 'choose: bgs, pts, cats, sigs, roch, acc, rms, bgc, table, ptrats, accTable'
