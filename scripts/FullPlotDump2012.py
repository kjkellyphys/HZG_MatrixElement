#!/usr/bin/env python
import sys
sys.argv.append('-b')
from ROOT import *
import os
gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
from PlotterV2 import *
import time
import math

MuLumiData2011 = 4.848
#MuLumiData2012 = 4.875
MuLumiData2012 = 5.135

ELumiData2011 = 4.634
#ELumiData2012 = 4.875
ELumiData2012 = 5.135

LumiData = 0




def FullPlotDump(savepath,rootfile,lep,period,scale = 'norm'):
  """Template for producing ZGamma histograms for different webpages and rootfiles"""

  if period == '2011':
    if lep.lower() == 'mumu': LumiData = MuLumiData2011
    elif lep.lower() == 'ee': LumiData = ELumiData2011
    HZG125ScaleX300 = LumiData*300/39956.7504536717
  elif period == '2012':
    if lep.lower() == 'mumu': LumiData = MuLumiData2012
    elif lep.lower() == 'ee': LumiData = ELumiData2012
    HZG125ScaleX300 = LumiData*300/31359.4988191727

  print HZG125ScaleX300

  selectionInfo = [
      #'MuMu2012 NWU nominal cut + Cos(theta) cut',
      #'scaling photon pT cut:'
      #'pT_g/llg_M > 15/110',
      #'pT_g > 15 GeV'
      #'PU test'
      #'No PU Weights',
      #'R9 > 0.9',
      #'New nominal',
      'PU obs Test'
      ]

  cutList = ['Initial',
      'Electron Presel',
      'HLT',
      'Primary VTX',
      'Data Qual Bits',
      '2 ID leptons',
      '2 Iso leptons',
      'Gamma ID/Iso/pT',
      'Gamma DR Cuts',
      'Good Z (eeGamma)',
      'Z Mass',
      '3rd lep veto',
      'bjet veto',
      'MET',
      'MET/QT',
      'Jet Multi',
      '3body Mass Low',
      '3body Mass High',
      'ZG pT diff',
      'Angle Cuts',
      'M12 Cut',
      'MVA Cuts',
      'VBF Cuts',
      ]

  initWebPage(savepath)

  f1 = TFile('../HiggsZGAnalyzer/batchHistos/higgsHistograms_'+rootfile+'.root', 'OPEN')
  h1_acceptanceByCut = getHist(f1, 'acceptanceByCut',directory = 'Misc')
  h1_acceptanceByCutRaw = getHist(f1, 'acceptanceByCutRaw',directory = 'Misc')
  h1_invMass = getHist(f1, 'InvariantMass',directory = 'ZGamma',clone=True)
  h1_diLeptonMass = getHist(f1, 'diLeptonMassPreSel', directory = 'PreSel')
  axis = h1_diLeptonMass[-1].GetXaxis()
  ZJetsFrac = h1_diLeptonMass[-1].Integral(axis.FindBin(81.2),axis.FindBin(101.2))/h1_diLeptonMass[0].Integral(axis.FindBin(81.2),axis.FindBin(101.2))
  #HiggsFrac = h1_diLeptonMass[-1].Integral(axis.FindBin(81.2),axis.FindBin(101.2))/h1_diLeptonMass[-2].Integral(axis.FindBin(81.2),axis.FindBin(101.2))
  HiggsFrac = 1


  acc = makeYields(h1_acceptanceByCut, h1_acceptanceByCutRaw, h1_invMass, cutList, lep, period,
      #,scale = {'Signal2011ggM125':1}
      #,scale = {'ZJets':ZJetsScale}
      #,scale = {'ZJets':ZJetsFrac,'HZG125Signal':1}
    )


  cname = ['ZJets',
      'ZGToLLG',
  #    'WZJets',
  #    'ZZJets',
  #    'WWJets',
  #    'ttbar',
      'ggHZG125x300',
      'vbfHZG125',
      'DATA']


  allFolders = getAllFolders(f1)
  #presPath = savepath+'/Presentation'
  #initWebPage(presPath)
  for folder in allFolders:
    if folder in ["TESTS","METCor"]: continue
    if folder.find('Iso') != -1: continue
    deepPath = savepath+'/'+folder
    initWebPage(deepPath)
    if folder == 'GenLvl':
      for key in allFolders[folder]:
        for i,hist in enumerate(allFolders[folder][key]):
          ppps(deepPath,key+'_'+str(i),hist,dopt='hist')

    else:
      for key in allFolders[folder]:
        if key in["TESTS"]: continue
        if allFolders[folder][key][0].ClassName()[0:3] in ['TH1']:
          print HZG125ScaleX300
          scaleHists(allFolders[folder][key],{'Signal2012ggM125':HZG125ScaleX300,'Signal2011ggM125':HZG125ScaleX300}, lep = lep, period = period)

          stackList = makeStackGroup(allFolders[folder][key],"125")
          #stackList = sorted(stackList,key = lambda x:x.Integral(),reverse=True)

          if len(stackList)>1:
            if key in ['pvMultPreSel']:
              stacksAndRatio(deepPath,key+'_stack',stackList,logToggle = True)
              #stacksAndRatio(deepPath,key+'_stack',stackList,cname = cname,logToggle = True)
            else:
              stacksAndRatio(deepPath,key+'_stack',stackList,cname = cname, logToggle = False)
              #stacksAndRatio(deepPath,key+'_stack',stackList,cname = cname,logToggle = True)
          '''
          for hist in allFolders[folder][key]:
            name = dataName(hist.GetName())
            if name in ['DATA','Signal2011ggM125','Signal2012ggM125']:
              ppps(deepPath,key+name,hist,legOn = True)

          '''
        elif allFolders[folder][key][0].ClassName()[0:3] == 'TH2':
          for hist in allFolders[folder][key]:
            name = dataName(hist.GetName())
            if name in ['DATA','ZJets','HZG125Signal']:
              ppps(deepPath,key+name,hist,legOn = True)

    makeWebPageBetter(deepPath)


  h1_invariantMassCut = getHist(f1, 'InvariantMassCut', dataType = 'DATA', directory = 'ZGamma')
  h1_invariantMassSignal = getHist(f1, 'InvariantMassCut', dataType = 'Signal2011ggM125', directory = 'ZGamma')
  fit = makeFit(h1_invariantMassCut,exclude = True, excludeRegion = [120,130],xmin=113)
  #fit = makeFit(h1_invariantMassCut,exclude = False, fitType ='expo', xmin=113)
#  fitHigs = makeFit(h1_invariantMassSignal)
  ppps(savepath+'/ZGamma','superFits',h1_invariantMassCut,legOn=False)
  gStyle.SetOptStat(0)
  #gStyle.SetOptFit(1)
  #ppps(savepath+'/ZGamma','signalStandAlone',h1_invariantMassSignal, legOn=False)
  makeWebPageBetter(savepath+'/ZGamma')

  #makeWebPageBetter(presPath)
  makeWebYeild(savepath+'/Yields',acc)
  makeWebPageBest(savepath,selectionInfo,'H->Zg->'+lep+' g')

  for hist in allFolders['ZGamma']['InvariantMass']:
    err = Double(0)
    print hist.GetName(), hist.IntegralAndError(hist.FindBin(122),hist.FindBin(127),err),err
    #for x in range(1, hist.GetNbinsX()):
    #  print x, hist.GetBinLowEdge(x)
    #print 'find 125!', hist.FindBin(125)
  print "fit result:",fit.Integral(122,128)/2.,fit.IntegralError(122,128)/2.


if __name__=="__main__":
  if len(sys.argv)>1:
    FullPlotDump(*sys.argv[1:-1])
  else:
    print 'WTF you need to specify a savepath and a goddamn root file!'

