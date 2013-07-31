#!/usr/bin/env python
# coding: utf-8
from ROOT import *
import os
import shutil
import time
import string
from collections import defaultdict
import datetime
import codecs

MuLumiData2011 = 0
MuLumiData2011 += 2.098
MuLumiData2011 += 2.750

#MuLumiData2012 = 0
#MuLumiData2012 += 0.682875
#MuLumiData2012 += 0.076248
#MuLumiData2012 += 4.116
#MuLumiData2012 = 4.875
MuLumiData2012 = 5.135


ELumiData2011 = 0
ELumiData2011 += 1.972
ELumiData2011 += 2.681

#ELumiData2012 = 0
#ELumiData2012 += 0.682917
#ELumiData2012 += 0.076248
#ELumiData2012 += 4.116
ELumiData2012 = 5.135


LumiData = 0

#'ZJets':36247791
iniEventDict = {'ZJets':25993422,'ZGMuMu':321534,'HZG135Signal':50000,'WWJets':1197558,'ZZJets':1103468,'WJets':81251000,'WZJets':1221134}
mainLumiDict = {'Aug05':370.9,'May10':215.1,'PromptV4':927.6,'PromptV6':663.0,'2011B':2511}


def getHist(tfile, var, dataType=None,directory=None,tcls = None,clone=False):
  '''Get histogram from ROOT file, or list of hists for multiple datatypes'''
  if tcls == None: tcls = 'h1'
  if dataType != None:
    if dataType != '': dataType = '_' + dataType
    histogramName = tcls + '_' + var + dataType
    if directory==None:
      if clone:
        hist = (tfile.Get(histogramName)).Clone()
      else: hist = tfile.Get(histogramName)
    else:
      if clone:
        hist = (tfile.GetDirectory(directory).Get(histogramName)).Clone()
      else:
        hist = tfile.GetDirectory(directory).Get(histogramName)

  elif dataType == None:
    var = '_'+var+'_'
    hist = []
    if directory == None:
      lok = tfile.GetListOfKeys()
      for i in range(0, lok.GetEntries()):
        name = lok.At(i).GetName()
        if string.find(name,var) != -1:
          if clone:
            hist.append((tfile.Get(name)).Clone())
          else:
            hist.append(tfile.Get(name))
    else:
      lok = tfile.GetDirectory(directory).GetListOfKeys()
      for i in range(0, lok.GetEntries()):
        name = lok.At(i).GetName()
        if string.find(name,var) != -1:
          if clone:
            hist.append((tfile.GetDirectory(directory).Get(name)).Clone())
          else:
            hist.append(tfile.GetDirectory(directory).Get(name))

  return hist

def getFolder(tfile, directory):
  '''Get a directory from ROOT file, return dictionary of all histos in that directory, grouped by hist type'''
  lok = tfile.GetDirectory(directory).GetListOfKeys()
  folder = defaultdict(list)
  for i in range(0, lok.GetEntries()):
    name = lok.At(i).GetName()
    hist = tfile.GetDirectory(directory).Get(name)
    key = (hist.GetName().strip(string.letters+'0123456789')).strip('_')
    folder[key].append(hist)
  return folder

def getAllFolders(tfile):
  '''Return a dictionary of folder dictionaries'''
  allFolders = {}
  lok = tfile.GetListOfKeys()
  for i in range(lok.GetEntries()):
    name = lok.At(i).GetName()
    if name[0:2] not in ['h1','h2']:
      print name
      allFolders[name] = getFolder(tfile,name)
  return allFolders

def getPUHists(f1,pathname,suffix, dirc = None):
  h1_PU = []
  h1_PU_F = []
  for i in range (1,20):
    histTemp = getHist(f1,pathname+'_'+str(i),suffix,directory = dirc)
    if histTemp != None:
      histTempF = histTemp.Clone()
      h1_PU.append(histTemp)
      h1_PU_F.append(histTempF)
      histTemp.Scale(histTemp.Integral()**(-1))
      print h1_PU[i-1].GetName()
    else:
      h1_PU.append(None)
      h1_PU_F.append(None)
  return [h1_PU,h1_PU_F]

def scaleHists(histList,userScale = {},lep = 'mumu',period = '2012'):

  if period == '2012':
    if lep.lower() == 'mumu': LumiData=MuLumiData2012
    elif lep.lower() == 'ee': LumiData=ELumiData2012
  elif period == '2011':
    if lep.lower() == 'mumu': LumiData=MuLumiData2011
    elif lep.lower() == 'ee': LumiData=ELumiData2011
  else: print 'OH SHIT'

  if type(histList) == type(dict()): None
  for hist in histList:
    name = dataName(hist.GetName())
    if name in userScale:
      hist.Scale(userScale[dataName(hist.GetName())])
  for hist in histList:
    name = dataName(hist.GetName())
    if name not in userScale:
      #if name.find('DYMuMu') > -1:          hist.Scale(LumiData/17.827)
      if name.find('DYMuMu') > -1:          hist.Scale(LumiData/11.9)
      elif name.find('ZGMuMu2j') > -1:   hist.Scale(LumiData/38.55)
      elif name.find('ZGMuMu') > -1:   hist.Scale(LumiData/7.11)
      elif name.find('ZgammaMuMu') > -1:   hist.Scale(LumiData/27.05)
      elif name.find('ZZJets') > -1:   hist.Scale(LumiData/3678.2)
      elif name.find('WZJets') > -1:    hist.Scale(LumiData/1426.6)
      #elif name.find('ZJets2011') > -1:      hist.Scale(LumiData/11.900)
      elif name.find('ZJetsBG2012') > -1:      hist.Scale(LumiData/7.41882804227519)
      elif name.find('ZJetsBG2011') > -1:      hist.Scale(LumiData/8.52802559055118)
      elif name.find('ZGToLLG') > -1:      hist.Scale(LumiData/37.394479638009)
      #elif name.find('ZJets') > -1:      hist.Scale(LumiData/9.520)
      elif name.find('WWJets') > -1:   hist.Scale(LumiData/250.38)
      elif name.find('WJets') > -1:   hist.Scale(LumiData/2.595)
      elif name.find('TTJets') > -1:   hist.Scale(LumiData/LumiData)
      elif name.find('ttbar') > -1:   hist.Scale(LumiData/618.754)
      elif name.find('tW') > -1:   hist.Scale(LumiData/103.2)
      elif name.find('GluGluWW') > -1:   hist.Scale(LumiData/730.831)
      elif name.find('GammaJet15-30') > -1:   hist.Scale(LumiData/0.012)
      elif name.find('GammaJet30-50') > -1:   hist.Scale(LumiData/0.131)
      elif name.find('GammaJet50-80') > -1:   hist.Scale(LumiData/0.748)
      elif name.find('GammaJet80-120') > -1:   hist.Scale(LumiData/4.577)
      elif name.find('GammaJet120-170') > -1:   hist.Scale(LumiData/24.810)
      elif name.find('GammaJet170-300') > -1:   hist.Scale(LumiData/91.394)
      elif name.find('GammaJet300-470') > -1:   hist.Scale(LumiData/1391.078)
      elif name.find('GammaJet470-800') > -1:   hist.Scale(LumiData/15812.212)
      elif name.find('OFFggHZG125') > -1:     hist.Scale(LumiData/40114.2)
      elif name.find('ggHZG125') > -1:     hist.Scale(LumiData/40454.268)
      elif name.find('OFFvbfHZG125') > -1:     hist.Scale(LumiData/532164.83)
      elif name.find('vbfHZG125') > -1:     hist.Scale(LumiData/574187.037)
      elif name.find('HZG125') > -1:     hist.Scale(LumiData/124532.73)
      elif name.find('HZG125') > -1 and name.find('OLD')>-1:     hist.Scale(LumiData/31447.658)
      elif name.find('HZG135') > -1 and name.find('v2')==-1:        hist.Scale((LumiData/1632.90))
      elif name.find('HZG135') > -1 and name.find('v2')>-1:        hist.Scale((LumiData/13553.056))
      elif name.find('Signal2012ggM120') > -1:     hist.Scale(LumiData/42218.1040760702)
      elif name.find('Signal2012ggM125') > -1:     hist.Scale(LumiData/31359.4988191727)
      elif name.find('Signal2012ggM130') > -1:     hist.Scale(LumiData/28101.993226726)
      elif name.find('Signal2012ggM135') > -1:     hist.Scale(LumiData/25981.6282816437)
      elif name.find('Signal2012ggM140') > -1:     hist.Scale(LumiData/25752.4427312811)
      elif name.find('Signal2012vbfM120') > -1:     hist.Scale(LumiData/540407.382343672)
      elif name.find('Signal2012vbfM125') > -1:     hist.Scale(LumiData/407060.467517419)
      elif name.find('Signal2012vbfM130') > -1:     hist.Scale(LumiData/335758.254809396)
      elif name.find('Signal2012vbfM135') > -1:     hist.Scale(LumiData/301017.825345686)
      elif name.find('Signal2012vbfM140') > -1:     hist.Scale(LumiData/289477.298318026)

      elif name.find('Signal2011ggM120') > -1:     hist.Scale(LumiData/53577.6900376794)
      elif name.find('Signal2011ggM125') > -1:     hist.Scale(LumiData/39956.7504536717)
      elif name.find('Signal2011ggM130') > -1:     hist.Scale(LumiData/35861.7950287386)
      elif name.find('Signal2011ggM135') > -1:     hist.Scale(LumiData/33274.7169221051)
      elif name.find('Signal2011ggM140') > -1:     hist.Scale(LumiData/33046.8538497475)
      elif name.find('Signal2011vbfM120') > -1:     hist.Scale(LumiData/702231.499987955)
      elif name.find('Signal2011vbfM125') > -1:     hist.Scale(LumiData/530422.310274556)
      elif name.find('Signal2011vbfM130') > -1:     hist.Scale(LumiData/439628.009546792)
      elif name.find('Signal2011vbfM135') > -1:     hist.Scale(LumiData/396248.919182322)
      elif name.find('Signal2011vbfM140') > -1:     hist.Scale(LumiData/382209.094452222)

def makePretty(hist, title=None, xtitle=None,ytitle=None):
  '''Make the histo presentable to the unwashed masses'''
  if title!=None: hist.SetTitle(title)
  if xtitle!=None: xtitle=hist.GetXaxis().SetTitle(xtitle)
  if ytitle!=None: ytitle=hist.GetYaxis().SetTitle(ytitle)
  pt = TPaveText(0.2,0.925,0.85,0.995,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(11)
  text = pt.AddText("")
  pt.Draw()
  if hist.ClassName()[0:3] == 'TH2': hist.SetOption("COLZ"); hist.SetStats(0)
  #if hist.ClassName()[0:3] == 'TH1': hist.SetMarkerStyle(20)
  return pt

def makeLegend(hist,leg = None, lstr = None, rstr = None, cname = None,
    x1 = None, x2 = None, y1 = None, y2 = None):
  '''Input a hist or list of hists, names or list of names (or none and it'll try to figure it out), and an old legend (or none and it'll make one), get a kickass legend in return'''
  if type(hist)==list:
    if leg == None:
      if x1 == None: x1 = 0.77
      if x2 == None: x2 = 1.0
      if y1 == None: y1 = 0.93-len(hist)*0.07
      if y2 == None: y2 = 0.93

      #leg = TLegend(0.77,0.93-len(hist)*0.05,1.0,0.93)
      leg = TLegend(x1,y1,x2,y2)
      leg.SetFillColor(0)
      leg.SetBorderSize(1)
      leg.SetShadowColor(0)
    else: leg.SetY1NDC(leg.GetY1NDC()-len(hist)*0.04)
    for i,x in enumerate(hist):
      if cname == None and (lstr == None or rstr == None):
        cname = x.GetName()
        cname = cname.lstrip(string.letters+'0123456789')
        cname = cname.lstrip('_')
        cname = cname.lstrip(string.letters+'0123456789')
        cname = cname.lstrip('_')
        if x !=None: leg.AddEntry(x, cname)
        cname = None
      elif cname == None:
        cname = x.GetName()
        cname = cname.lstrip(lstr)
        cname = cname.rstrip(rstr)
        if x !=None: leg.AddEntry(x, cname)
        cname = None
      else:
        cname[i] = cname[i].strip()
        if x !=None: leg.AddEntry(x, cname[i])
  else:
    if leg == None:
      if x1 == None: x1 = 0.77
      if x2 == None: x2 = 1.0
      if y1 == None: y1 = 0.88
      if y2 == None: y2 = 0.93
      # leg = TLegend(0.51,0.89,0.76,0.93)
      leg = TLegend(x1,y1,x2,y2)
      leg.SetFillColor(0)
      leg.SetBorderSize(1)
      leg.SetShadowColor(0)
    else:
      newSize = leg.GetY1NDC()-0.04
      leg.SetY1NDC(newSize)

    if cname == None and (lstr == None or rstr == None):
      cname = hist.GetName()
      cname = cname.lstrip(string.letters+'0123456789')
      cname = cname.lstrip('_')
      cname = cname.lstrip(string.letters+'0123456789')
      cname = cname.lstrip('_')
    elif cname == None:
      cname = hist.GetName()
      cname = cname.lstrip(lstr)
      cname = cname.rstrip(rstr)
    leg.AddEntry(hist, cname)
  return leg




def makeFit(hist,fitType = None, xmin = None, xmax = None, exclude = False,excludeRegion = [0,0]):
  '''Fit a histo! choose from our selection of many excellent fits!'''
  def fline(x, par):
    if (exclude and x[0] > excludeRegion[0] and x[0] < excludeRegion[1]):
      TF1.RejectPoint()
      return 0
    if fitType == 'expo' or fitType == None:
      return exp(par[0]+par[1]*x[0])
    elif fitType == 'pol3':
      return par[0]+par[1]*x[0]+par[2]*x[0]**2+par[3]*x[0]**3


  fitTypeFull = ''
  if xmin == None: xmin = hist.GetXaxis().GetBinLowEdge(1)
  if xmax == None: xmax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
  if exclude:
    fit = TF1(hist.GetName()+"_excFit",fline,xmin,xmax,2)
    fit.SetParameters(8,-0.03)
    hist.Fit(fit,"0")
    exclude = False
    fitLeft = TF1("fitLeft",fline,xmin,excludeRegion[0],2)
    fitLeft.SetParameters(fit.GetParameters())
    hist.GetListOfFunctions().Add(fitLeft)
    fitRight = TF1("fitRight",fline,excludeRegion[1],xmax,2)
    fitRight.SetParameters(fit.GetParameters())
    hist.GetListOfFunctions().Add(fitRight)
    return fit







  if fitType == None or fitType == 'gaus':
    fitTypeFull = 'gaus'
    fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
    #return fit1
#    fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),6)
  elif fitType == 'novo':
    fitTypeFull = "[0]*exp(-0.5*(log(1+[3]*((x-[1])/[2])*sinh([3]*sqrt(log(4)))/([3]*sqrt(log(4))))**2/[3]**2+[3]**2))"
    #fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
    fit = TF1('novo',fitTypeFull,xmin,xmax)
    #print hist.GetName()
    #time.sleep(2)
    if hist.GetName().find('_ZZJetsTo2L2Nu')!=-1:
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),25,0.5)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('TrackMET0') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),14,0.5)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('CorTrackMET2vsVTX_2') !=-1):
      fit.SetParameters(34,18,0.6)
      fit.FixParameter(2,18)
      fit.FixParameter(1,34)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('CorTrackMET2') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),22,0.6)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('CorReducedMET1') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),15,0.7)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('CorProjLep') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),10,0.5)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('Cor') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),15,0.6)
    elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('Proj') !=-1):
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),15,0.5)
    elif hist.GetName().find('ggH')!=-1:
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),25,0.5)
    else:
      fit.SetParameters(hist.GetBinContent(hist.GetMaximumBin()),hist.GetMaximumBin(),6,0.5)
    fit.SetParNames("Height","Peak","Width","Tail");
  elif fitType == 'peak':
    fitTypeFull = "[0]+(x-1)*[1]"
    fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
    fit.SetParNames('a_{0}','a_{PU}')
    fit.SetLineColor(3)
    fit.SetLineWidth(2)
  elif fitType == 'width':
    fitTypeFull = "sqrt([0]^2+(x-1)*[1]^2)"
    fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
    fit.SetParNames('#sigma_{0}','#sigma_{PU}')
    fit.SetLineColor(3)
    fit.SetLineWidth(2)
  elif fitType =="pol1Again":
    fitTypeFull = 'pol1'
    fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
  else:
    fitTypeFull = fitType
    fit = TF1(hist.GetName()+'_fit',fitTypeFull,xmin,xmax)
  if (fitType == 'peak' or fitType == 'width'):
    hist.Fit(fit.GetName(),'R','sames')
  elif fitType == 'pol1Again':
    hist.Fit(fit.GetName(),'R+')
  elif hist.GetName().find('GluGlu')!=-1:
    hist.Fit(fit.GetName(),'RML')
  elif (hist.GetName().find('ggH')!=-1 and hist.GetName().find('CorTrackMET2vsVTX_2') !=-1):
    hist.Fit(fit.GetName(),'RME')
  else: hist.Fit(fit.GetName(),'RM')
  return fit


def makeStackGroup(histList,selection,alt = False,kill = ['asdfgadfgbadfbsdfb']):
  '''feed in a histlist, discard unwanted histograms with irrelevant selections, put all the BG in front, signal second to last, data last'''
  stackGroup = []
  if alt == False:
    for hist in histList:
      name = dataName(hist.GetName())
      if name in kill: continue
      if hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find('Alt') == -1:
        stackGroup.append(hist)
  elif alt == True:
    for hist in histList:
      if name in kill: continue
      if hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find('ZGMuMu') != -1:
        stackGroup.append(hist)
  else:
    for hist in histList:
      if name in kill: continue
      #print dataName(hist.GetName())
      if hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find('_'+alt) != -1:
        stackGroup.append(hist)
  for hist in histList:
    if name in kill: continue
    if hist.GetName().find("Signal") != -1 and hist.GetName().find(selection) != -1:
      stackGroup.append(hist)
  for hist in histList:
    if name in kill: continue
    if hist.GetName().find('DATA') != -1:
      stackGroup.append(hist)

  return stackGroup

def ppps(path, filename, hist, title=None, xtitle=None,ytitle=None, xmin=None, xmax=None,ymin=None, ymax=None, lstr = None, rstr = None, extraObj = None,legOn = True, noLine = True,dopt = '', comp = False,cname = None, logy = False):
  '''Power Plot Pro: Stacks, one stop shop for all your plotting needs'''
  tf1Bool = False
  canvas = TCanvas('canvas','canvas',800,600)
  canvas.SetGrid(1,1)
  if logy: canvas.SetLogy(kTRUE)
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  if not os.path.isdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path):
    os.mkdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path)

  # single hist function

  if type(hist)!=list:
    if (xmin and xmax) !=None: hist.SetAxisRange(xmin,xmax)
    if (ymin and ymax) !=None: hist.GetYaxis().SetRange(ymin,ymax)
    hist.SetMarkerStyle(20)
    hist.SetOption('EP')
    hist.Draw(dopt)
    pt = makePretty(hist, title, xtitle, ytitle)
    if not noLine:
      if hist.GetListOfFunctions().GetEntries()>0 and hist.GetListOfFunctions().At(0).ClassName()=='TF1':
          hist.GetListOfFunctions().At(0).SetLineWidth(2)
          hist.GetListOfFunctions().At(0).SetLineColor(kBlack)
      else: hist.Draw('C hist same')
    if legOn:
      leg = makeLegend(hist,lstr = lstr, rstr = rstr,y1=0.75,y2=0.80)
      leg.Draw()
    if extraObj != None: extraObj.Draw()
    canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'.png')

  # list of hists function

  else:
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    histcopy = []
    for x in hist:
      clone = x.Clone()
      if clone.Integral() != 0: clone.Scale(1/clone.Integral())
      histcopy.append(clone)
    maxbinval = 0
    maxbinhist = 0
    for i,x in enumerate(histcopy):
      #print x.GetName()
      if x!=None:
        if x.GetBinContent(x.GetMaximumBin())>maxbinval:
          maxbinval = x.GetBinContent(x.GetMaximumBin())
          maxbinhist = i
          print x.GetName()
          name = dataName(x.GetName())
    if (xmin and xmax) !=None: histcopy[maxbinhist].SetAxisRange(xmin,xmax)
    if (ymin and ymax) !=None: histcopy[maxbinhist].GetYaxis().SetRangeUser(ymin,ymax)
    if maxbinhist == len(histcopy)-1 and comp:
      pt = makePretty(histcopy[maxbinhist], title, xtitle, ytitle)
      histcopy[maxbinhist].Draw(dopt)
      #histcopy[maxbinhist].SetLineWidth(4)
      if histcopy[maxbinhist].GetListOfFunctions().GetEntries()>0:
        if histcopy[maxbinhist].GetListOfFunctions().At(0).ClassName()=='TF1':
          histcopy[maxbinhist].GetListOfFunctions().At(0).SetLineWidth(4)
          histcopy[maxbinhist].GetListOfFunctions().At(0).SetLineColor(kBlack)
        else: histcopy[maxbinhist].Draw('C hist same')
    else:
      #histcopy[maxbinhist].SetFillColor(ColorDict[maxbinhist])
      #histcopy[maxbinhist].SetFillStyle(FillList[maxbinhist])
      histcopy[maxbinhist].Draw(dopt)
      pt = makePretty(histcopy[maxbinhist], title, xtitle, ytitle)
      histcopy[maxbinhist].SetLineColor(ColorDict[name])
    #  histcopy[maxbinhist].SetLineStyle(9)
      histcopy[maxbinhist].SetLineWidth(2)
      if dopt.find('P') != -1:
        histcopy[maxbinhist].SetMarkerStyle(MarkerList[maxbinhist])
        histcopy[maxbinhist].SetMarkerColor(ColorDict[name])
        histcopy[maxbinhist].SetMarkerSize(2)
      if histcopy[maxbinhist].GetListOfFunctions().GetEntries()>0:
        if histcopy[maxbinhist].GetListOfFunctions().At(0).ClassName()=='TF1':
          histcopy[maxbinhist].GetListOfFunctions().At(0).SetLineWidth(2)
          histcopy[maxbinhist].GetListOfFunctions().At(0).SetLineColor(ColorDict[maxbinhist])
        else: histcopy[maxbinhist].Draw('C hist same')
    for i,x in enumerate(histcopy):
      name = dataName(x.GetName())
      if i==len(histcopy)-1 and i!=maxbinhist and comp:
        pt = makePretty(x, title, xtitle, ytitle)
        x.Draw('sames'+dopt)
        if dopt.find('P') != -1:
          x.SetMarkerColor(ColorDict[i])
          x.SetMarkerStyle(MarkerList[i])
          x.SetMarkerSize(2)
        if x.GetListOfFunctions().GetEntries()>0:
          if x.GetListOfFunctions().At(0).ClassName()=='TF1':
            x.GetListOfFunctions().At(0).SetLineWidth(2)
            x.GetListOfFunctions().At(0).SetLineColor(ColorDict[i])
          else: x.Draw('C hist same')
      elif x!=None and i!=maxbinhist:
        #x.SetFillColor(ColorDict[i])
        #x.SetFillStyle(FillList[i])
        x.Draw('sames'+dopt)
        pt = makePretty(x, title, xtitle, ytitle)
        x.SetLineColor(ColorDict[name])
        #x.SetLineStyle(9)
        x.SetLineWidth(2)
        if dopt.find('P') != -1:
          x.SetMarkerColor(ColorDict[i])
          x.SetMarkerStyle(MarkerList[i])
          x.SetMarkerSize(2)
        if x.GetListOfFunctions().GetEntries()>0:
          if x.GetListOfFunctions().At(0).ClassName()=='TF1':
            x.GetListOfFunctions().At(0).SetLineWidth(2)
            x.GetListOfFunctions().At(0).SetLineColor(ColorDict[i])
          else: x.Draw('C hist same')
    if legOn:
      if cname != None:
        leg = makeLegend(histcopy,cname = cname,y1 =0.70,y2=0.75 )
      else:
        leg = makeLegend(histcopy,lstr = lstr, rstr = rstr)
    for x in histcopy: x.GetYaxis().SetTitleOffset(1.1)
    leg.Draw()
    canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'.png')
    return histcopy


def makeYields(histlistWeight, histlistRaw, histlistMass, cutlist, lep, period, scale = {}):
  dataDictWeight = {}
  bgDictWeight = {}
  signalDictWeight = {}
  extraDictWeight = {}
  dataDictRaw = {}
  bgDictRaw = {}
  signalDictRaw = {}
  extraDictRaw = {}
  dataDictMass = {}
  bgDictMass = {}
  signalDictMass = {}
  extraDictMass = {}
  dataDictMassRaw = {}
  bgDictMassRaw = {}
  signalDictMassRaw = {}
  extraDictMassRaw = {}

  totalbg = 0
  totalbgRaw = 0
  totalbgerr = 0
  totalbgCol = []
  totalbgRawCol = []
  totalbgerrCol = []

  totalsig = 0
  totalsigRaw = 0
  totalsigerr = 0
  totalsigCol = []
  totalsigRawCol = []
  totalsigerrCol = []

  oldcut = []
  newcut = []
  tmpList = []

  bgcol = True
  binomialerror = 0

  for hist in histlistMass:
    tmpList.append(hist.Clone())
  scaleHists(histlistMass, scale, lep, period)
  scaleHists(histlistWeight, scale, lep, period)

  table = '<table border = "10"\ncellpadding="5">\n<tr>\n'
  table = table +'<th></th>\n'

  ### Fill the dictionaries with DATA, BG, Signal, or Extra
  for hist in histlistWeight:
    if dataName(hist.GetName()).find("DATA") != -1: dataDictWeight[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("BG") != -1: bgDictWeight[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("Signal") != -1: signalDictWeight[dataName(hist.GetName())] = hist
    else: extraDictWeight[dataName(hist.GetName())] = hist
  for hist in histlistRaw:
    if dataName(hist.GetName()).find("DATA") != -1: dataDictRaw[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("BG") != -1: bgDictRaw[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("Signal") != -1: signalDictRaw[dataName(hist.GetName())] = hist
    else: extraDictRaw[dataName(hist.GetName())] = hist
  for hist in histlistMass:
    if dataName(hist.GetName()).find("DATA") != -1: dataDictMass[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("BG") != -1: bgDictMass[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("Signal") != -1: signalDictMass[dataName(hist.GetName())] = hist
    else: extraDictMass[dataName(hist.GetName())] = hist
  for hist in tmpList:
    if dataName(hist.GetName()).find("DATA") != -1: dataDictMassRaw[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("BG") != -1: bgDictMassRaw[dataName(hist.GetName())] = hist
    elif dataName(hist.GetName()).find("Signal") != -1: signalDictMassRaw[dataName(hist.GetName())] = hist
    else: extraDictMassRaw[dataName(hist.GetName())] = hist

  for name,value in sorted(bgDictWeight.iteritems(), key = lambda (k,v):(v,k)):
    table = table +'<th>{0}</th>\n'.format(name)
  for name,value in sorted(signalDictWeight.iteritems(), key = lambda (k,v):(v,k)):
    table = table +'<th>{0}</th>\n'.format(name)
  for name in sorted(dataDictWeight.iteritems(), key = lambda (k,v):(v,k)):
    table = table +'<th>{0}</th>\n'.format(name)


  # Do we need to sum some columns at the end?
  doTotalBG = False
  doTotalSig = False
  if len(bgDictWeight) > 0:
    table = table +'<th>{0}</th>\n'.format("Total BG")
    doTotalBG = True
  if len(signalDictWeight) > 1:
    table = table +'<th>{0}</th>\n'.format("Total Signal")
    doTotalSig = True

  table = table +'</tr>\n<tr style="background-color:LightGrey">\n'

  for i in range(1,len(cutlist)+1):
    print cutlist[i-1],i-1
    # Setting shit up for the first row of values
    if i == 1:
      for hist in histlistWeight: oldcut.append(hist.GetBinContent(i))
      table = table +'<th>{0}</th>\n'.format(cutlist[i-1])
      for name,hist in sorted(bgDictWeight.iteritems(), key = lambda (k,v):(v,k)):
        table = table +u'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(hist.GetBinContent(i),hist.GetBinError(i))
        totalbg+= hist.GetBinContent(i)
        totalbgRaw += bgDictRaw[name].GetBinContent(i)
        totalbgerr += hist.GetBinError(i)**2
      for name,hist in sorted(signalDictWeight.iteritems(), key = lambda (k,v):(v,k)):
        table = table +u'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(hist.GetBinContent(i),hist.GetBinError(i))
        totalsig+= hist.GetBinContent(i)
        totalsigRaw += signalDictRaw[name].GetBinContent(i)
        totalsigerr += hist.GetBinError(i)**2
      print signalDictRaw
      for name,hist in sorted(dataDictWeight.iteritems(), key = lambda (k,v):(v,k)):
        table = table +u'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(hist.GetBinContent(i),hist.GetBinError(i))
      if (doTotalBG):
        table = table +'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(totalbg, sqrt(totalbgerr))
        totalbgCol.append(totalbg)
        totalbgRawCol.append(totalbgRaw)
        totalbgerrCol.append(totalbgerr)
        totalbg = 0
        totalbgRaw = 0
        totalbgerr = 0
      if (doTotalSig):
        table = table +'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(totalsig, sqrt(totalsigerr))
        totalsigCol.append(totalsig)
        totalsigRawCol.append(totalsigRaw)
        totalsigerrCol.append(totalsigerr)
        totalsig = 0
        totalsigRaw = 0
        totalsigerr = 0
      table = table +'</tr>\n'

    # Doing the rest of the rows
    else:
      skipCut = False
      for hist in histlistWeight: newcut.append(hist.GetBinContent(i))
      if newcut == oldcut: skipCut = True
      elif newcut != oldcut:
        del oldcut[:]
        del newcut[:]
        for hist in histlistWeight: oldcut.append(hist.GetBinContent(i))
      if not skipCut:
        if bgcol==False: table = table +'<tr style="background-color:LightGrey">\n'; bgcol = True
        else:      table = table +'<tr>\n';                                       bgcol = False
        table = table +'<th>{0}</th>\n'.format(cutlist[i-1])

        for name,hist in sorted(bgDictWeight.iteritems(), key = lambda (k,v):(v,k)):
          binomialerror = (bgDictWeight[name].GetBinContent(1)*sqrt(bgDictRaw[name].GetBinContent(i)*
            (1-bgDictRaw[name].GetBinContent(i)/bgDictRaw[name].GetBinContent(1)))/bgDictRaw[name].GetBinContent(1))
          table = table +u'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p>\n'.format(hist.GetBinContent(i),binomialerror,hist.GetBinContent(i)*100/hist.GetBinContent(1))
          totalbg+= hist.GetBinContent(i)
          totalbgRaw += bgDictRaw[name].GetBinContent(i)
          totalbgerr += binomialerror**2
        for name,hist in sorted(signalDictWeight.iteritems(), key = lambda (k,v):(v,k)):
          binomialerror = (signalDictWeight[name].GetBinContent(1)*sqrt(signalDictRaw[name].GetBinContent(i)*
            (1-signalDictRaw[name].GetBinContent(i)/signalDictRaw[name].GetBinContent(1)))/signalDictRaw[name].GetBinContent(1))
          table = table +u'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p>\n'.format(hist.GetBinContent(i),binomialerror, hist.GetBinContent(i)*300/hist.GetBinContent(1))
          totalsig+= hist.GetBinContent(i)
          totalsigRaw += signalDictRaw[name].GetBinContent(i)
          totalsigerr += binomialerror**2
        for name,hist in sorted(dataDictWeight.iteritems(), key = lambda (k,v):(v,k)):
          binomialerror = (dataDictWeight[name].GetBinContent(1)*sqrt(dataDictRaw[name].GetBinContent(i)*
            (1-dataDictRaw[name].GetBinContent(i)/dataDictRaw[name].GetBinContent(1)))/dataDictRaw[name].GetBinContent(1))
          table = table +u'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p>\n'.format(hist.GetBinContent(i),binomialerror, hist.GetBinContent(i)*100/hist.GetBinContent(1))
        if (doTotalBG):
          table = table +'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p>\n'.format(totalbg, sqrt(totalbgerr), totalbg*100/totalbgCol[0])
          totalbgCol.append(totalbg)
          totalbgRawCol.append(totalbgRaw)
          totalbgerrCol.append(totalbgerr)
          totalbg = 0
          totalbgRaw = 0
          totalbgerr = 0
        if (doTotalSig):
          table = table +'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p>\n'.format(totalsig, sqrt(totalsigerr), totalsig*100/totalsigCol[0])
          totalsigCol.append(totalsig)
          totalsigRawCol.append(totalsigRaw)
          totalsigerrCol.append(totalsigerr)
          totalsig = 0
          totalsigRaw = 0
          totalsigerr = 0
        table = table +'</tr>\n'

  # add the yields in the signal range of interest
  table = table +'<th>{0}</th>\n'.format('Signal Range: 115-180')
  for name,hist in sorted(bgDictMass.iteritems(), key = lambda (k,v):(v,k)):
    err = Double(0)
    area = hist.IntegralAndError(hist.FindBin(115),hist.FindBin(181),err)
    areaRaw = bgDictMassRaw[name].Integral(bgDictMassRaw[name].FindBin(115),bgDictMassRaw[name].FindBin(181))
    binomialerror = bgDictWeight[name].GetBinContent(1)*sqrt(areaRaw*(1-areaRaw/bgDictRaw[name].GetBinContent(1)))/bgDictRaw[name].GetBinContent(1)
    table = table +'<td>{0:.4f} &#177 {1:.4f} <p><i>{2:.4f}%</i></p></td>\n'.format(area,binomialerror, area*100/bgDictWeight[name].GetBinContent(1))
    totalbg+= area
    totalbgerr += binomialerror**2
  for name,hist in sorted(signalDictMass.iteritems(), key = lambda (k,v):(v,k)):
    err = Double(0)
    area = hist.IntegralAndError(hist.FindBin(115),hist.FindBin(181),err)
    areaRaw = signalDictMassRaw[name].Integral(signalDictMassRaw[name].FindBin(115),signalDictMassRaw[name].FindBin(181))
    binomialerror = signalDictWeight[name].GetBinContent(1)*sqrt(areaRaw*(1-areaRaw/signalDictRaw[name].GetBinContent(1)))/signalDictRaw[name].GetBinContent(1)
    table = table +'<td>{0:.4f} &#177 {1:.4f} <p><i>{2:.4f}%</i></p></td>\n'.format(area,binomialerror, area*300/signalDictWeight[name].GetBinContent(1))
    totalsig+= area
    totalsigerr += binomialerror**2
  for name,hist in sorted(dataDictMass.iteritems(), key = lambda (k,v):(v,k)):
    err = Double(0)
    area = hist.IntegralAndError(hist.FindBin(115),hist.FindBin(181),err)
    areaRaw = dataDictMassRaw[name].Integral(dataDictMassRaw[name].FindBin(115),dataDictMassRaw[name].FindBin(181))
    #binomialerror = dataDictWeight[name].GetBinContent(1)*sqrt(areaRaw*(1-areaRaw/dataDictRaw[name].GetBinContent(1)))/dataDictRaw[name].GetBinContent(1)
    table = table +'<td>{0:.4f} &#177 {1:.4f} <p><i>{2:.4f}%</i></p></td>\n'.format(area,err, area*100/dataDictWeight[name].GetBinContent(1))

  if (doTotalBG): table = table +u'<td>{0:.4f} &#177 {1:.4f}</td>\n'.format(totalbg,sqrt(totalbgerr))
  if (doTotalSig): table = table +u'<td>{0:.4f} &#177 {1:.4f}</td>\n'.format(totalsig,sqrt(totalsigerr))
  if bgcol==False: table = table +'<tr style="background-color:LightGrey">\n'; bgcol = True
  else:      table = table +'<tr>\n';                                       bgcol = False

  table = table +'</table>\n'
  return table



def makeAcceptance(histlistWeight, histlistRaw, histlistMass, cutlist,scale = {}):
  oldcut = []
  newcut = []
  initial = []
  initialRaw = []
  totalbg = 0
  totalbgRaw = 0
  totalbgerr=0
  totalsig = 0
  totalsigRaw = 0
  totalsigerr=0
  binomialerror = 0
  multiSig = False
  sigCount = 0
  lumiDict = {}
  dataPeriods = []
  histlistMassRaw = []
  bgcol=True
  table = '<table border = "10"\ncellpadding="5">\n<tr>\n'
  table = table +'<th></th>\n'
  #scaledClones = []
  scaleHists(histlistWeight,scale)
  for hist in histlistMass:
    histlistMassRaw.append(hist.Clone())
  scaleHists(histlistMass,scale)

  for hist in histlistWeight:
    table = table +'<th>{0}</th>\n'.format(dataName(hist.GetName()))
    if dataName(hist.GetName()).find("Signal")!=-1: sigCount+=1

  if sigCount > 1: multiSig = True

  table = table +'<th>{0}</th>\n'.format("Total BG")
  if multiSig: table = table +'<th>{0}</th>\n'.format("Total Signal")
  table = table +'</tr>\n<tr style="background-color:LightGrey">\n'

  print histlistWeight[0].GetNbinsX()
  for i in range(1,len(cutlist)+1):
    print cutlist[i-1],i-1
    if i == 1:
      for hist in histlistWeight: oldcut.append(hist.GetBinContent(i))
      print oldcut
      table = table +'<th>{0}</th>\n'.format(cutlist[i-1])
      for j,hist in enumerate(histlistWeight):
        table = table +u'<td>{0:.1f} &#177 {1:.1f}</td>\n'.format(hist.GetBinContent(i),hist.GetBinError(i))
        initial.append(hist.GetBinContent(i))
        initialRaw.append(histlistRaw[j].GetBinContent(i))
        if (hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find("Alt") == -1):
          totalbg+= hist.GetBinContent(i)
          totalbgRaw += histlistRaw[j].GetBinContent(i)
        if (multiSig and hist.GetName().find("Signal") != -1):
          totalsig+= hist.GetBinContent(i)
          totalsigRaw += histlistRaw[j].GetBinContent(i)


      print totalbg, totalbgRaw
      table = table +'<td>{0:.1f}</td>\n'.format(totalbg)
      initial.append(totalbg)
      initialRaw.append(totalbgRaw)
      totalbg = 0
      totalbgRaw= 0
      if multiSig:
        table = table +'<td>{0:.1f}</td>\n'.format(totalsig)
        initial.append(totalsig)
        initialRaw.append(totalsigRaw)
        totalsig = 0
        totalsigRaw= 0
      table = table +'</tr>\n'
    else:
      for hist in histlistWeight: newcut.append(hist.GetBinContent(i))
      print oldcut, newcut
      if newcut == oldcut:
        del oldcut[:]
        del newcut[:]
        for hist in histlistWeight: oldcut.append(hist.GetBinContent(i))
      elif newcut != oldcut:
        del oldcut[:]
        del newcut[:]
        for hist in histlistWeight: oldcut.append(hist.GetBinContent(i))
        if bgcol==False: table = table +'<tr style="background-color:LightGrey">\n'; bgcol = True
        else:      table = table +'<tr>\n';                                       bgcol = False
        table = table +'<th>{0}</th>\n'.format(cutlist[i-1])
        for j,hist in enumerate(histlistWeight):
          print hist, initialRaw
          binomialerror = initial[j]*sqrt(histlistRaw[j].GetBinContent(i)*(1-histlistRaw[j].GetBinContent(i)/initialRaw[j]))/initialRaw[j]
          print 'binomial error for '+dataName(hist.GetName())+':',binomialerror
          if dataName(hist.GetName()).find('DATA') != -1: table = table +'<td>{0:.1f} &#177 {1:.0f} <p><i>{2:.3f}%</i></p></td>\n'.format(hist.GetBinContent(i),hist.GetBinError(i), hist.GetBinContent(i)*100/initial[j])
          elif dataName(hist.GetName()).find('Signal') != -1: table = table +'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p></td>\n'.format(hist.GetBinContent(i),binomialerror, hist.GetBinContent(i)*100/initial[j])
          else: table = table +'<td>{0:.1f} &#177 {1:.1f} <p><i>{2:.3f}%</i></p></td>\n'.format(hist.GetBinContent(i),binomialerror, hist.GetBinContent(i)*100/initial[j])
          if (hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find("Alt") == -1):
            totalbg+= hist.GetBinContent(i)
            totalbgerr += binomialerror**2

          if (multiSig and hist.GetName().find("Signal") != -1):
            totalsig+= hist.GetBinContent(i)
            totalsigerr += binomialerror**2

        table = table +'<td>{0:.0f} &#177 {1:.0f} <p><i>{2:.3f}%</i></p></td>\n'.format(totalbg, totalbgerr**0.5, totalbg*100/initial[-1])
        totalbg=0
        totalbgerr=0

        if (multiSig):
          table = table +'<td>{0:.2f} &#177 {1:.2f} <p><i>{2:.3f}%</i></p></td>\n'.format(totalsig, totalsigerr**0.5, totalsig*100/initial[-1])
          totalsig=0
          totalsigerr=0

        table = table +'</tr>\n'

  # add the yields in the signal range of interest
  table = table +'<th>{0}</th>\n'.format('Signal Range: 124-128')
  for j,hist in enumerate(histlistMass):
    err = Double(0)
    area = hist.IntegralAndError(hist.FindBin(124),hist.FindBin(127),err)
    areaRaw = histlistMassRaw[j].Integral(histlistMassRaw[j].FindBin(124),histlistMassRaw[j].FindBin(127))
    binomialerror = initial[j]*sqrt(areaRaw*(1-areaRaw/initialRaw[j]))/initialRaw[j]
    print 'binomial error for '+dataName(hist.GetName())+':',binomialerror
    #time.sleep(1)
    if dataName(hist.GetName()).find('DATA') != -1: table = table +'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(area,err)
    else: table = table +'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(area,binomialerror)
    if (hist.GetName().find("DATA") == -1 and hist.GetName().find("Signal") == -1 and hist.GetName().find("Extra") == -1 and hist.GetName().find("Alt") == -1):
      totalbg+= area
      totalbgerr += binomialerror**2

    if (multiSig and hist.GetName().find("Signal") != -1):
      totalsig+= area
      totalsigerr += binomialerror**2

  print totalbg, totalbgRaw
  table = table +u'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(totalbg,sqrt(totalbgerr))
  if multiSig: table = table +u'<td>{0:.2f} &#177 {1:.2f}</td>\n'.format(totalsig,sqrt(totalsigerr))
  if bgcol==False: table = table +'<tr style="background-color:LightGrey">\n'; bgcol = True
  else:      table = table +'<tr>\n';                                       bgcol = False

  table = table +'</table>\n'
  return table


def makeWebPage(path):
  if not os.path.isfile(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/writeIndex.csh'):
    shutil.copy(os.getenv("HOME")+'/nobackup/TempWeb/writeIndex.csh',os.getenv("HOME")+'/nobackup/TempWeb/'+path)
  os.system('cd ~/nobackup/TempWeb/'+path+'; ./writeIndex.csh')

def initWebPage(path):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  if not os.path.isdir(os.getenv("HOME")+'/nobackup/TempWeb/'+year+'/'+month+'/'+day+'/'+path):
    os.makedirs(os.getenv("HOME")+'/nobackup/TempWeb/'+year+'/'+month+'/'+day+'/'+path)

def makeWebYeild(path,acc):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  if not os.path.isdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path):
    os.mkdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path)
  title = "<html><head><title>%s</title></head><body>\n"
  upPath = 'http://home.fnal.gov/~bpollack/'+year+'/'+month+'/'+day
  parentDir = '<a href="'+upPath+'">Parent Directory</a>\n'
  header = "<h1>%s</h1>\n"
  index = codecs.open(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/index.html',encoding='utf-8',mode='w')
  index.write(acc)
  index.close()


def makeWebPageBetter(path,acc=None):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  if not os.path.isdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path):
    os.mkdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path)
  title = "<html><head><title>%s</title></head><body>\n"
  upPath = 'http://home.fnal.gov/~bpollack/'+year+'/'+month+'/'+day
  parentDir = '<a href="'+upPath+'">Parent Directory</a>\n'
  header = "<h1>%s</h1>\n"
  image  = "<img src={0} width=47%>\n"
  os.system('cd ~/nobackup/TempWeb/'+path+';ls -1 *png > temp.txt')
  f = open(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/temp.txt', 'r')
  dic = defaultdict(list)
  for line in f:
    value = line.strip()
    key = (value.rstrip(string.letters+'.0123456789')).rstrip('_')
    dic[key].append(value)
  f.close()
  index = open(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/index.html','w')
  index.write(title % path)
  #index.write(parentDir)
  if acc != None: index.write(acc)
  for key in dic:
    #index.write(header % key)
    for val in dic[key]:
      index.write(image.format(val))
    #index.write('<hr />\n')

  index.write("</body></html>")
  index.close()

def makeWebPageBest(path, selectionInfo=[''], title = 'Test'):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  if not os.path.isdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path):
    os.mkdir(os.getenv("HOME")+'/nobackup/TempWeb/'+path)
  #title = "<html><head><title>%s</title></head><body>\n"
  upPath = 'http://home.fnal.gov/~bpollack/'+year+'/'+month+'/'+day
  parentDir = '<a href="'+upPath+'">Parent Directory</a>\n'
  header = "<h1>%s</h1>\n"
  image  = "<img src=%s >\n"
  os.system('cd ~/nobackup/TempWeb/'+path+';ls -d */ -1 > temp.txt')
  f = open(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/temp.txt', 'r')
  dic = dict()
  for line in f:
    #key = (line.strip()).lstrip('\/')
    key = (line.strip()).strip('/')
    value = key
    dic[key] = value
  f.close()
  index = open(os.getenv("HOME")+'/nobackup/TempWeb/'+path+'/index.html','w')
  #index.write(title % path)
  #index.write(parentDir)
  #if acc != None: index.write(acc)
  #for key in dic:
    #index.write(header % key)
    #for val in dic[key]:
     # index.write(image % val)
    #index.write('<hr />\n')
  #index.write("</body></html>")
########## OLD NAV STUFF #########

    #nav li { float: left; width: 10%; text-align:center; padding: 0.5%; color: #fff; font-size: 150%;}
    #nav ul { list-style: none; overflow: hidden; padding: 0; margin: 0; background-color: #424242; }
    #a { color: #65ac1c; }
    #nav li a { color: #fff; }
##################################

  index.write('''
<!DOCTYPE html>
<html>
  <head>
  <style>
    body {position:absolute; top:0px; bottom:0px; left:0px; right:0px; margin: 0 auto; width: 99%; height: 100%; background-color: #6F8C8F;
      color: #EAF1FB;}
    section { display: block; }
    section#plots { width: 80%; height: 100%; float: left; padding: 0.5%; background-color: #FFFFFF; }
    article > header { text-decoration: underline; margin-bottom: 0.5%; }
    aside { float: left; width: 15%; padding: 0.5%; font-size: 110%; }
    footer { overflow: hidden; clear: both; text-align: center; padding: 1.5%; }
    nav  {
        display: block;
        width: 100%;
        overflow: hidden;
    }

    nav ul {
            margin: 1% auto;
            padding: 1%;
            float: left;
            list-style: none;
            background: #444;
            background: rgba(0,0,0,.2);
            -moz-border-radius: .5em;
            -webkit-border-radius: .5em;
            border-radius: .5em;
            -moz-box-shadow: 0 1px 0 rgba(255,255,255,.2), 0 2px 1px rgba(0,0,0,.8) inset;
            -webkit-box-shadow: 0 1px 0 rgba(255,255,255,.2), 0 2px 1px rgba(0,0,0,.8) inset;
            box-shadow: 0 1px 0 rgba(255,255,255,.2), 0 2px 1px rgba(0,0,0,.8) inset;
    }

    nav li {
            float: left;

    }

    nav a {
            float:left;
            padding: .8em 1.5em;
            text-decoration: none;
            color: #555;
            text-shadow: 0 1px 0 rgba(255,255,255,.5);
            font: bold 1.1em/1 'trebuchet MS', Arial, Helvetica;
            letter-spacing: 1px;
            text-transform: none;
            border-width: 1px;
            border-style: solid;
            border-color: #fff #ccc #999 #eee;
            background: #c1c1c1;
            background: -moz-linear-gradient(#f5f5f5, #c1c1c1);
            background: -webkit-gradient(linear, left top, left bottom, from(#f5f5f5), to(#c1c1c1));
            background: -webkit-linear-gradient(#f5f5f5, #c1c1c1);
            background: -o-linear-gradient(#f5f5f5, #c1c1c1);
            background: -ms-linear-gradient(#f5f5f5, #c1c1c1);
            background: linear-gradient(#f5f5f5, #c1c1c1);
     }

    nav a:hover, nav a:focus {
            outline: 0;
            color: #fff;
            text-shadow: 0 1px 0 rgba(0,0,0,.2);
            background: #fac754;
            background: -moz-linear-gradient(#fac754, #f8ac00);
            background: -webkit-gradient(linear, left top, left bottom, from(#fac754), to(#f8ac00));
            background: -webkit-linear-gradient(#fac754, #f8ac00);
            background: -o-linear-gradient(#fac754, #f8ac00);
            background: -ms-linear-gradient(#fac754, #f8ac00);
            background: linear-gradient(#fac754, #f8ac00);
    }

    nav a:active {
            -moz-box-shadow: 0 0 2px 2px rgba(0,0,0,.3) inset;
            -webkit-box-shadow: 0 0 2px 2px rgba(0,0,0,.3) inset;
            box-shadow: 0 0 2px 2px rgba(0,0,0,.3) inset;
    }

    nav li:first-child a {
            border-left: 0;
            -moz-border-radius: 4px 0 0 4px;
            -webkit-border-radius: 4px 0 0 4px;
            border-radius: 4px 0 0 4px;
    }

    nav li:last-child a {
            border-right: 0;
            -moz-border-radius: 0 4px 4px 0;
            -webkit-border-radius: 0 4px 4px 0;
            border-radius: 0 4px 4px 0;
    }
  </style>
  </head>
  <body style="font: 95% Lucida Grande, Trebuchet MS">
    <header>
      <hgroup>
      ''')
  index.write('''
        <h1>{0}</h1>
        '''.format(title))
  index.write('''
        <h2>{0} {1}, {2}</h2>
        '''.format(today.strftime('%B'),today.strftime('%d'),today.strftime('%Y')))
  index.write('''
      </hgroup>
    </header>
    <nav>
      <ul>
      ''')
  index.write('''        <li><a href="{0}" target="iframe_a">{1}</a></li>'''.format(dic['Yields'],'Yields'))
  for link in dic:
    #  print link,dic[link]
    if link != 'Yields':
      index.write('''        <li><a href="{0}" target="iframe_a">{1}</a></li>'''.format(dic[link],link))
  index.write('''
      </ul>
    </nav>
    <section id="plots">
    <iframe src="Yields" name="iframe_a" frameborder="0" width = "100%" height = "100%"></iframe>
    </section>
    <aside>
      <h2>Selection/Cut Info</h2>
      <ul>
      ''')
  for item in selectionInfo:
    index.write('''       <li>{0}</li>'''.format(item))
  index.write('''
      </ul>
    </aside>

    <footer>
      Copyright BPollack Industries, MMXI
    </footer>
  </body>
</html>
''')
  index.close()



def MeanWidthVTX(hist = [],fitType = None,label = 'placeholder', metType = 'TrackMET', xMin = None, xMax = None,tag = 'MC',maintitle = ''):
  corrections = None
  if tag == 'Data':
    h1_Means = TH1F('h1_Means','',9,0.5,9.5)
    h1_Widths = TH1F('h1_Widths','',9,0.5,9.5)
  if tag == 'MC':
    h1_Means = TH1F('h1_Means','',12,0.5,12.5)
    h1_Widths = TH1F('h1_Widths','',12,0.5,12.5)
  if (fitType == None or fitType == 'gaus'):
    for i,x in enumerate(hist):
      if x != None:
        name = hist[i].GetName()
        name.lstrip("h1_")
        name.rstrip("muPromptV4DYToMuMu_")
        fit = makeFit(x,fitType,xmin = xMin, xmax = xMax)
        h1_Means.Fill(i+1,fit.GetParameter(1))
        h1_Means.SetBinError(i+1,fit.GetParError(1))
        h1_Widths.Fill(i+1,fit.GetParameter(2))
        h1_Widths.SetBinError(i+1,fit.GetParError(2))
        ppps(metType+'/VTXs/'+label.lstrip(metType+"_"),label+'_VTX_'+str(i+1).zfill(2),x,title = label+" VTX "+str(i+1).zfill(2)+tag,xmin = 0, xmax = 150,lstr = 'h1_', rstr = 'muPromptV4DYToMuMuJet012')
    makeWebPage(metType+'/VTXs/'+label.lstrip(metType+"_"))
    fit = makeFit(h1_Means,'pol1')
    fit = makeFit(h1_Widths,'pol1')
    ppps(metType,label+'_Means',h1_Means,title = label+', Means')
    ppps(metType,label+'_Widths',h1_Widths,title = label+', Widths')

  if fitType == 'novo':
    for i,x in enumerate(hist):
      if x != None and i<12:
        name = hist[i].GetName()
        name.lstrip("h1_")
        name.rstrip("muPromptV4DYToMuMu_")
        fit = makeFit(x,fitType,xmin = xMin, xmax = xMax)
        h1_Means.Fill(i+1,fit.GetParameter(1))
        h1_Means.SetBinError(i+1,fit.GetParError(1))
        h1_Widths.Fill(i+1,fit.GetParameter(2))
        h1_Widths.SetBinError(i+1,fit.GetParError(2))
        ppps(metType+'/VTXs/'+label.lstrip(metType+"_"),label+'_VTX_'+str(i+1).zfill(2),x,title = label+" VTX "+str(i+1).zfill(2)+tag,xmin = 0, xmax = 150,lstr = 'h1_', rstr = 'muPromptV4DYToMuMuJet012')
    makeWebPage(metType+'/VTXs/'+label.lstrip(metType+"_"))
    h1_WidthsClone = h1_Widths.Clone()
### shitty toggle for cor values
    #gStyle.SetOptStat(1)
    #h1_Widths.SetStats(1)
### tooooooggggggggle
    fit = makeFit(h1_Widths,'pol1')
    tempCan = TCanvas('canvas','canvas',800,600)
    h1_Widths.Draw()
    gPad.Update()
    st1 = h1_Widths.FindObject("stats")
    st2 = st1.Clone()
    h1_Widths.GetListOfFunctions().Remove(st1)
    #h1_Widths.SetStats(1)
    st2.SetName("stPol1")
    st2.SetX1NDC(0.2) #new x start position
    st2.SetX2NDC(0.4) #new x end position
    st2.SetY1NDC(0.8) #new y start position
    st2.SetY2NDC(0.9) #new y end position
    st2.Draw()
    gPad.Update()
    h1_Widths.GetListOfFunctions().Add(st2)
    gStyle.SetOptStat(1)
    fitCor1 = makeFit(h1_Means,'peak')
    fitCor2 = makeFit(h1_Widths,'width')
    fit = makeFit(h1_Widths,'pol1Again')
    st3 = h1_Widths.FindObject("stats")
    st3.SetName("stWidth")
    st3.SetX1NDC(0.71) #new x start position
    st3.SetX2NDC(0.91) #new x end position
    st3.SetY1NDC(0.8) #new y start position
    st3.SetY2NDC(1.0) #new y end position
    gPad.Update()
    h1_Widths.GetListOfFunctions().Add(st3)
    gPad.Update()

    ppps(metType,label+'_Peaks',h1_Means,title = maintitle+' Peaks',xtitle = 'nVtx', ytitle = 'GeV', legOn=False)
    ppps(metType,label+'_Widths',h1_Widths,title = maintitle+' Widths',xtitle = 'nVtx',ytitle = 'GeV', legOn=False)
    corrections = [fitCor1.GetParameter(0), fitCor1.GetParameter(1), fitCor2.GetParameter(0), fitCor2.GetParameter(1)]
    #corrections =[h1_Means,h1_Widths]

  if fitType == 'novoCor':
    for i,x in enumerate(hist):
      if x != None and i < 12:
        name = hist[i].GetName()
        name.lstrip("h1_")
        name.rstrip("muPromptV4DYToMuMu_")
        fit = makeFit(x,'novo',xmin = xMin, xmax = xMax)
        h1_Means.Fill(i+1,fit.GetParameter(1))
        h1_Means.SetBinError(i+1,fit.GetParError(1))
        h1_Widths.Fill(i+1,fit.GetParameter(2))
        h1_Widths.SetBinError(i+1,fit.GetParError(2))
        ppps(metType+'/VTXs/'+label.lstrip(metType+"_"),'Cor_'+label+'_VTX_'+str(i+1).zfill(2),x,title = 'PUCor '+label+" VTX "+str(i+1).zfill(2)+tag,xmin = 0, xmax = 150,lstr = 'h1_', rstr = 'muPromptV4DYToMuMuJet012')
    makeWebPage(metType+'/VTXs/'+label.lstrip(metType+"_"))
    fit = makeFit(h1_Widths,'pol1')
    fitCor1 = makeFit(h1_Means,'pol1')

    ppps(metType,label+'_Peaks',h1_Means,title = 'Corrected '+maintitle+' Peaks',xtitle = 'nVtx', ytitle = 'GeV', legOn=False)
    ppps(metType,label+'_Widths',h1_Widths,title = 'Corrected '+maintitle+' Widths',xtitle = 'nVtx',ytitle = 'GeV', legOn=False)
    #corrections = [fitCor1.GetParameter(0), fitCor1.GetParameter(1), fitCor2.GetParameter(0), fitCor2.GetParameter(1)]
    #corrections =[h1_Means,h1_Widths]

  if fitType == 'expo':
    for i,x in enumerate(hist):
      if x != None:
        name = hist[i].GetName()
        name.lstrip("h1_")
        name.rstrip("muPromptV4DYToMuMu_")
        fit = makeFit(x,fitType,xmin = xMin, xmax = xMax)
        h1_Means.Fill(i+1,fit.GetParameter(0))
        h1_Means.SetBinError(i+1,fit.GetParError(0))
        h1_Widths.Fill(i+1,fit.GetParameter(1))
        h1_Widths.SetBinError(i+1,fit.GetParError(1))
        ppps(metType+'/VTXs/'+label.lstrip(metType+"_"),label+'_VTX_'+str(i+1).zfill(2),x,title = label+" VTX "+str(i+1).zfill(2)+tag,xmin = 0, xmax = 80,lstr = 'h1_', rstr = 'muPromptV4DYToMuMuJet012')
    makeWebPage(metType+'/VTXs/'+label.lstrip(metType+"_"))
    fit = makeFit(h1_Means,'pol1')
    fit = makeFit(h1_Widths,'pol1')
    ppps(metType,label+'_Consts',h1_Means,title = label+', Constants')
    ppps(metType,label+'_Slopes',h1_Widths,title = label+', Slopes')

  if fitType == 'novoNW':
    for i,x in enumerate(hist):
      if x != None:
        name = hist[i].GetName()
        name.lstrip("h1_")
        name.rstrip("muPromptV4DYToMuMu_")
        fit = makeFit(x,'novo',xmin = xMin, xmax = xMax)
        h1_Means.Fill(i+1,fit.GetParameter(1))
        h1_Means.SetBinError(i+1,fit.GetParError(1))
        h1_Widths.Fill(i+1,fit.GetParameter(2))
        h1_Widths.SetBinError(i+1,fit.GetParError(2))
    fit = makeFit(h1_Widths,'pol1')
    fitCor1 = makeFit(h1_Means,'pol1')

    corrections =[h1_Means,h1_Widths]
  return corrections

def PUCompare(hist,title,ytitle,tag,legnames,saveas,path = 'HWTalk',table = dict()):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  canvas = TCanvas('canvas','canvas',800,600)
  if (title.find('ZZ') != -1) or (title.find('H200') != -1):
    h2_Temp = TH2F("h2_Temp","",9,0.5,9.5,41,5,45)
  else:
    h2_Temp = TH2F("h2_Temp","",9,0.5,9.5,16,0,15)
  #elif tag == 'Data': h2_Temp = TH2F("h2_Temp","",9,0.5,9.5,16,3,19)
  #else: h2_Temp = TH2F("h2_Temp","",12,0.5,12.5,16,3,19)
  h2_Temp.Draw()
  h2_Temp.SetStats(0)
  pt = makePretty(h2_Temp,title,'nVtx',ytitle)
  legd = None
  gStyle.SetOptFit(1)
  dataName = saveas.lstrip(string.letters+'0123456789')
  dataName = dataName.lstrip('_')
  statBox = []
  for i,x in enumerate(hist):
    if i <10 and x.GetFunction(x.GetName()+'_fit')!= None :
      if legnames[i] not in table.keys():
        table[legnames[i]] = defaultdict(list)
      x.Draw("sames")
      x.SetLineColor(ColorDict[i])
      x.SetMarkerColor(ColorDict[i])
      x.SetMarkerStyle(MarkerList[i])
      x.SetMarkerSize(2)
      x.GetFunction(x.GetName()+'_fit').SetLineColor(ColorDict[i])
      x.GetFunction(x.GetName()+'_fit').SetLineStyle(9)
      x.GetFunction(x.GetName()+'_fit').SetLineWidth(2)
      ch2 = x.GetFunction(x.GetName()+'_fit').GetChisquare()
      ndof = x.GetFunction(x.GetName()+'_fit').GetNDF()
      p0 = x.GetFunction(x.GetName()+'_fit').GetParameter(0)
      p1 = x.GetFunction(x.GetName()+'_fit').GetParameter(1)
      p0err = x.GetFunction(x.GetName()+'_fit').GetParError(0)
      p1err = x.GetFunction(x.GetName()+'_fit').GetParError(1)
      table[legnames[i]][dataName].append(p0)
      table[legnames[i]][dataName].append(p0err)
      table[legnames[i]][dataName].append(p1)
      table[legnames[i]][dataName].append(p1err)
      gPad.Update()
      statBox.append(x.FindObject("stats"))
      x.GetListOfFunctions().Remove(statBox[i])
      x.SetStats(0)
      lineList = statBox[i].GetListOfLines()
      txt = TLatex(0,0,legnames[i])
      lineList.AddFirst(txt)
      #time.sleep(5)
      statBox[i].SetX1NDC(0.80) #new x start position
      statBox[i].SetX2NDC(1.0) #new x end position
      statBox[i].SetY1NDC(0.8-i*0.14) #new y start position
      statBox[i].SetY2NDC(0.93-i*0.14) #new y end position
      statBox[i].Paint()
      gPad.Update()
      legd = makeLegend(x,leg = legd,cname = legnames[i])
      legd.Draw()
      canvas.Update()
      gPad.Update()
  #for box in statBox:
    #    box.Draw()
  #table = makeLatexTables(len(hist),legnames,ch2,ndof,p0,p0err,p1,p1err,title)
  canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+saveas+'.png')
  #return table

def makeLatexTables(numOfFits,fitNames,chi2,ndof,p0,p0err,p1,p1err,title):
  table = '\\begin{tabular}{|l|l|'+'|l|l|'*(numOfFits-1)+'}\n \hline\n'
  table = table + ' \multicolumn{'+str(2*numOfFits)+'}{|c|}{%s} \\\\ \n \hline\n' % title
  table = table + ' \multicolumn{2}{|c||}{%s}' % fitNames[0]
  for i in range(1,numOfFits):
    if i < numOfFits-1:
      table = table + ' & \multicolumn{2}{|c||}{%s}' % fitNames[i]
    else:
      table = table + ' & \multicolumn{2}{|c|}{%s}' % fitNames[i]
  table = table + '\\\\ \n \hline\n'
  table = table + ' $\\frac{\chi^2}{ndof}$ & %.3f/%d' % (chi2[0],ndof[0])
  for i in range(1,numOfFits):
    table = table + ' & $\\frac{\chi^2}{ndof}$ & %.3f/%d' % (chi2[i],ndof[i])
  table = table + ' \\\\ \n'
  table = table + ' y$_{0}$ & %.3f $\pm$ %.3f' % (p0[0],p0err[0])
  for i in range(1,numOfFits):
    table = table + ' & y$_{0}$ & %.3f $\pm$ %.3f' % (p0[i],p0err[i])
  table = table + ' \\\\ \n'
  table = table + ' slope & %.3f $\pm$ %.3f' % (p1[0],p1err[0])
  for i in range(1,numOfFits):
    table = table + ' & slope & %.3f $\pm$ %.3f' % (p1[i],p1err[i])
  table = table + ' \\\\ \n'
  table = table + ' \hline\n'
  table = table + '\end{tabular}'
  return table


def makeLatexTablesLarge(megaDict):
  finalList=[]
  metList = ['Pf MET','Type 1 Corrected MET','Projected (lep) MET','Projected (Z) MET',
    'Track MET (all Jets)','Track MET (0 Jets)','Track MET (1 Jet)','Track MET (2+ Jets)',
    'Reduced MET v1','Reduced MET v2']

  dataList = ['muData','DYJetsToLL','ZZJetsTo2L2Nu','H200']

  table = '\\begin{tabular}{|p{2cm}|l|l|l|l||l|l|l|l|}\n'
  table = table + '\multicolumn{9}{c}{Peak Values (Intercept)} \\\\ \n \hline\n'
  table = table + 'Intercept & \multicolumn{4}{|c||}{No Correction} & \multicolumn{4}{|c|}{PU Corrected} \\\\ \n \hline\n'
  for key in dataList:
    table = table + ' & %s' % key
  for key in dataList:
    table = table + ' & %s' % key
  table = table + ' \\\\ \n \hline\n'

  for met in metList:
    table = table + met
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][0],megaDict[met][data][1])
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][4],megaDict[met][data][5])
    table = table + ' \\\\ \n \hline\n'
  table = table + '\end{tabular}'

  finalList.append(table) ; table = ''

  table = '\\begin{tabular}{|p{2cm}|l|l|l|l||l|l|l|l|}\n'
  table = table + '\multicolumn{9}{c}{Peak Values (Slope)} \\\\ \n \hline\n'
  table = table + 'Slope & \multicolumn{4}{|c||}{No Correction} & \multicolumn{4}{|c|}{PU Corrected} \\\\ \n \hline\n'
  for key in dataList:
    table = table + ' & %s' % key
  for key in dataList:
    table = table + ' & %s' % key
  table = table + ' \\\\ \n \hline\n'

  for met in metList:
    table = table + met
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][2],megaDict[met][data][3])
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][6],megaDict[met][data][7])
    table = table + ' \\\\ \n \hline\n'
  table = table + '\end{tabular}'

  finalList.append(table) ; table = ''

  table = '\\begin{tabular}{|p{2cm}|l|l|l|l||l|l|l|l|}\n'
  table = table + '\multicolumn{9}{c}{Width Values (Intercept)} \\\\ \n \hline\n'
  table = table + 'Intercept & \multicolumn{4}{|c||}{No Correction} & \multicolumn{4}{|c|}{PU Corrected} \\\\ \n \hline\n'
  for key in dataList:
    table = table + ' & %s' % key
  for key in dataList:
    table = table + ' & %s' % key
  table = table + ' \\\\ \n \hline\n'

  for met in metList:
    table = table + met
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][8],megaDict[met][data][9])
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][12],megaDict[met][data][13])
    table = table + ' \\\\ \n \hline\n'
  table = table + '\end{tabular}'

  finalList.append(table) ; table = ''

  table = '\\begin{tabular}{|p{2cm}|l|l|l|l||l|l|l|l|}\n'
  table = table + '\multicolumn{9}{c}{Width Values (Slope)} \\\\ \n \hline\n'
  table = table + 'Slope & \multicolumn{4}{|c||}{No Correction} & \multicolumn{4}{|c|}{PU Corrected} \\\\ \n \hline\n'
  for key in dataList:
    table = table + ' & %s' % key
  for key in dataList:
    table = table + ' & %s' % key
  table = table + ' \\\\ \n \hline\n'

  for met in metList:
    table = table + met
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][10],megaDict[met][data][11])
    for data in dataList:
      table = table + ' & %.3f $\pm$ %.3f' % (megaDict[met][data][14],megaDict[met][data][15])
    table = table + ' \\\\ \n \hline\n'
  table = table + '\end{tabular}'

  finalList.append(table) ; table = ''

  return finalList

def dataName(histName):
  dataName = histName.lstrip(string.letters+'0123456789')
  dataName = dataName.lstrip('_')
  dataName = dataName.lstrip(string.letters+'0123456789')
  dataName = dataName.lstrip('_')
  return dataName

def plotName(histName):
  plotName = histName.strip(string.letters+'0123456789')
  plotName = plotName.strip('_')
  return plotName

def doublePlotRatio(path,filename,histList,logToggle = False,axisLabels = None):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path

  canvas = TCanvas('canvas','canvas',900,700)
  pad1 = TPad('pad1', '', 0.02, 0.34, 0.90, 0.98, 0)
  pad1.SetBottomMargin(0.07)
  pad1.SetGrid(1,1)
  if logToggle: pad1.SetLogy(kTRUE)
  pad2 = TPad('pad2', '', 0.02, 0.02, 0.90, 0.35, 0)
  pad2.SetTopMargin(0.)
  pad2.SetBottomMargin(0.2)
  pad2.SetGrid(1,1)
  pad1.Draw()
  pad2.Draw()
  maxbinval = max(histList[0].GetBinContent(histList[0].GetMaximumBin()),histList[1].GetBinContent(histList[1].GetMaximumBin()))
  name = plotName(histList[0].GetName())
  leg = makeLegend(histList,x1 = 0.2,x2 = 0.43)
  pad1.cd()
  h2_Temp = TH2F("h2_Temp",histList[0].GetTitle()+';;'+histList[0].GetYaxis().GetTitle(),histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1),50,0,maxbinval+maxbinval*(0.1))
  h2_Temp.SetStats(0)
  if axisLabels != None:
    for i in range(len(axisLabels)): h2_Temp.GetXaxis().SetBinLabel(i+1,axisLabels[i])
  h2_Temp.Draw()
  #print histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1)
  histList[0].SetMarkerStyle(20)
  histList[0].SetOption('EP')
  histList[0].SetStats(0)
  histList[0].Draw('same EP')
  histList[1].SetMarkerStyle(21)
  histList[1].SetMarkerColor(kRed+1)
  histList[1].SetLineColor(kRed+1)
  histList[1].SetOption('EP')
  histList[1].SetStats(0)
  histList[1].Draw('same EP')
  leg.Draw()
  pt = TPaveText(0.2,0.930,0.85,1,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(11)
  #pt.SetTextSize(2)
  text = pt.AddText("")
  pt.Draw()



  pad2.cd()
  histMC = histList[0]
  for hist in histList[1:-2]:
    histMC = histMC + hist
  nBins = histMC.GetNbinsX()
  xAxisName = histMC.GetXaxis().GetTitle()
  hRatio = TH1D('hRatio', ';'+xAxisName+';LumiCalc/ZCount', nBins, histMC.GetBinLowEdge(1), histMC.GetBinLowEdge(nBins + 1))
  hRatio.Divide(histList[1],histList[0])

  hRatio.SetStats(0)
  hRatio.SetLineWidth(2);
  hRatio.SetLineColor(1);
  hRatio.SetMarkerStyle(20);
  hRatio.GetYaxis().SetLabelSize(0.08);
  hRatio.GetYaxis().SetTitleSize(0.09);
  hRatio.GetYaxis().SetTitleOffset(0.40);
  hRatio.GetYaxis().SetNdivisions(5);
  hRatio.GetYaxis().CenterTitle();
  hRatio.GetXaxis().SetLabelSize(0.08);
  hRatio.GetXaxis().SetTitleSize(0.09);
  hRatio.GetXaxis().SetTitleOffset(0.90);
  hRatio.SetMaximum(min(hRatio.GetBinContent(hRatio.GetMaximumBin())*1.2,4));
  hRatio.SetMinimum(0);
  fit = TF1("ratioFit","pol0")
  fit.SetLineWidth(3)
  fit.SetLineStyle(10)
  hRatio.Fit("ratioFit")
  if axisLabels != None:
    for i in range(len(axisLabels)): hRatio.GetXaxis().SetBinLabel(i+1,axisLabels[i])
  hRatio.Draw()

  canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'.png')

  canvas2 = TCanvas('canvas','canvas',900,700)
  canvas2.cd()
  canvas2.SetGrid(1,1)
  h2_Temp.SetStats(0)
  h2_Temp.GetXaxis().SetTitle("Run Period")
  if axisLabels != None:
    for i in range(len(axisLabels)): h2_Temp.GetXaxis().SetBinLabel(i+1,axisLabels[i])
  h2_Temp.Draw()
  #print histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1)
  histList[0].SetMarkerStyle(20)
  histList[0].SetOption('EP')
  histList[0].SetStats(0)
  histList[0].Draw('same EP')
  histList[1].SetMarkerStyle(21)
  histList[1].SetMarkerColor(kRed+1)
  histList[1].SetLineColor(kRed+1)
  histList[1].SetOption('EP')
  histList[1].SetStats(0)
  histList[1].Draw('same EP')
  leg.Draw()
  pt = TPaveText(0.2,0.930,0.85,1,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(21)
  #pt.SetTextSize(2)
  text = pt.AddText("")
  pt.Draw()

  canvas2.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'_1.png')

  canvas3 = TCanvas('canvas','canvas',900,700)
  canvas3.cd()
  canvas3.SetGrid(1,1)
  histMC = histList[0]
  for hist in histList[1:-2]:
    histMC = histMC + hist
  nBins = histMC.GetNbinsX()
  xAxisName = histMC.GetXaxis().GetTitle()
  hRatio = TH1D('hRatio', ';'+xAxisName+';LumiCalc/ZCountLumi', nBins, histMC.GetBinLowEdge(1), histMC.GetBinLowEdge(nBins + 1))
  hRatio.Divide(histList[1],histList[0])

  hRatio.SetStats(1)
  hRatio.SetLineWidth(2);
  hRatio.SetLineColor(1);
  hRatio.SetMarkerStyle(20);
  hRatio.GetYaxis().SetNdivisions(5);
  hRatio.GetYaxis().CenterTitle();
  hRatio.SetMaximum(1.3);
  hRatio.SetMinimum(0.7);
  fit = TF1("ratioFit","pol0")
  fit.SetLineWidth(3)
  fit.SetLineStyle(10)
  hRatio.Fit("ratioFit")
  if axisLabels != None:
    for i in range(len(axisLabels)): hRatio.GetXaxis().SetBinLabel(i+1,axisLabels[i])
  hRatio.SetTitle("(LumiCalc / Calculated Lumi)")
  hRatio.Draw()
  pt = TPaveText(0.2,0.930,0.85,1,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(11)
  #pt.SetTextSize(2)
  text = pt.AddText("")
  pt.Draw()
  line = TLine(histMC.GetBinLowEdge(1),1.07,histMC.GetBinLowEdge(nBins + 1),1.07)
  line.SetLineWidth(3)
  line.SetLineColor(kGreen+1)
  line.SetLineStyle(9)
  line.Draw()

  canvas3.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'_2.png')

def stacksAndRatio(path,filename,histList,logToggle = False,cname = None):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path
  totalDict = {}
  dataDict = {}
  bgDict = {}
  signalDict = {}
  for hist in histList:
    name = dataName(hist.GetName())
    if name.find('DATA') != -1:
      dataDict[name] = hist
      totalDict[name] = hist
    elif name.find('Signal') != -1:
      signalDict[name] = hist
      totalDict[name] = hist
    else:
      bgDict[name] = hist
      totalDict[name] = hist


  if len(dataDict) < 1 or len(bgDict) < 1: print 'data or bg DNE',plotName(histList[0].GetName());return
  canvas = TCanvas('canvas','canvas',1100,700)
  pad1 = TPad('pad1', '', 0.02, 0.34, 0.89, 0.98, 0)
  pad1.SetBottomMargin(0.07)
  pad1.SetGrid(1,1)
  if logToggle: pad1.SetLogy(kTRUE)
  pad2 = TPad('pad2', '', 0.02, 0.02, 0.89, 0.35, 0)
  pad2.SetTopMargin(0.)
  pad2.SetBottomMargin(0.2)
  pad2.SetGrid(1,1)
  pad1.Draw()
  pad2.Draw()
  stack = THStack('stack','')
  maxbinval = 0
  histmax = TH1F("temp","",histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1))
  errorHist = TH1F("errorHist","",histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1))
  for hist in bgDict:
    histmax.Add(bgDict[hist])

  maxbinval = histmax.GetBinContent(histmax.GetMaximumBin())
  #print maxbinval, histList[-1].GetBinContent(histmax.GetMaximumBin())
  if len(dataDict)>0:
    if dataDict.values()[0].GetBinContent(dataDict.values()[0].GetMaximumBin()) > maxbinval:
      maxbinval = dataDict.values()[0].GetBinContent(dataDict.values()[0].GetMaximumBin())

  if len(signalDict)>0:
    if signalDict.values()[0].GetBinContent(signalDict.values()[0].GetMaximumBin()) > maxbinval:
      maxbinval = signalDict.values()[0].GetBinContent(signalDict.values()[0].GetMaximumBin())

  tempList = sorted(bgDict.values(),key = lambda x: x.Integral())
  print tempList
  tempPairs = []
  for x in tempList:
    for y in bgDict:
      if bgDict[y] == x: tempPairs.append((y,x))
  print tempPairs
  bgDict = tempPairs

  for name, hist in bgDict:
    if len(bgDict) == 1:
      hist.SetOption('hist')
      hist.SetFillColor(ColorDict[name])
      hist.SetFillStyle(3001)
      hist.SetLineColor(ColorDict[name]+2)
      hist.SetLineStyle(2)
      hist.SetLineWidth(2)
      stack.Add(hist)
      stack.SetTitle(hist.GetTitle() + ';' + hist.GetXaxis().GetTitle() + ';' + hist.GetYaxis().GetTitle())
      for j in range(hist.GetNbinsX()+1):
        errorHist.SetBinContent(j,hist.GetBinContent(j))
        errorHist.SetBinError(j,hist.GetBinError(j))
        print hist.GetBinContent(j),hist.GetBinError(j),sqrt(hist.GetBinContent(j))
      errorHist.SetFillColor(ColorDict[dataName(hist.GetName())])
      errorHist.SetMarkerStyle(0)
    else:
      hist.SetOption('hist')
      hist.SetFillColor(ColorDict[name])
      print name
      hist.SetFillStyle(3001)
      stack.Add(hist)
  stack.SetTitle(bgDict[0][1].GetTitle() + ';' + bgDict[0][1].GetXaxis().GetTitle() + ';' + bgDict[0][1].GetYaxis().GetTitle())
  pname = plotName(histList[0].GetName())
  leg = makeLegend(histList, x1 = 0.89, x2 = 1.11, cname = cname)
  #else: leg = makeLegend([histList[0],histList[-1]])
  pad1.cd()
  h2_Temp = TH2F("h2_Temp",histList[0].GetTitle()+';'+histList[0].GetXaxis().GetTitle()+';'+histList[0].GetYaxis().GetTitle(),histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1),50,0,maxbinval+maxbinval*(0.1))
  h2_Temp.SetStats(0)
  if logToggle: h2_Temp.GetYaxis().SetLimits(0.01,stack.GetMaximum()*10) #print stack.GetMinimum(), stack.GetMaximum()
  h2_Temp.Draw()
  #print histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1)
  if len(bgDict) == 1: errorHist.Draw("same P E2")
  stack.Draw('same hist ][')
  dataDict.values()[0].SetMarkerStyle(20)
  dataDict.values()[0].SetOption('EP')
  dataDict.values()[0].Draw('same EP')
  if len(signalDict) > 0:
    signalDict.values()[0].SetLineColor(kRed+2)
    signalDict.values()[0].SetLineWidth(4)
    signalDict.values()[0].SetLineStyle(9)
    signalDict.values()[0].Draw('same hist')
  leg.Draw()
  pt = TPaveText(0.2,0.930,0.85,1,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(11)
  #pt.SetTextSize(2)
  text = pt.AddText("")
  pt.Draw()
  pad2.cd()
  nBins = histmax.GetNbinsX()
  xAxisName = histmax.GetXaxis().GetTitle()
  hRatio = TH1D('hRatio', ';'+histList[0].GetXaxis().GetTitle()+';Data/MC', nBins, histmax.GetBinLowEdge(1), histmax.GetBinLowEdge(nBins + 1))
  hRatio.Divide(histList[-1],histmax)

  hRatio.SetStats(0)
  hRatio.SetLineWidth(2);
  hRatio.SetLineColor(1);
  hRatio.SetMarkerStyle(20);
  hRatio.GetYaxis().SetLabelSize(0.08);
  hRatio.GetYaxis().SetTitleSize(0.09);
  hRatio.GetYaxis().SetTitleOffset(0.40);
  hRatio.GetYaxis().SetNdivisions(5);
  hRatio.GetYaxis().CenterTitle();
  hRatio.GetXaxis().SetLabelSize(0.08);
  hRatio.GetXaxis().SetTitleSize(0.09);
  hRatio.GetXaxis().SetTitleOffset(0.90);
  #hRatio.SetMaximum(min(hRatio.GetBinContent(hRatio.GetMaximumBin())*1.2,4));
  hRatio.SetMaximum(1.5);
  hRatio.SetMinimum(hRatio.GetBinContent(hRatio.GetMinimumBin())*0.8);
  fit = TF1("ratioFit","pol0")
  fit.SetLineWidth(3)
  fit.SetLineStyle(10)
  hRatio.Fit("ratioFit")
  hRatio.Draw()



  canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'.png')

def stacksAndSub(path,filename,histList,logToggle = False):
  today = datetime.date.today()
  year = today.strftime('%Y')
  month = today.strftime('%m').zfill(2)+today.strftime('-%B')
  day = today.strftime('%d').zfill(2)+today.strftime('-%A')
  path = year+'/'+month+'/'+day+'/'+path

  canvas = TCanvas('canvas','canvas',900,700)
  pad1 = TPad('pad1', '', 0.02, 0.34, 0.90, 0.98, 0)
  pad1.SetBottomMargin(0.07)
  pad1.SetGrid(1,1)
  if logToggle: pad1.SetLogy(kTRUE)
  pad2 = TPad('pad2', '', 0.02, 0.02, 0.90, 0.35, 0)
  pad2.SetTopMargin(0.)
  pad2.SetBottomMargin(0.2)
  pad2.SetGrid(1,1)
  pad1.Draw()
  pad2.Draw()
  stack = THStack('stack','')
  maxbinval = 0
  histmax = TH1F("temp","",histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1))
  for i,x in enumerate(histList[:-2]):
    if x!=None: histmax.Add(x)

  maxbinval = histmax.GetBinContent(histmax.GetMaximumBin())
  print maxbinval, histList[-1].GetBinContent(histmax.GetMaximumBin())
  if histList[-1].GetBinContent(histList[-1].GetMaximumBin()) > maxbinval:
    maxbinval = histList[-1].GetBinContent(histList[-1].GetMaximumBin())

  for i,hist in enumerate(histList[0:-2]):
    hist.SetOption('hist')
    hist.SetFillColor(ColorDict[dataName(hist.GetName())])
    #hist.SetFillStyle(FillList[i])
    stack.Add(hist)
    if i == 0:
      stack.SetTitle(hist.GetTitle() + ';' + hist.GetXaxis().GetTitle() + ';' + hist.GetYaxis().GetTitle())
  leg = makeLegend(histList)
  pad1.cd()
  h2_Temp = TH2F("h2_Temp","",histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1),50,0,maxbinval+maxbinval*(0.1))
  h2_Temp.SetStats(0)
  h2_Temp.Draw()
  #print histList[0].GetNbinsX(),histList[0].GetBinLowEdge(1),histList[0].GetBinLowEdge(histList[0].GetNbinsX()+1)
  stack.Draw('same hist')
  histList[-1].SetMarkerStyle(20)
  histList[-1].SetOption('EP')
  histList[-1].Draw('same EP')
  histList[-2].SetLineColor(kRed+2)
  histList[-2].SetLineWidth(4)
  histList[-2].SetLineStyle(9)
  histList[-2].Draw('same hist')
  leg.Draw()
  pt = TPaveText(0.2,0.930,0.85,1,"blNDC")
  pt.SetName("title")
  pt.SetBorderSize(0)
  pt.SetFillColor(0)
  pt.SetFillStyle(0)
  pt.SetTextAlign(11)
  #pt.SetTextSize(2)
  text = pt.AddText("")
  pt.Draw()
  pad2.cd()
  histMC = histList[0]
  for hist in histList[1:-2]:
    histMC = histMC + hist
  nBins = histMC.GetNbinsX()
  xAxisName = histMC.GetXaxis().GetTitle()
  hRatio = TH1D('hRatio', ';'+xAxisName+';Data/BG', nBins, histMC.GetBinLowEdge(1), histMC.GetBinLowEdge(nBins + 1))
  hRatio = histList[-1] - histMC

  hRatio.SetTitle('')
  hRatio.SetStats(0)
  hRatio.SetLineWidth(2);
  hRatio.SetLineColor(1);
  hRatio.SetMarkerStyle(20);
  hRatio.GetYaxis().SetLabelSize(0.08);
  hRatio.GetYaxis().SetTitleSize(0.09);
  hRatio.GetYaxis().SetTitle('Data - BG')
  hRatio.GetYaxis().SetTitleOffset(0.40);
  hRatio.GetYaxis().SetNdivisions(5);
  hRatio.GetYaxis().CenterTitle();
  hRatio.GetXaxis().SetLabelSize(0.08);
  hRatio.GetXaxis().SetTitleSize(0.09);
  hRatio.GetXaxis().SetTitleOffset(0.90);
  hRatio.SetMaximum(10);
  hRatio.SetMinimum(-10);
  hRatio.Draw()

  canvas.SaveAs('~/nobackup/TempWeb/'+path+'/'+filename+'sub.png')


ColorDict = {'ZJets':kGreen+2, 'ZJetsBG':kGreen+2, 'ZJetsBG2012':kGreen+2, 'ZJetsBG2011':kGreen+2,
    'ZZJets':kMagenta+2,'WZJets':kBlue,'WWJets':kRed-6,'ZGMuMu':kCyan-2,'ZGToLLGBG2012':kCyan-2,'ZGToLLGBG2011':kCyan-2,'ZGMuMu2j':kCyan-1,
    'GluGluWW':kOrange-8,'ZGFrac':kViolet+1,'DYMuMu':kTeal+1,'tW':kYellow+2,'ttbar':kTeal-5,'ZGEE':kCyan-2,
    'TTJets':kViolet,'WJets':kGray+2,'HZG125Signal':kRed,'ggHZG125Signal':kRed, 'vbfHZG125Signal':kRed,
    'OFFggHZG125Signal':kRed, 'OFFvbfHZG125Signal':kRed,
    'HZG135Signal':kRed,'DATA':kBlack,
    'GammaJet15-30':kBlue-8,
    'GammaJet30-50':kBlue-8,
    'GammaJet50-80':kBlue-8,
    'GammaJet80-120':kBlue-8,
    'GammaJet120-170':kBlue-8,
    'GammaJet170-300':kBlue-8,
    'GammaJet300-470':kBlue-8,
    'GammaJet470-800':kBlue-8,
    }
MarkerList = [20,21,22,23,33,34,29,3,32,27]
FillList = [3004,3005,3006,3007,3020,3244,3023]

if __name__ == "__main__":

  import sys, string

  funcname, args = sys.argv[1], sys.argv[2:]

  arglist = string.join(args, ", ")
  cmd = "%s('%s')" % (funcname, arglist)
  print "Executing:", cmd

  result = eval(cmd)

  print "Result:", result

