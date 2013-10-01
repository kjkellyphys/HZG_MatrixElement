#! /usr/bin/env python

from ROOT import *

myFile = TFile('pp_zg_background_AnglePlots.root','r')

histoList = ['WD', 'WD_A', 'WD_B', 'WD_C', 'WD_D']

canvas = TCanvas()

for i, hist in enumerate(histoList):
    myHist = myFile.Get(hist)
    myHist.SetLineColor(i+2)
    if i is 0:
        myHist.Draw()
    else:
        myHist.Draw('same')

canvas.SaveAs('kevinsFirstPYProg.pdf')
