#!/usr/bin/env python
from ROOT import *
import sys

f = open(sys.argv[1], 'r')
numTotal = 0
fileNum = 0
for line in f:
  fileNum += 1
  tmp =  line.lstrip('file:')
  tmp = tmp.rstrip('\n')
  x = TFile(tmp,'OPEN')
  h1_nEvents = x.Get("ntupleProducer/numOfEvents")
  num = h1_nEvents.GetBinContent(1)
  numTotal +=num
  print 'file number:',fileNum,num
  x.Close()
print numTotal
