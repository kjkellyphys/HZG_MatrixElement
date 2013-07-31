#!/usr/bin/env python
import sys, os
from math import *
from operator import mul
import string

def stringFinder(txtFile):
  f = open(txtFile, 'r')
  allNum = set()
  dupNum = []
  for line in f:
    num=line.strip().lstrip('/uscms/home/bpollack/dCache/Apr17/MC/ZJets/nuTuple').lstrip('_').rstrip('.root').rstrip(string.ascii_lowercase+string.ascii_uppercase+'0123456789').rstrip('_').rstrip('0123').rstrip('_')
    if num != '':
      num = int(num)
      print num
    else: continue
    if num in allNum:
      dupNum.append(num)

    allNum.add(num)
  return dupNum



if __name__=="__main__":
  if len(sys.argv)==2:
    answer = stringFinder(sys.argv[1])
  else:
    answer = stringFinder()
  print answer
