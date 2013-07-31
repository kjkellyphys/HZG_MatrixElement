#!/usr/bin/env python
import sys, os
import subprocess

def checkFiles(selection):
  if selection == 'mumu':
    proc = subprocess.Popen(["bash -c 'wc -l <(ls -1 ~/nobackup/BatchOutput/mumuGamma_Combined)'"],stdout=subprocess.PIPE, shell=True)
    (out,err) = proc.communicate()
    proc2 = subprocess.Popen(["bash -c 'wc -l <(ls -1 ~/nobackup/BatchOutput/mumuGamma_Combined/res)'"],stdout=subprocess.PIPE, shell=True)
    (out2,err2) = proc2.communicate()
  elif selection == 'ee':
    proc = subprocess.Popen(["bash -c 'wc -l <(ls -1 ~/nobackup/BatchOutput/eeGamma_Combined)'"],stdout=subprocess.PIPE, shell=True)
    (out,err) = proc.communicate()
    proc2 = subprocess.Popen(["bash -c 'wc -l <(ls -1 ~/nobackup/BatchOutput/eeGamma_Combined/res)'"],stdout=subprocess.PIPE, shell=True)
    (out2,err2) = proc2.communicate()
  else: print 'choose ee or mumu'; return
  num1 = int(out.rsplit(" ")[0])-1
  num2 = int(out2.rsplit(" ")[0])
  print num1/7, num2/3
  if num1/7 == num2/3: print 'ALL JOBS FINISHED'
  else: print "SOME JOBS HAVEN'T FUCKING FINISHED GODDAMNIT!"






if __name__=="__main__":
  if len(sys.argv)==2:
    checkFiles(sys.argv[1])
  else:
    #answer = hadd()
    print "You need 1 argument, asshole"
