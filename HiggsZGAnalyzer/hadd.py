#!/usr/bin/env python
import sys, os

# use it like this: ./hadd.py data/higgsHistograms_EE_allBG-05-23-12.root eeGamma Histograms ZJets ZGEE WWJets ZZJets WZJets HZG125 Run2011A Run2011B


def hadd(output = 'data/temp.root', inSelection = 'mumuGamma', inFileType = 'Histogram', *inputs):
  batchPath = '~/nobackup/BatchOutput/'
  inpaths = ''
  for i in inputs:
    if inFileType == 'm_llgFile':
      inpaths += batchPath+inSelection+'_Combined/'+inFileType+'_'+i+'_* '
    elif 'ele' in inFileType:
      inpaths += batchPath+inSelection+'_Combined/'+inFileType+'File_'+i+'_* '
    else:
      inpaths += batchPath+inSelection+'_Combined/'+'higgs'+inFileType+'_'+i+'_* '
  os.system('hadd -f '+output+' '+inpaths)
  #print output, inputs

#hadd -f otherHistos/eleSmear2011.root ~/nobackup/BatchOutput/eeGamma_Combined/eleSmearFile_*




if __name__=="__main__":
  if len(sys.argv)>=2:
    answer = hadd(*sys.argv[1:])
  else:
    #answer = hadd()
    print "You need more arguments, asshole"
