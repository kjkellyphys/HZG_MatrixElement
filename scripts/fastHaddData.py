#!/usr/bin/env python
import sys, os

#useage:
# ./fastHadd.py year lepton LowR9(HighR9)
if len(sys.argv) !=4:
  print 'you are doing this wrong'
else:

  year = sys.argv[1]
  if sys.argv[2].lower() == 'mumu':
    leptonA = 'MuMu'
    leptonB = 'mumuGamma'
  elif sys.argv[2].lower() == 'ee':
    leptonA = 'EE'
    leptonB = 'eeGamma'
  R9 = sys.argv[3]

  selectionString = ''
  if 'A' in year:
    selectionString += ' Run2012A_Jul13 Run2012A_Aug06rec'
  if 'B' in year:
    selectionString += ' Run2012B_Jul13'
  if 'C' in year:
    selectionString += ' Run2012C_Aug24 Run2012C_Dec11 Run2012C_PromptV2'
  if 'D' in year:
    selectionString += ' Run2012D_PromptV1'

  if '2012' in year:
    os.system('./hadd.py ForStoyan/higgsHistograms_'+leptonA+year+'_Lite_'+R9+'.root '+leptonB+' File'+selectionString)
    os.system('./hadd.py ~/afsHome/public/m_llgFile_'+leptonA+year+'_'+R9+'.root '+leptonB+' m_llgFile'+selectionString)
    os.system('./hadd.py batchHistos/higgsHistograms_'+leptonA+year+'_'+R9+'.root '+leptonB+' Histograms'+selectionString)
  else:
    os.system('./hadd.py ForStoyan/higgsHistograms_'+leptonA+year+'_Limits_'+R9+'.root '+leptonB+' File Run'+year+'A Run'+year+'B ggHZG_M120 ggHZG_M125 ggHZG_M130 ggHZG_M135 ggHZG_M140 ggHZG_M145 ggHZG_M150 vbfHZG_M120 vbfHZG_M125 vbfHZG_M130 vbfHZG_M135 vbfHZG_M140 vbfHZG_M145 vbfHZG_M150 tthHZG_M120 tthHZG_M125 tthHZG_M130 tthHZG_M135 tthHZG_M140 tthHZG_M145 tthHZG_M150 whHZG_M120 whHZG_M125 whHZG_M130 whHZG_M135 whHZG_M140 whHZG_M145 whHZG_M150 zhHZG_M120 zhHZG_M125 zhHZG_M130 zhHZG_M135 zhHZG_M140 zhHZG_M145 zhHZG_M150')
    os.system('./hadd.py ~/afsHome/public/m_llgFile_'+leptonA+year+'_'+R9+'.root '+leptonB+' m_llgFile Run'+year+'A Run'+year+'B ggHZG_M120 ggHZG_M125 ggHZG_M130 ggHZG_M135 ggHZG_M140 ggHZG_M145 ggHZG_M150 vbfHZG_M120 vbfHZG_M125 vbfHZG_M130 vbfHZG_M135 vbfHZG_M140 vbfHZG_M145 vbfHZG_M150 tthHZG_M120 tthHZG_M125 tthHZG_M130 tthHZG_M135 tthHZG_M140 tthHZG_M145 tthHZG_M150 whHZG_M120 whHZG_M125 whHZG_M130 whHZG_M135 whHZG_M140 whHZG_M145 whHZG_M150 zhHZG_M120 zhHZG_M125 zhHZG_M130 zhHZG_M135 zhHZG_M140 zhHZG_M145 zhHZG_M150')
    os.system('./hadd.py batchHistos/higgsHistograms_'+leptonA+year+'_'+R9+'.root '+leptonB+' Histograms Run'+year+'A Run'+year+'B ggHZG_M120 ggHZG_M125 ggHZG_M130 ggHZG_M135 ggHZG_M140 ggHZG_M145 ggHZG_M150 vbfHZG_M120 vbfHZG_M125 vbfHZG_M130 vbfHZG_M135 vbfHZG_M140 vbfHZG_M145 vbfHZG_M150 tthHZG_M120 tthHZG_M125 tthHZG_M130 tthHZG_M135 tthHZG_M140 tthHZG_M145 tthHZG_M150 whHZG_M120 whHZG_M125 whHZG_M130 whHZG_M135 whHZG_M140 whHZG_M145 whHZG_M150 zhHZG_M120 zhHZG_M125 zhHZG_M130 zhHZG_M135 zhHZG_M140 zhHZG_M145 zhHZG_M150')


#hadd -f otherHistos/eleSmear2011.root ~/nobackup/BatchOutput/eeGamma_Combined/eleSmearFile_*
#./hadd.py ForStoyan/higgsHistograms_MuMu2011_Limits_Nominal_LowR9.root mumuGamma File Run2011A Run2011B ggHZG_M120 ggHZG_M125 ggHZG_M130 ggHZG_M135 ggHZG_M140 vbfHZG_M120 vbfHZG_M125 vbfHZG_M130 vbfHZG_M135 vbfHZG_M140

