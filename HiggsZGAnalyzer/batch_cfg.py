#! /usr/bin/env python
import BatchMaster as b
import os

dCache      = '/pnfs/cms/WAX/11/store/user'
EOS         = '/eos/uscms/store/user'
FuckEOS     = '/uscms_data/d2/bpollack/FuckEOS'
outputPath  = '/uscms/home/bpollack/nobackup/BatchOutput'
#outputPath  = '/eos/uscms/store/user/bpollack/BatchOutput'

#selection   = 'mumuGamma'
selection   = 'eeGamma'

period      = '2012'
#period      = '2011'

#period2     = 'AB'
#period2     = 'CD'
period2     = 'ABCD'

doBG     = False
doSignal = True
doData   = True



if selection == 'mumuGamma':
    data = []
    if period in ['2011', 'Combined']:
      data.append(b.JobConfig('Run2011A', EOS+'/lusito/V01_02_7TeV/DoubleMu_Run2011A', 40, 'DATA A  mumuGamma 2011A', selection+'_'+'Combined'))
      data.append(b.JobConfig('Run2011B', EOS+'/lusito/V01_02_7TeV/DoubleMu_Run2011B', 40, 'DATA B mumuGamma 2011B', selection+'_'+'Combined'))
    if period in ['2012', 'Combined']:
      if period2 is 'AB':
        data.append(b.JobConfig('Run2012A_Jul13',          EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012A_Jul13', 20, 'DATA A  mumuGamma 2012A_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012A_Aug06rec',       EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012A_Aug06rec', 10, 'DATA A  mumuGamma 2012A_Aug06rec', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012B_Jul13',          EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012B_Jul13', 50, 'DATA B mumuGamma 2012B_Jul13', selection+'_'+'Combined'))
      elif period2 is 'CD':
        data.append(b.JobConfig('Run2012C_Aug24',    EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012C_Aug24', 20,    'DATA C  mumuGamma 2012C_Aug24', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Dec11',    EOS +'/bpollack/V02_01_8TeV/DoubleMu_Run2012C_Dec11', 10,    'DATA C  mumuGamma 2012C_Dec11', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_PromptV2', EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012C_PromptV2', 50, 'DATA C mumuGamma 2012C_PromptV2', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012D_PromptV1', EOS +'/lusito/V02_01_8TeV/DoubleMu_Run2012D_PromptV1', 50,   'DATA D mumuGamma 2012D_PromptV1', selection+'_'+'Combined'))
      elif period2 is 'ABCD':
        data.append(b.JobConfig('Run2012A_Jul13',          EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012A_Jul13', 20, 'DATA A  mumuGamma 2012A_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012A_Aug06rec',       EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012A_Aug06rec', 10, 'DATA A  mumuGamma 2012A_Aug06rec', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012B_Jul13',          EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012B_Jul13', 50, 'DATA B mumuGamma 2012B_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Aug24',    EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012C_Aug24', 20,    'DATA C  mumuGamma 2012C_Aug24', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Dec11',    EOS +'/bpollack/V02_01_8TeV/DoubleMu_Run2012C_Dec11', 10,    'DATA C  mumuGamma 2012C_Dec11', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_PromptV2', EOS +'/devildog/V02_01_8TeV/DoubleMu_Run2012C_PromptV2', 50, 'DATA C mumuGamma 2012C_PromptV2', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012D_PromptV1', EOS +'/lusito/V02_01_8TeV/DoubleMu_Run2012D_PromptV1', 50,   'DATA D mumuGamma 2012D_PromptV1', selection+'_'+'Combined'))


if selection == 'eeGamma':
    data = []
    if period in ['2011', 'Combined']:
      data.append(b.JobConfig('Run2011A', EOS+'/devildog/V01_02_7TeV/DoubleEl_Run2011A', 60, 'DATA 17,18 eeGamma 2011A', selection+'_'+'Combined'))
      data.append(b.JobConfig('Run2011B', EOS+'/devildog/V01_02_7TeV/DoubleEl_Run2011B', 40, 'DATA 17,18 eeGamma 2011B', selection+'_'+'Combined'))
    if period in ['2012', 'Combined']:
      if period2 is 'AB':
        data.append(b.JobConfig('Run2012A_Jul13',          EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012A_Jul13', 20, 'DATA A  eeGamma 2012A_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012A_Aug06rec',       EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012A_Aug06rec', 10,     'DATA A  eeGamma 2012A_Aug06rec', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012B_Jul13',          EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012B_Jul13', 50,        'DATA B eeGamma 2012B_Jul13', selection+'_'+'Combined'))
      if period2 is 'CD':
        data.append(b.JobConfig('Run2012C_Aug24',    EOS +'/devildog/V02_01_8TeV/DoubleEl_Run2012C_Aug24', 20,    'DATA C eeGamma 2012C_Aug24', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Dec11',    EOS +'/lusito/V02_01_8TeV/DoubleEl_Run2012C_Dec11', 10,    'DATA C eeGamma 2012C_Dec11', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_PromptV2', EOS +'/lusito/V02_01_8TeV/DoubleEl_Run2012C_PromptV2', 50,   'DATA C eeGamma 2012C_PromptV2', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012D_PromptV1', EOS +'/bpollack/V02_01_8TeV/DoubleEl_Run2012D_PromptV1', 50, 'DATA D eeGamma 2012D_PromptV1', selection+'_'+'Combined'))
      if period2 is 'ABCD':
        data.append(b.JobConfig('Run2012A_Jul13',          EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012A_Jul13', 20, 'DATA A  eeGamma 2012A_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012A_Aug06rec',       EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012A_Aug06rec', 10,     'DATA A  eeGamma 2012A_Aug06rec', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012B_Jul13',          EOS +'/lusito/V02_01_8TeV/DoublEl_Run2012B_Jul13', 50,        'DATA B eeGamma 2012B_Jul13', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Aug24',    EOS +'/devildog/V02_01_8TeV/DoubleEl_Run2012C_Aug24', 20,    'DATA C eeGamma 2012C_Aug24', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_Dec11',    EOS +'/lusito/V02_01_8TeV/DoubleEl_Run2012C_Dec11', 10,    'DATA C eeGamma 2012C_Dec11', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012C_PromptV2', EOS +'/lusito/V02_01_8TeV/DoubleEl_Run2012C_PromptV2', 50,   'DATA C eeGamma 2012C_PromptV2', selection+'_'+'Combined'))
        data.append(b.JobConfig('Run2012D_PromptV1', EOS +'/bpollack/V02_01_8TeV/DoubleEl_Run2012D_PromptV1', 50, 'DATA D eeGamma 2012D_PromptV1', selection+'_'+'Combined'))

bg = []
if period == '2011':
  bg.extend([
    b.JobConfig('ZJets', dCache+'/bpollack/V01_01_7TeV/DYJetsToLL', 60, 'ZJetsBG'+period+' '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
  ])
if period == '2012':
  bg.extend([
    b.JobConfig('ZJets', dCache+'/bpollack/V01_10_8TeV/DYJetsToLL', 60, 'ZJetsBG'+period+' '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    #b.JobConfig('ZGToLLG', EOS+'/devildog/V01_10_8TeV/ZGToLLG', 20, 'ZGToLLGBG'+period+triggers+selection+' '+period, selection+'_'+'Combined'),
  ])


signalZG = []
if period == '2011':
  signalZG.extend([
    b.JobConfig('ggHZG_M120', EOS+'/devildog/V01_02_7TeV/ggHZG_M120_7TeV', 2, 'Signal'+period+'ggM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M125', EOS+'/devildog/V01_02_7TeV/ggHZG_M125_7TeV', 2, 'Signal'+period+'ggM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M130', EOS+'/devildog/V01_02_7TeV/ggHZG_M130_7TeV', 2, 'Signal'+period+'ggM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M135', EOS+'/devildog/V01_02_7TeV/ggHZG_M135_7TeV', 2, 'Signal'+period+'ggM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M140', EOS+'/devildog/V01_02_7TeV/ggHZG_M140_7TeV', 2, 'Signal'+period+'ggM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M145', EOS+'/devildog/V01_02_7TeV/ggHZG_M145_7TeV', 2, 'Signal'+period+'ggM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M150', EOS+'/devildog/V01_02_7TeV/ggHZG_M150_7TeV', 2, 'Signal'+period+'ggM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('vbfHZG_M120', EOS+'/lusito/V01_02_7TeV/vbfHZG_M120_7TeV', 2, 'Signal'+period+'vbfM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M125', EOS+'/lusito/V01_02_7TeV/vbfHZG_M125_7TeV', 2, 'Signal'+period+'vbfM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M130', EOS+'/lusito/V01_02_7TeV/vbfHZG_M130_7TeV', 2, 'Signal'+period+'vbfM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M135', EOS+'/lusito/V01_02_7TeV/vbfHZG_M135_7TeV', 2, 'Signal'+period+'vbfM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M140', EOS+'/lusito/V01_02_7TeV/vbfHZG_M140_7TeV', 2, 'Signal'+period+'vbfM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M145', EOS+'/lusito/V01_02_7TeV/vbfHZG_M145_7TeV', 2, 'Signal'+period+'vbfM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M150', EOS+'/lusito/V01_02_7TeV/vbfHZG_M150_7TeV', 2, 'Signal'+period+'vbfM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('tthHZG_M120', EOS+'/devildog/V01_02_7TeV/tthHZG_M120_7TeV', 2, 'Signal'+period+'tthM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M125', EOS+'/devildog/V01_02_7TeV/tthHZG_M125_7TeV', 2, 'Signal'+period+'tthM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M130', EOS+'/devildog/V01_02_7TeV/tthHZG_M130_7TeV', 2, 'Signal'+period+'tthM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M135', EOS+'/devildog/V01_02_7TeV/tthHZG_M135_7TeV', 2, 'Signal'+period+'tthM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M140', EOS+'/devildog/V01_02_7TeV/tthHZG_M140_7TeV', 2, 'Signal'+period+'tthM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M145', EOS+'/devildog/V01_02_7TeV/tthHZG_M145_7TeV', 2, 'Signal'+period+'tthM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M150', EOS+'/devildog/V01_02_7TeV/tthHZG_M150_7TeV', 2, 'Signal'+period+'tthM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('whHZG_M120', EOS+'/devildog/V01_02_7TeV/whzhHZG_M120_7TeV', 6, 'Signal'+period+'whM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M125', EOS+'/devildog/V01_02_7TeV/whzhHZG_M125_7TeV', 6, 'Signal'+period+'whM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M130', EOS+'/devildog/V01_02_7TeV/whzhHZG_M130_7TeV', 6, 'Signal'+period+'whM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M135', EOS+'/devildog/V01_02_7TeV/whzhHZG_M135_7TeV', 6, 'Signal'+period+'whM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M140', EOS+'/devildog/V01_02_7TeV/whzhHZG_M140_7TeV', 6, 'Signal'+period+'whM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M145', EOS+'/devildog/V01_02_7TeV/whzhHZG_M145_7TeV', 6, 'Signal'+period+'whM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M150', EOS+'/devildog/V01_02_7TeV/whzhHZG_M150_7TeV', 6, 'Signal'+period+'whM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('zhHZG_M120', EOS+'/devildog/V01_02_7TeV/whzhHZG_M120_7TeV', 6, 'Signal'+period+'zhM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M125', EOS+'/devildog/V01_02_7TeV/whzhHZG_M125_7TeV', 6, 'Signal'+period+'zhM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M130', EOS+'/devildog/V01_02_7TeV/whzhHZG_M130_7TeV', 6, 'Signal'+period+'zhM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M135', EOS+'/devildog/V01_02_7TeV/whzhHZG_M135_7TeV', 6, 'Signal'+period+'zhM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M140', EOS+'/devildog/V01_02_7TeV/whzhHZG_M140_7TeV', 6, 'Signal'+period+'zhM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M145', EOS+'/devildog/V01_02_7TeV/whzhHZG_M145_7TeV', 6, 'Signal'+period+'zhM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M150', EOS+'/devildog/V01_02_7TeV/whzhHZG_M150_7TeV', 6, 'Signal'+period+'zhM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

  ])
if period == '2012':
  signalZG.extend([
    b.JobConfig('ggHZG_M120', EOS+'/bpollack/V02_01_8TeV/ggHZG_M120', 2, 'Signal'+period+'ggM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M125', EOS+'/bpollack/V02_01_8TeV/ggHZG_M125', 2, 'Signal'+period+'ggM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M130', EOS+'/bpollack/V02_01_8TeV/ggHZG_M130', 2, 'Signal'+period+'ggM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M135', EOS+'/bpollack/V02_01_8TeV/ggHZG_M135', 2, 'Signal'+period+'ggM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M140', EOS+'/bpollack/V02_01_8TeV/ggHZG_M140', 2, 'Signal'+period+'ggM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M145', EOS+'/bpollack/V02_01_8TeV/ggHZG_M145', 2, 'Signal'+period+'ggM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('ggHZG_M150', EOS+'/bpollack/V02_01_8TeV/ggHZG_M150', 2, 'Signal'+period+'ggM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('vbfHZG_M120', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M120', 2, 'Signal'+period+'vbfM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M125', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M125', 2, 'Signal'+period+'vbfM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M130', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M130', 2, 'Signal'+period+'vbfM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M135', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M135', 2, 'Signal'+period+'vbfM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M140', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M140', 2, 'Signal'+period+'vbfM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M145', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M145', 2, 'Signal'+period+'vbfM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('vbfHZG_M150', EOS+'/bpollack/V02_01_8TeV/vbfHZG_M150', 2, 'Signal'+period+'vbfM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('tthHZG_M120', EOS+'/bpollack/V02_01_8TeV/tthHZG_M120', 2, 'Signal'+period+'tthM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M125', EOS+'/bpollack/V02_01_8TeV/tthHZG_M125', 2, 'Signal'+period+'tthM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M130', EOS+'/bpollack/V02_01_8TeV/tthHZG_M130', 2, 'Signal'+period+'tthM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M135', EOS+'/bpollack/V02_01_8TeV/tthHZG_M135', 2, 'Signal'+period+'tthM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M140', EOS+'/bpollack/V02_01_8TeV/tthHZG_M140', 2, 'Signal'+period+'tthM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M145', EOS+'/bpollack/V02_01_8TeV/tthHZG_M145', 2, 'Signal'+period+'tthM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('tthHZG_M150', EOS+'/bpollack/V02_01_8TeV/tthHZG_M150', 2, 'Signal'+period+'tthM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('whHZG_M120', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M120', 2, 'Signal'+period+'whM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M125', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M125', 2, 'Signal'+period+'whM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M130', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M130', 2, 'Signal'+period+'whM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M135', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M135', 2, 'Signal'+period+'whM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M140', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M140', 2, 'Signal'+period+'whM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M145', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M145', 2, 'Signal'+period+'whM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('whHZG_M150', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M150', 2, 'Signal'+period+'whM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

    b.JobConfig('zhHZG_M120', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M120', 2, 'Signal'+period+'zhM120 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M125', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M125', 2, 'Signal'+period+'zhM125 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M130', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M130', 2, 'Signal'+period+'zhM130 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M135', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M135', 2, 'Signal'+period+'zhM135 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M140', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M140', 2, 'Signal'+period+'zhM140 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M145', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M145', 2, 'Signal'+period+'zhM145 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),
    b.JobConfig('zhHZG_M150', EOS+'/bpollack/V02_01_8TeV/whzhHZG_M150', 2, 'Signal'+period+'zhM150 '+period2+' '+selection+' '+period, selection+'_'+'Combined'),

  ])

os.system('rm -r '+outputPath+'/'+selection+'_Combined/*')

#os.system('rm outputBatch/*')
os.system('rm dumps/dataDump_'+selection+'.txt')

os.system('tar -zcvf stageball.tar.gz higgsAnalyzer_Template.C higgsAnalyzer.h ../src/*.cc ../src/*.h otherHistos/*.root plugins/*.C plugins/*.cc plugins/*.h* plugins/PHOSPHOR_NUMBERS_EXPFIT.txt plugins/PHOSPHOR_NUMBERS_EXPFIT_ERRORS.txt plugins/PHOSPHOR_Corr plugins/MVAWeights/*')
if doData:
  batcher = b.BatchMaster(data, outputPath,'execBatch.csh')
  batcher.SubmitToLPC()
if doBG:
  batcher = b.BatchMaster(bg, outputPath,'execBatch.csh')
  batcher.SubmitToLPC()
if doSignal:
  batcher = b.BatchMaster(signalZG, outputPath,'execBatch.csh')
  batcher.SubmitToLPC()








