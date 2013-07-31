#! /usr/bin/env python
import BatchMaster as b
import os

dCache      = '/pnfs/cms/WAX/11/store/user'
EOS         = '/eos/uscms/store/user'
outputPath  = '/uscms/home/bpollack/nobackup/BatchOutput'

selection   = 'mumuGamma'
period      = '2012'

bg = []
bg.extend([
  b.JobConfig('ZJets', EOS+'/devildog/V01_10_8TeV/DYJetsToLLV2', 60, 'ZJetsBG'+period+' 3,4 '+selection+' '+period, 'test'),
])

os.system('rm outputBatch/*')
batcher = b.BatchMaster(bg, outputPath,'execBatch.csh')
batcher.SubmitToLPC()
