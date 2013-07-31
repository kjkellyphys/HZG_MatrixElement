import sys, os, subprocess, fileinput, pickle, math, tempfile, shutil, glob

class JobConfig():
    '''Class for storing configuration for each dataset'''
    def __init__(self, dataName = 'TEST', inDir = '/uscms/home/naodell/nobackup/TEST', nJobs = 1, arguments = 'TEST 0 muon', selection = 'muon'):
        self._dataName  = dataName
        self._inDir     = inDir
        self._nJobs     = nJobs
        self._args      = arguments
        self._selection = selection

    def GetMembers(self, verbose = False):
        '''Returns data members of class instance'''
        if verbose:
            print 'dataset: %s \t number of jobs: %i \t triggers: %s' % (self._dataName, self._nJobs, self._triggers)
        return (self._dataName, self._inDir, self._nJobs, self._triggers)

class BatchMaster():
    '''A tool for submitting batch jobs'''
    def __init__(self, configList, outDir = '.',execution = 'execBatch.csh'):
        self._current    = os.path.abspath('.')
        self._outDir     = outDir
        self._configList = configList
        self._executable = execution

    def MakeDirectory(self, filePath, clear = True):
        '''Create save path in case it doesn't already exist'''
        if not os.path.exists(filePath):
            os.system('mkdir -p '+filePath)
        if not os.path.exists(filePath+'/res'):
            os.system('mkdir -p '+filePath+'/res')
        elif clear:
            os.system('rm '+filePath+'/*')

    def SplitJobs(self, directory, nJobs):
        '''Split jobs by dataset'''
        fileList = os.listdir(directory)
        nFiles = len(fileList)

        if nJobs > nFiles:
            nJobs = nFiles
        nFilesPerJob = int(math.ceil(nFiles/nJobs))
        fileSplit = [fileList[i:i+nFilesPerJob] for i in range(0, len(fileList), nFilesPerJob)]

        return fileSplit

    def MakeExecutable(self, config, sourceFiles, count):
        '''Writes config into executable'''
        infile   = open(self._executable, 'r')
        exec_tmp = tempfile.NamedTemporaryFile(prefix = config._dataName+'_'+config._selection+'_', delete=False)

        for i, line in enumerate(infile.readlines()):
            if i != 6:
                exec_tmp.write(line)
            else:
                path = ''
                if config._inDir[:5] == '/pnfs':
                  path = 'dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr'+config._inDir[5:]
                else:
                  path = config._inDir
                for file in sourceFiles:
                    exec_tmp.write('echo '+repr(path+'/'+file)+' >> input.txt\n')

        exec_tmp.seek(0)
        infile.close()
        return exec_tmp

    def MakeBatchConfig(self, config, count, exec_tmp, sourceFiles):
        '''Write batch configuration file'''
        batch_tmp = 'Arguments  = %s %s %s' % (str(count+1), config._dataName, config._args)
        batch_tmp = batch_tmp+'\nExecutable = %s' % exec_tmp.name
        batch_tmp = batch_tmp+'\nShould_Transfer_Files = YES'
        batch_tmp = batch_tmp+'\nWhenToTransferOutput = ON_EXIT'
        batch_tmp = batch_tmp+'\nTransfer_Input_Files = stageball.tar.gz'

        batch_tmp = batch_tmp+'''
universe              = vanilla
Requirements          = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
+LENGTH               = "LONG"
Output                = res/report_$(Cluster)_$(Process).stdout
Error                 = res/report_$(Cluster)_$(Process).stderr
Log                   = res/report_$(Cluster)_$(Process).log
notify_user           = Nathaniel.Jay.Odell@cern.ch
Queue
'''
        batch_submit_file = tempfile.NamedTemporaryFile(prefix = 'batch_submit_file', delete=False)
        batch_submit_file.write(batch_tmp)
        batch_submit_file.seek(0)
        return batch_submit_file

    def SubmitToLPC(self):
        '''Submits batch jobs to lpc batch'''
        self.MakeDirectory(self._outDir, clear = False)
        for j,cfg in enumerate(self._configList):
            self.MakeDirectory(self._outDir+'/'+cfg._selection, clear = False)

            # copy files to staging so they cant be hurt during submission
            sourceFiles = self.SplitJobs(cfg._inDir, cfg._nJobs)

            if j ==0:

              if not os.path.exists(self._outDir+'/'+cfg._selection+'/stageball.tar.gz'):
                os.system('cp stageball.tar.gz '+self._outDir+'/'+cfg._selection+'/.' )
              if not os.path.exists(self._outDir+'/'+cfg._selection+'/execBatch.csh'):
                os.system('cp execBatch.csh '+self._outDir+'/'+cfg._selection+'/.' )
              os.chdir(self._outDir+'/'+cfg._selection)
              #os.system('tar -zxf stageball.tar.gz')
            #os.system('cd '+self._outDir+'/'+cfg._selection)

            for i, source in enumerate(sourceFiles):
                exec_tmp  = self.MakeExecutable(cfg, source, i)
                batch_tmp = self.MakeBatchConfig(cfg, i, exec_tmp, source)
                subprocess.call('condor_submit ' + batch_tmp.name , shell=True)
                exec_tmp.close()
                batch_tmp.close()
