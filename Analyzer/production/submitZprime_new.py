#!/usr/bin/env python

import sys, commands, os, fnmatch
from optparse import OptionParser


###### need to choose the json
#    jsonPrompt = "$PWD/../data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
#    jsonRereco = "$PWD/../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
#    jsonPrompt17 = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
#    jsonRereco17 = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
json = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
    
##### normalize MC to 1pb ?
xsec = 1

        
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--test',dest="test",default=False,action='store_true',help="Just print out commands to run")    
    parser.add_option('-s','--sample',dest="sample", default="GluGluHToBB_M125_13TeV_powheg_pythia8",help="sample to produce")
    parser.add_option('-e','--executable',dest="executable", default="runZprime", help = "executable name")
    parser.add_option('-t','--tag',dest="tag", default = "zprimebits-v12.08", help = "tag, which is the same as folder")
    parser.add_option('--sub',dest="sub",default=False,action='store_true',help="submit condor jobs")    
    (options,args) = parser.parse_args()


    isMC = True
    if "Run20" in options.sample:
        isMC = False

    ###### intput root file list
    inputFiles = os.environ['CMSSW_BASE']+"/src/BaconAnalyzer/Analyzer/lists/production14/"+ options.sample +".txt"
    print inputFiles

    ##### EOS  Output directories 
    # eosOutDir = '/store/user/lpcbacon/dazsle'
    eosOutDir = '/store/user/benitezj/ggHbb/bits/' + options.tag + '/' + options.sample
    print eosOutDir

    #### Log directory
    logDir = '/uscms_data/d3/benitezj/production/ggHbb/bits/' + options.tag + '/' + options.sample
    print logDir

    #### choose algorithm options
    if isMC:
        execOptions = "%s Output.root --passSumEntries 5:Events -a 6:subjob_i -a 7:1 -a 2:mc -a 3:none -a 4:%f  -n 8000 -q 2nw4cores "%(options.executable,xsec)
    else:
        execOptions = "%s Output.root -a 5:1 -a 6:subjob_i -a 7:1 -a 2:data -a 3:%s -a 4:%f -n 8000 -q 1nd "%(options.executable,json,xsec)


    #### set the command to be executed
    command="python BaconAnalyzer/Analyzer/production/baconCondor.py %s  --njobs-per-file 1 --list 1:%s --outdir %s --eosoutdir %s"%(execOptions,inputFiles,logDir,eosOutDir)

    if options.sub:
        command=command+" --monitor sub"

    print command  

    #### run
    if not options.test:
        if not options.sub:
            os.system('mkdir -p /eos/uscms/%s'%(eosOutDir))
            os.system('mkdir -p %s'%(logDir))
        os.system(command)
