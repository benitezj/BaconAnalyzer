
## Set these options before running 
TEST=0  ## 0= all samples,  1= only run one signal sample
SUBMIT=$1 ## 0 : only print command, 1 : create the job scripts, 2 : condor submission, 3: check log files, 4: check log file and resubmit 


INPUTPROD=BaconAnalyzer/Analyzer/lists/production14

LOGDIR=/uscms_data/d3/benitezj/production/ggHbb/bits/ggHbits-b14-01  

EOSOUTPUTDIR=/store/user/benitezj/ggHbb/bits/ggHbits-b14-01 

echo 'input prod= ' $INPUTPROD
echo 'local log dir = ' $LOGDIR
echo 'eos output dir = ' $EOSOUTPUTDIR



##############################################
### need to tar the CMSSW to submit condor job
if [ "${SUBMIT}" == "2" ]; then
    /bin/rm -f ${CMSSW_BASE}.tar
    /bin/tar -cvf ${CMSSW_BASE}.tar -C $CMSSW_BASE/../ $CMSSW_VERSION
    if [ ! -f ${CMSSW_BASE}.tar ]; then
	echo "CMSSW tar not created"
	return 0
    fi
fi

### job submision function needed below
submit()
{
    eval merge="$1"
    eval sample="$2"

    echo "$sample:"

    local TOTFILES=`cat ${CMSSW_BASE}/src/${INPUTPROD}/${sample}.txt | wc -l`
    local NJOBS=`echo "$TOTFILES $merge" | awk -F' ' '{print int($1/$2)}'`
    local   REM=`echo "$TOTFILES $merge" | awk -F' ' '{print $1%$2}'`
    if [ "$REM" != "0" ]; then 
	NJOBS=`echo $NJOBS | awk '{print $1+1}'`
    fi

    ## condor submission
    local COUNTER=0
    for f in `cat ${CMSSW_BASE}/src/${INPUTPROD}/${sample}.txt`; do 
	
	##truncate loop at the number of desired jobs
	if [ "$COUNTER" == "$NJOBS" ];then return; fi

	###define job parameters
	local FIRST=`echo "$merge $COUNTER" | awk -F' ' '{print $1*$2 + 1}'`
	local LAST=`echo "$merge $COUNTER" | awk -F' ' '{print $1*($2+1) }'`
	local command="runGGHbb $FIRST $LAST ${INPUTPROD}/${sample}.txt"
	local outfile=${sample}_${COUNTER}
	local condorsub="/usr/bin/condor_submit $LOGDIR/${outfile}.sub"

	### only print the command
	if [ "${SUBMIT}" == "0" ]; then
	    echo $command
	fi

	#####create the job and submit
	if [ "${SUBMIT}" == "1" ]; then

	############
	## clean out the output
	#########
	    /bin/rm -f $LOGDIR/${outfile}.sh
	    /bin/rm -f $LOGDIR/${outfile}.sub

	######################
	### create the execution script
	#######################
	    touch $LOGDIR/${outfile}.sh
	    echo "pwd"  >> $LOGDIR/${outfile}.sh
	    echo "mount"  >> $LOGDIR/${outfile}.sh
            echo "/bin/tar -xf ${CMSSW_VERSION}.tar"  >> $LOGDIR/${outfile}.sh
	    echo "ls ."  >> $LOGDIR/${outfile}.sh
	    echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"  >> $LOGDIR/${outfile}.sh
            echo "export SCRAM_ARCH=slc6_amd64_gcc530"  >> $LOGDIR/${outfile}.sh 
	    echo "cd ./${CMSSW_VERSION}/src"  >> $LOGDIR/${outfile}.sh
	    echo "scramv1 b ProjectRename "  >> $LOGDIR/${outfile}.sh
	    echo "eval \`scramv1 runtime -sh\` "  >> $LOGDIR/${outfile}.sh
	    echo "env"  >> $LOGDIR/${outfile}.sh
	    echo "${command}" >> $LOGDIR/${outfile}.sh 
	    echo "xrdcp Output.root root://cmseos.fnal.gov/${EOSOUTPUTDIR}/${outfile}.root" >> $LOGDIR/${outfile}.sh 
	    
	################
	### create condor jdl
	################
	    touch $LOGDIR/${outfile}.sub
	    echo "Universe   = vanilla" >> $LOGDIR/${outfile}.sub 
	    echo "Executable = /bin/bash" >> $LOGDIR/${outfile}.sub 
	    echo "Log        = $LOGDIR/${outfile}.log" >> $LOGDIR/${outfile}.sub
	    echo "Output     = $LOGDIR/${outfile}.log" >> $LOGDIR/${outfile}.sub
	    echo "Error      = $LOGDIR/${outfile}.log" >> $LOGDIR/${outfile}.sub
	    echo "Arguments  = ${outfile}.sh" >> $LOGDIR/${outfile}.sub
	    echo "Should_Transfer_Files = YES" >> $LOGDIR/${outfile}.sub
	    echo "WhenToTransferOutput = ON_EXIT" >> $LOGDIR/${outfile}.sub
	    echo "Transfer_Input_Files = ${LOGDIR}/${outfile}.sh, ${CMSSW_BASE}.tar" >> $LOGDIR/${outfile}.sub
	    echo "Queue" >> $LOGDIR/${outfile}.sub
	   
	fi  

	#### condor submission
	if [ "${SUBMIT}" == "2" ] ; then
	    /bin/rm -f $LOGDIR/${outfile}.log
	    /bin/rm -f /eos/uscms/${EOSOUTPUTDIR}/${outfile}.root
	    echo $condorsub
	    `${condorsub}`
	fi
	    
	####check if the job was successfull
	if [ "${SUBMIT}" == "3" ] || [ "${SUBMIT}" == "4" ] ; then
	    local selected=`cat $LOGDIR/${outfile}.log | grep Selected`
	    if [ ! -f  $LOGDIR/${outfile}.log ] || [ "{selected}" == "" ] || [ ! -f /eos/uscms/${EOSOUTPUTDIR}/${outfile}.root ] ; then
		echo "$LOGDIR/${outfile}.log"
		if [ "${SUBMIT}" == "4" ] ; then
		    ## resubmit failed job
		    /bin/rm -f $LOGDIR/${outfile}.log
		    /bin/rm -f /eos/uscms/${EOSOUTPUTDIR}/${outfile}.root
		    `${condorsub}`
		fi
	    fi
	fi


	COUNTER=`echo $COUNTER | awk '{print $1+1}'`
     done


}


#################################
### submit the jobs

#submit 100 GluGluHToBB_M125_13TeV_powheg_pythia8
#if [ "$TEST" == "1" ]; then return 1 ; fi
#submit 100 VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix
#submit 100 WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8
#submit 100 WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8
#submit 100 ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8
#submit 100 ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8
#submit 100 ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_herwigpp
#submit 100 ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8
#submit 100 ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8
#
#submit 100 DYJetsToLL_M_50_HT_100to200_TuneCP5_13TeV
#submit 100 DYJetsToLL_M_50_HT_200to400_TuneCP5_13TeV
#submit 100 DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV
#submit 100 DYJetsToLL_M_50_HT_600to800_TuneCP5_13TeV
#submit 100 DYJetsToLL_M_50_HT_800to1200_TuneCP5_13TeV
#
#submit 100 WJetsToLNu_HT_100To200_TuneCP5_13TeV
#submit 100 WJetsToLNu_HT_1200To2500_TuneCP5_13TeV
#submit 100 WJetsToLNu_HT_200To400_TuneCP5_13TeV
#submit 100 WJetsToLNu_HT_400To600_TuneCP5_13TeV
#submit 100 WJetsToLNu_HT_600To800_TuneCP5_13TeV
#submit 100 WJetsToLNu_HT_800To1200_TuneCP5_13TeV
#
#submit 100 WW_TuneCP5_13TeV_pythia8
#submit 100 WZ_TuneCP5_13TeV_pythia8
#submit 100 ZZ_TuneCP5_13TeV_pythia8
#
#submit 100 ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8
#submit 100 ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8
#submit 100 ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8
#submit 100 ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8
#
#submit 100 TTJets_TuneCP5_13TeV_amcatnloFXFX_pythia8
#submit 100 TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8
#submit 100 TTToSemiLeptonic_WspTgt150_TuneCUETP8M2T4_13TeV_powheg_pythia8
#
#submit 10 QCD_HT1000to1500_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT100to200_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT1500to2000_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT2000toInf_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT200to300_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT300to500_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT500to700_TuneCP5_13TeV_madgraph_pythia8
#submit 10 QCD_HT700to1000_TuneCP5_13TeV_madgraph_pythia8

submit 50 JetHTRun2017B_17Nov2017_v1
submit 50 JetHTRun2017C_17Nov2017_v1
submit 50 JetHTRun2017D_17Nov2017_v1
submit 50 JetHTRun2017E_17Nov2017_v1
submit 50 JetHTRun2017F_17Nov2017_v1

submit 50 SingleMuonRun2017B_17Nov2017_v1
submit 50 SingleMuonRun2017C_17Nov2017_v1
submit 50 SingleMuonRun2017D_17Nov2017_v1
submit 50 SingleMuonRun2017E_17Nov2017_v1
submit 50 SingleMuonRun2017F_17Nov2017_v1


#### show the jobs in the queue
if [ "${SUBMIT}" == "2" ]; then
    /usr/bin/condor_q
fi
