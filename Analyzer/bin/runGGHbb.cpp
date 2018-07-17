//================================================================================================
//________________________________________________________________________________________________

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/ElectronLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PhotonLoader.hh"
#include "../include/TauLoader.hh"
#include "../include/JetLoader.hh"
#include "../include/VJetLoader.hh"
#include "../include/RunLumiRangeMap.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <TError.h>
#include <string>
#include <iostream>
#include <fstream>


const float PTCUT = 400;



#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int parseLine(char* line){
  // This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}

int getValue(){ //Note: this value is in KB!
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}



TChain* loadChain(std::string inputFile, unsigned int firstfile, unsigned int lastfile) { 

  std::ifstream infile(inputFile.c_str());
  if(! infile.is_open()){
    std::cout<<"failed to open input file list:"<< inputFile.c_str()<<std::endl;
    return NULL;
  }

  TChain * lTree=new TChain("Events"); 
  std::string line;
  unsigned int counterInput=0;
  while (std::getline(infile, line)){
    counterInput++;
    
    if(counterInput>=firstfile && counterInput<=lastfile){
      
      //cout<<line<<endl;
      TFile *lFile = TFile::Open(line.c_str());
      if(lFile->IsZombie()){
	std::cout<<"failed to open root file "<<line.c_str()<<std::endl;
	return NULL;
      }

      TTree *lT = (TTree*) lFile->FindObjectAny("Events");
      if(lT==NULL){
	std::cout<<"failed to get Tree from "<<line.c_str()<<std::endl;
	return NULL;
      }
      
      std::cout<<line.c_str()<<std::endl;
      lTree->Add(line.c_str());
    }
  }

  return lTree;
}


// === MAIN =======================================================================================================
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");//why is this here?
  const int      firstfile       = atoi(argv[1]);
  const int      lastfile        = atoi(argv[2]);
  const std::string lName        = argv[3]; //input file (w.r.t. src/  dir)

  std::string cmssw=getenv("CMSSW_BASE");

  bool isData=false;
  if(lName.find("Run20")!=std::string::npos) isData = true;

  ///set the json for data
  RunLumiRangeMap fRangeMap;
  if(isData ){
    std::string lJson=cmssw;
    lJson.append("/src/BaconAnalyzer/Analyzer/data/"); 
    if(lName.find("Run2016")!=std::string::npos) 
      lJson.append("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt");
    else if(lName.find("Run2017")!=std::string::npos) 
      lJson.append("Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt");
    else {
      std::cout<<" no json"<<std::endl;
      return 0;
    }
    std::cout << "json " << lJson << std::endl;
    fRangeMap.AddJSONFile(lJson.c_str());
  }
  


  //////get the input
  std::string linput=cmssw;
  linput.append("/src/");
  linput.append(lName);//lName must be w.r.t. src/
  TChain *lTree = loadChain(lName,firstfile,lastfile); 


  // Object Processors
  GenLoader       fGen(lTree); 
  EvtLoader       fEvt(lTree,lName); 
  MuonLoader      fMuon(lTree); 
  ElectronLoader  fElectron(lTree); 
  TauLoader  fTau(lTree);          
  PhotonLoader    fPhoton(lTree); 
  JetLoader     fJet4(lTree, isData);  
  VJetLoader      fVJet8(lTree,"AK8Puppi","AddAK8Puppi",3, isData);
  //VJetLoader      fVJet15(lTree,"CA15Puppi","AddCA15Puppi",3, isData);
  

  // Output file and tree
  TFile lFile("Output.root","RECREATE");
  gROOT->cd();

  //Setup histograms containing total number of processed events (for normalization)
  TH1F NEvents("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F SumWeights("SumWeights", "SumWeights", 1, 0.5, 1.5);
   

  // Setup Tree
  TTree lOut("Events","Events");
  lOut.SetDirectory(&lFile);
  fEvt.setupTree(&lOut); 
  fVJet8.setupTree(&lOut,"AK8Puppijet"); 
  fVJet8.setupTreeZprime(&lOut,"AK8Puppijet");
  //fVJet15.setupTree(&lOut,"CA15Puppijet");
  //fVJet15.setupTreeZprime(&lOut,"CA15Puppijet");
  fJet4.setupTree(&lOut,"AK4Puppijet");
  fMuon.setupTree(&lOut); 
  fElectron.setupTree(&lOut); 
  fTau.setupTree(&lOut); 
  fPhoton.setupTree(&lOut); 
  if(!isData)fGen.setupTree (&lOut,1.);

  //Add the triggers we want
  fEvt.addTrigger("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*");
  fEvt.addTrigger("HLT_PFHT800_v*");//pre-scaled in 2017
  fEvt.addTrigger("HLT_PFHT900_v*");//pre-scaled in 2017
  fEvt.addTrigger("HLT_PFHT1050_v*")  ; 
  fEvt.addTrigger("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*");
  fEvt.addTrigger("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*");
  fEvt.addTrigger("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v*");
  fEvt.addTrigger("HLT_PFJet450_v*");
  fEvt.addTrigger("HLT_PFJet500_v*");

  fEvt.addTrigger("HLT_IsoMu24_v*"); // W(munu)H(bb)
  fEvt.addTrigger("HLT_IsoTkMu24_v*"); // W(munu)H(bb)
  fEvt.addTrigger("HLT_Mu50_v*");
  fEvt.addTrigger("HLT_TkMu50_v*");
    
  fEvt.addTrigger("HLT_Ele27_WPTight_Gsf_v*"); // W(enu)H(bb)
  fEvt.addTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //Z(ee)H(bb)
  fEvt.addTrigger("HLT_Ele45_WPLoose_v*");
  fEvt.addTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"); 

  fEvt.addTrigger("HLT_PFMET110_PFMHT110_IDTight_v*"); // Z(nunu)H(bb) 
  fEvt.addTrigger("HLT_PFMET120_PFMHT120_IDTight_v*"); // Z(nunu)H(bb)
  fEvt.addTrigger("HLT_PFMET170_NoiseCleaned_v*"); // Z(nunu)H(bb)
  fEvt.addTrigger("HLT_PFMET170_HBHECleaned_v*"); // Z(nunu)H(bb)
  fEvt.addTrigger("HLT_PFMET170_HBHE_BeamHaloCleaned_v*"); // Z(nunu)H(bb)

  //ggH(bb)
  fEvt.addTrigger("HLT_AK8PFJet360_TrimMass30_v*");//pre-scaled in 2017
  fEvt.addTrigger("HLT_AK8PFJet380_TrimMass30_v*");//pre-scaled in 2017
  fEvt.addTrigger("HLT_AK8PFJet400_TrimMass30_v*");
  fEvt.addTrigger("HLT_AK8PFJet420_TrimMass30_v*");
  fEvt.addTrigger("HLT_AK8PFHT800_TrimMass50_v*");
  fEvt.addTrigger("HLT_AK8PFHT850_TrimMass50_v*");
  fEvt.addTrigger("HLT_AK8PFHT900_TrimMass50_v*");
  fEvt.addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*");
  fEvt.addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*");
  fEvt.addTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
  fEvt.addTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
  fEvt.addTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
  fEvt.addTrigger("HLT_CaloJet500_NoJetID_v*");
  fEvt.addTrigger("HLT_CaloJet550_NoJetId_v*");
  fEvt.addTrigger("HLT_AK8PFJet500_v*");
  fEvt.addTrigger("HLT_AK8PFJet550_v*");



  // Loop over events i0 = iEvent
  int neventstest = 0;
  unsigned int totalEvents = lTree->GetEntries();
  for(unsigned int i0 = 0; i0 < totalEvents; i0++) {
    //if(i0<115000) continue;
    if (i0%1000 == 0) std::cout << i0 <<"/"<<totalEvents<< " events processed , "<<neventstest<<" , mem="<< getValue()<< std::endl;

    ///load all branches
    fEvt.load(i0);
    if(!isData) fGen.load(i0);
    fMuon.load(i0);
    fElectron.load(i0);
    fTau.load(i0);
    fPhoton.load(i0);
    //fVJet15.load(i0);
    fVJet8.load(i0);
    fJet4.load(i0); 

    lTree->GetEntry(i0);
 

    // Check GenInfo
    float lWeight = 1;
    unsigned int passJson = 0;
    if(!isData){
  
      if(lName.find("TTJets")!=std::string::npos || lName.find("TT_")!=std::string::npos || lName.find("TTTo")!=std::string::npos){
	float ttbarPtWeight = fGen.computeTTbarCorr();
	//fEvt.fevtWeight *= ttbarPtWeight; //not needed because gets added in fGen.fWeight
	fGen.fWeight *= ttbarPtWeight;
	fGen.saveTTbarType();
      }
      
      
      lWeight = fGen.fWeight;
      passJson = 1;
    }else{
      RunLumiRangeMap::RunLumiPairType lRunLumi(fEvt.fRun,fEvt.fLumi);
      if(fRangeMap.HasRunLumi(lRunLumi)) passJson = 1;
    }


    ///Save the initial sum of weights
    NEvents.SetBinContent(1, NEvents.GetBinContent(1)+lWeight);
    SumWeights.Fill(1.0, lWeight);


    /// apply json otherwise some jetEC give errors creating large log files
    if(!passJson) continue;

 
    // Triggerbits
    unsigned int trigbits=1;   
    if(isData){
      if(
	 fEvt.passTrigger("HLT_AK8PFJet400_TrimMass30_v*") ||
	 fEvt.passTrigger("HLT_AK8PFJet420_TrimMass30_v*") ||
	 fEvt.passTrigger("HLT_AK8PFHT800_TrimMass50_v*") ||
	 fEvt.passTrigger("HLT_AK8PFHT850_TrimMass50_v*") ||
	 fEvt.passTrigger("HLT_AK8PFHT900_TrimMass50_v*") ||
	 fEvt.passTrigger("HLT_PFJet500_v*") ||
	 fEvt.passTrigger("HLT_PFHT1050_v*") ||
	 fEvt.passTrigger("HLT_AK8PFJet500_v*") ||
	 fEvt.passTrigger("HLT_AK8PFJet550_v*") ||
	 fEvt.passTrigger("HLT_CaloJet500_NoJetID_v*") ||
	 fEvt.passTrigger("HLT_CaloJet550_NoJetId_v*") ||
	 fEvt.passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*") ||
	 fEvt.passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*")
	 )  trigbits = trigbits | 2;  // hadronic signal region
      if( fEvt.passTrigger("HLT_Mu50_v*") ||
	  fEvt.passTrigger("HLT_TkMu50_v*")
	  ) trigbits = trigbits | 4; // single muon control region
      if( fEvt.passTrigger("HLT_Ele45_WPLoose_v*") ||
	  fEvt.passTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")
	  ) trigbits = trigbits | 8; // single electron control region 
      if(
	 fEvt.passTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt.passTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt.passTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") 
	 )   trigbits = trigbits | 16; // AK4maxDeta + doubleBTag  
    }


    fEvt.fillEvent(trigbits,lWeight,passJson);
    


    // TTbar, EWK and kFactor correction 
    // Note: computeCorr sets the  fEvt.fkfactor  which is set to 1 in fEvt.fillEvent() above, so keep this order. 
    // Does not change the eventweight "scale1fb" which is set to  lWeight above)
    if(lName.find("ZJets")!=std::string::npos || lName.find("DYJets")!=std::string::npos){
      fGen.findBoson(23,0);
      if(fGen.fBosonPt>0)      fEvt.computeCorr(fGen.fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
    }
    if(lName.find("WJets")!=std::string::npos){
      fGen.findBoson(24,1);
      if(fGen.fBosonPt>0)      fEvt.computeCorr(fGen.fBosonPt,"WJets_012j_NLO/nominal","WJets_LO/inv_pt","EWKcorr/W","WJets_012j_NLO");
    }
    if(lName.find("ZPrime")!=std::string::npos || lName.find("VectorDiJet")!=std::string::npos){
      fGen.findBoson(55,1);
      if(fGen.fBosonPt>0)      fEvt.computeCorr(fGen.fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
    }
    if(lName.find("Spin0")!=std::string::npos){
      fGen.findBoson(55,1);
    }
    if(lName.find("HToBB")!=std::string::npos || lName.find("HTobb")!=std::string::npos){
      fGen.findBoson(25,1);
    }



    // Primary vertex requirement
    if(!fEvt.PV()) continue;


    // Objects
    gErrorIgnoreLevel=kError;
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    fMuon.selectMuons(cleaningMuons,fEvt.fMet,fEvt.fMetPhi);
    fElectron.selectElectrons(fEvt.fRho,fEvt.fMet,cleaningElectrons);
    fPhoton.selectPhotons(fEvt.fRho,cleaningElectrons,cleaningPhotons);

    fTau.selectTaus(cleaningElectrons, cleaningMuons);

    // CA15Puppi Jets
    //fVJet15.selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,1.5,fEvt.fRho,fEvt.fRun);
    //if(fVJet15.selectedVJets.size()>0) fEvt.fselectBits =  fEvt.fselectBits | 4;

    // AK8Puppi Jets    
    fVJet8.selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt.fRho,fEvt.fRun);
    if(fVJet8.fLooseVJets.size()>0) fEvt.fselectBits =  fEvt.fselectBits | 2;

    // Match leading AK8 Puppi jet with CA15 Puppi jet within dR = 0.4 (to get pT ratio)
    //if(fVJet8.selectedVJets.size()>0) fVJet8 .matchJet15(fVJet15.selectedVJets,fVJet8.selectedVJets[0],0.4);
    

    
    // Select at least one AK8 or one CA15 jet
    if(!(fEvt.fselectBits & 2))continue;
    if(fVJet8.fLooseVJets[0]->pt < PTCUT)continue;
    //if((!(fEvt.fselectBits & 2)||(fVJet8.fLooseVJets[0]->pt < PTCUT)) && (!(fEvt.fselectBits & 4)||(fVJet15.fLooseVJets[0]->pt < PTCUT))) continue;


    TLorentzVector vJet8;
    if(fVJet8.fLooseVJets.size()>0) 
      vJet8.SetPtEtaPhiM(fVJet8.fLooseVJets[0]->pt,fVJet8.fLooseVJets[0]->eta,fVJet8.fLooseVJets[0]->phi,fVJet8.fLooseVJets[0]->mass);

    std::vector<TLorentzVector> vJet8List;
    vJet8List.push_back(vJet8);


    // TLorentzVector vJet15;
    // if(fVJet15.fLooseVJets.size()>0) 
    //   vJet15.SetPtEtaPhiM(fVJet15.fLooseVJets[0]->pt,fVJet15.fLooseVJets[0]->eta,fVJet15.fLooseVJets[0]->phi,fVJet15.fLooseVJets[0]->mass);

    // std::vector<TLorentzVector> vJet15List;
    // vJet15List.push_back(vJet15);
    // if(fVJet8.fLooseVJets.size()>0) fVJet8.matchJet15(vJet15List,vJet8,0.4);
    

    fJet4.selectJets(cleaningElectrons,cleaningMuons,cleaningPhotons,vJet8List,fEvt.fRho,fEvt.fRun);


    // truth match
    if(lName.find("ZJets")!=std::string::npos || lName.find("DYJets")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,23);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,23);
    }
    if(lName.find("WJets")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,24);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,24);
    }
    if(lName.find("ZPrime")!=std::string::npos || lName.find("VectorDiJet")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,10031);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,10031);
    }
    if(lName.find("Spin0")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,55);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,55);
    }
    if(lName.find("TTJets")!=std::string::npos || lName.find("TT_")!=std::string::npos || lName.find("TTTo")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,624);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,624);
    }    
    if(lName.find("ST_")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,624);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,624);
    }
    if(lName.find("HToBB")!=std::string::npos || lName.find("HTobb")!=std::string::npos){
      if(fVJet8.fLooseVJets.size()>0) fVJet8.fisHadronicV = fGen.ismatchedJet(vJet8,0.8,fVJet8.fvMatching,fVJet8.fvSize,25);
      //if(fVJet15.fLooseVJets.size()>0)  fVJet15.fisHadronicV = fGen.ismatchedJet(vJet15,1.5,fVJet15.fvMatching,fVJet15.fvSize,25);
    }



    lOut.Fill();
    neventstest++;
  }


  lFile.cd();
  lOut.Write();  
  NEvents.Write();
  SumWeights.Write();
  lFile.ls();
  lFile.Close();


  //delete lTree;

  std::cout <<"SelectedEvents : "<< neventstest << std::endl;
}
