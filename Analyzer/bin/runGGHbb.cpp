//================================================================================================
//
//
// Input arguments
//   argv[1] => lName = input bacon file name
//   argv[2] => lOption = dataset type: "mc", "data"
//   argv[3] => lJSON = JSON file for run-lumi filtering of data, specify "none" for MC or no filtering
//   argv[4] => lXS = cross section (pb), ignored for data 
//   argv[5] => weight = total weight, ignored for data
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


// Object Processors
GenLoader       *fGen        = 0; 
EvtLoader       *fEvt        = 0; 
MuonLoader      *fMuon       = 0; 
ElectronLoader  *fElectron   = 0; 
TauLoader       *fTau        = 0; 
PhotonLoader    *fPhoton     = 0; 
JetLoader       *fJet4       = 0;
VJetLoader      *fVJet8      = 0;
VJetLoader      *fVJet15     = 0;
RunLumiRangeMap *fRangeMap   = 0; 

const float PTCUT = 400;


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


// For Json 
bool passEvent(unsigned int iRun,unsigned int iLumi) { 
  RunLumiRangeMap::RunLumiPairType lRunLumi(iRun,iLumi);
  return fRangeMap->HasRunLumi(lRunLumi);
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
  fRangeMap = new RunLumiRangeMap();
  if(isData ){
    std::string lJson=cmssw;
    lJson.append("/src/BaconAnalyzer/Analyzer/data/"); 
    if(lName.find("Run2016")!=std::string::npos) lJson.append("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt");
    else if(lName.find("Run2017")!=std::string::npos) lJson.append("Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt");
    else {std::cout<<" no json"<<std::endl; return 0;}
    std::cout << "json " << lJson << std::endl;
    fRangeMap->AddJSONFile(lJson.c_str());
  }
  


  //////get the input
  std::string linput=cmssw;
  linput.append("/src/");
  linput.append(lName);//lName must be w.r.t. src/
  TChain *lTree = loadChain(lName,firstfile,lastfile); 


  // Declare Readers 
  fEvt       = new EvtLoader     (lTree,lName);    
  fMuon      = new MuonLoader    (lTree);          
  fElectron  = new ElectronLoader(lTree);          
  fTau       = new TauLoader     (lTree);          
  fPhoton    = new PhotonLoader  (lTree);          
  fJet4      = new JetLoader     (lTree, isData);  
  fVJet8     = new VJetLoader    (lTree,"AK8Puppi","AddAK8Puppi",3, isData);  
  fVJet15    = new VJetLoader    (lTree,"CA15Puppi","AddCA15Puppi",3, isData);
  fGen       = new GenLoader     (lTree);                 


  // Output file and tree
  TFile *lFile = TFile::Open("Output.root","RECREATE");
  gROOT->cd();

  //Setup histograms containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 0.5, 1.5);
  TH1F *SumWeights = new TH1F("SumWeights", "SumWeights", 1, 0.5, 1.5);
  TH1F *SumScaleWeights = new TH1F("SumScaleWeights", "SumScaleWeights", 6, -0.5, 5.5);
  

  // Setup Tree
  TTree *lOut  = new TTree("Events","Events");
  lOut->SetDirectory(lFile);
  fEvt      ->setupTree      (lOut); 
  fVJet8    ->setupTree      (lOut,"AK8Puppijet"); 
  fVJet8    ->setupTreeZprime(lOut,"AK8Puppijet");
  fVJet15   ->setupTree      (lOut,"CA15Puppijet");
  fVJet15   ->setupTreeZprime(lOut,"CA15Puppijet");
  fJet4     ->setupTree      (lOut,"AK4Puppijet");
  fMuon     ->setupTree      (lOut); 
  fElectron ->setupTree      (lOut); 
  fTau      ->setupTree      (lOut); 
  fPhoton   ->setupTree      (lOut); 
  if(isData==false) fGen ->setupTree (lOut,1.);


  // Loop over events i0 = iEvent
  int neventstest = 0;
  unsigned int totalEvents = lTree->GetEntries();
  for(unsigned int i0 = 0; i0 < totalEvents; i0++) {
    if (i0%1000 == 0) std::cout << i0 <<"/"<<totalEvents<< " events processed " << std::endl;

    ///load all branches
    fEvt->load(i0);
    if(isData==false) fGen->load(i0);
    fMuon     ->load(i0);
    fElectron ->load(i0);
    fTau      ->load(i0);
    fPhoton   ->load(i0);
    fVJet15   ->load(i0);
    fVJet8    ->load(i0);
    fJet4     ->load(i0);
 
    lTree->GetEntry(i0);
    
    // Check GenInfo
    float lWeight = 1;
    unsigned int passJson = 0;
    if(isData==false){
      //fGen->load(i0);
      lWeight = fGen->fWeight;
      passJson = 1;
    }
    else{
      if(passEvent(fEvt->fRun,fEvt->fLumi)) { passJson = 1;}
    }

    
    NEvents->SetBinContent(1, NEvents->GetBinContent(1)+lWeight);
    SumWeights->Fill(1.0, lWeight);

    // Primary vertex requirement
    if(!fEvt->PV()) continue;
    
    // Triggerbits
    unsigned int trigbits=1;   
    if(isData){
      if(
	 fEvt ->passTrigger("HLT_AK8PFJet400_TrimMass30_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet420_TrimMass30_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT800_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT850_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFHT900_TrimMass50_v*") ||
	 fEvt ->passTrigger("HLT_PFJet500_v*") ||
	 fEvt ->passTrigger("HLT_PFHT1050_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet500_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet550_v*") ||
	 fEvt ->passTrigger("HLT_CaloJet500_NoJetID_v*") ||
	 fEvt ->passTrigger("HLT_CaloJet550_NoJetId_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*") ||
	 fEvt ->passTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*")
	 )  trigbits = trigbits | 2;  // hadronic signal region
      if( fEvt ->passTrigger("HLT_Mu50_v*") ||
	  fEvt ->passTrigger("HLT_TkMu50_v*")
	  ) trigbits = trigbits | 4; // single muon control region
      if( fEvt ->passTrigger("HLT_Ele45_WPLoose_v*") ||
	  fEvt ->passTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*")
	  ) trigbits = trigbits | 8; // single electron control region 
      if(
	 fEvt ->passTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt ->passTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") ||
	 fEvt ->passTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*") 
	)   trigbits = trigbits | 16; // AK4maxDeta + doubleBTag  
    }
    // More trigger bits

    //Add the triggers we want
    fEvt ->addTrigger("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*");
    fEvt ->addTrigger("HLT_PFHT800_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_PFHT900_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_PFHT1050_v*")  ; 
    fEvt ->addTrigger("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*");
    fEvt ->addTrigger("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*");
    fEvt ->addTrigger("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v*");
    fEvt ->addTrigger("HLT_PFJet450_v*");
    fEvt ->addTrigger("HLT_PFJet500_v*");

    fEvt ->addTrigger("HLT_IsoMu24_v*"); // W(munu)H(bb)
    fEvt ->addTrigger("HLT_IsoTkMu24_v*"); // W(munu)H(bb)
    fEvt ->addTrigger("HLT_Mu50_v*");
    fEvt ->addTrigger("HLT_TkMu50_v*");
    
    fEvt ->addTrigger("HLT_Ele27_WPTight_Gsf_v*"); // W(enu)H(bb)
    fEvt ->addTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"); //Z(ee)H(bb)
    fEvt ->addTrigger("HLT_Ele45_WPLoose_v*");
    fEvt ->addTrigger("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"); 

    fEvt ->addTrigger("HLT_PFMET110_PFMHT110_IDTight_v*"); // Z(nunu)H(bb) 
    fEvt ->addTrigger("HLT_PFMET120_PFMHT120_IDTight_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_NoiseCleaned_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_HBHECleaned_v*"); // Z(nunu)H(bb)
    fEvt ->addTrigger("HLT_PFMET170_HBHE_BeamHaloCleaned_v*"); // Z(nunu)H(bb)

    //ggH(bb)
    fEvt ->addTrigger("HLT_AK8PFJet360_TrimMass30_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_AK8PFJet380_TrimMass30_v*");//pre-scaled in 2017
    fEvt ->addTrigger("HLT_AK8PFJet400_TrimMass30_v*");
    fEvt ->addTrigger("HLT_AK8PFJet420_TrimMass30_v*");
    fEvt ->addTrigger("HLT_AK8PFHT800_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFHT850_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFHT900_TrimMass50_v*");
    fEvt ->addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p17_v*");
    fEvt ->addTrigger("HLT_AK8PFJet330_PFAK8BTagCSV_p1_v*");
    fEvt ->addTrigger("HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagCSV_p33_v*");
    fEvt ->addTrigger("HLT_CaloJet500_NoJetID_v*");
    fEvt ->addTrigger("HLT_CaloJet550_NoJetId_v*");
    fEvt ->addTrigger("HLT_AK8PFJet500_v*");
    fEvt ->addTrigger("HLT_AK8PFJet550_v*");

    fEvt      ->fillEvent(trigbits,lWeight,passJson);
    
    // Objects
    gErrorIgnoreLevel=kError;
    std::vector<TLorentzVector> cleaningMuons, cleaningElectrons, cleaningPhotons; 
    fMuon     ->selectMuons(cleaningMuons,fEvt->fMet,fEvt->fMetPhi);
    fElectron ->selectElectrons(fEvt->fRho,fEvt->fMet,cleaningElectrons);
    fTau      ->selectTaus(cleaningElectrons, cleaningMuons);
    fPhoton   ->selectPhotons(fEvt->fRho,cleaningElectrons,cleaningPhotons);
        
    // CA15Puppi Jets
    fVJet15   ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,1.5,fEvt->fRho,fEvt->fRun);
    if(fVJet15->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 4;
      
    // AK8Puppi Jets    
    fVJet8    ->selectVJets(cleaningElectrons,cleaningMuons,cleaningPhotons,0.8,fEvt->fRho,fEvt->fRun);
    if(fVJet8->selectedVJets.size()>0) fEvt->fselectBits =  fEvt->fselectBits | 2;

    // Match leading AK8 Puppi jet with CA15 Puppi jet within dR = 0.4 (to get pT ratio)
    if(fVJet8->selectedVJets.size()>0) fVJet8 ->matchJet15(fVJet15->selectedVJets,fVJet8->selectedVJets[0],0.4);
    
    // AK4Puppi Jets
    fJet4     ->selectJets(cleaningElectrons,cleaningMuons,cleaningPhotons,fVJet8->selectedVJets,fEvt->fRho,fEvt->fRun);

    // Select at least one AK8 or one CA15 jet
    if(!(fEvt->fselectBits & 2) || !(fEvt->fselectBits & 4)) continue;
    if((fVJet8->selectedVJets[0].Pt() < PTCUT) || (fVJet15->selectedVJets[0].Pt() < PTCUT)) continue;

    // TTbar, EWK and kFactor correction
    if(lName.find("ZJets")!=std::string::npos || lName.find("DYJets")!=std::string::npos){
      fGen->findBoson(23,0);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0)  fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,23);
      if(fVJet15->selectedVJets.size()>0)  fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,23);
    }
    if(lName.find("WJets")!=std::string::npos){
      fGen->findBoson(24,1);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"WJets_012j_NLO/nominal","WJets_LO/inv_pt","EWKcorr/W","WJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0)  fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,24);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,24);
    }
    if(lName.find("ZPrime")!=std::string::npos || lName.find("VectorDiJet")!=std::string::npos){
      // 55 for newer, public samples (to be consistent with DM samples) 
      //fGen->findBoson(10031,0);
      fGen->findBoson(55,1);
      if(fGen->fBosonPt>0)      fEvt->computeCorr(fGen->fBosonPt,"ZJets_012j_NLO/nominal","ZJets_LO/inv_pt","EWKcorr/Z","ZJets_012j_NLO");
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,10031);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,10031);
    }
    if(lName.find("Spin0")!=std::string::npos){
      // 9900032 for older, private samples
      // 55 for newer, public samples (to be consistent with DM samples)
      fGen->findBoson(55,1);
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,55);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,55);
    }
    if(lName.find("TTJets")!=std::string::npos || lName.find("TT_")!=std::string::npos || lName.find("TTTo")!=std::string::npos){
      float ttbarPtWeight = fGen->computeTTbarCorr();
      fEvt->fevtWeight *= ttbarPtWeight;
      fGen->fWeight *= ttbarPtWeight;
      fGen->saveTTbarType();
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,624);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,624);
    }    
    if(lName.find("ST_")!=std::string::npos){
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,624);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,624);
    }
    if(lName.find("HToBB")!=std::string::npos || lName.find("HTobb")!=std::string::npos){
      fGen->findBoson(25,1);
      if(fVJet8->selectedVJets.size()>0) fVJet8->fisHadronicV = fGen->ismatchedJet(fVJet8->selectedVJets[0],0.8,fVJet8->fvMatching,fVJet8->fvSize,25);
      if(fVJet15->selectedVJets.size()>0) fVJet15->fisHadronicV = fGen->ismatchedJet(fVJet15->selectedVJets[0],1.5,fVJet15->fvMatching,fVJet15->fvSize,25);
    }

    lOut->Fill();
    neventstest++;
  }


  lFile->cd();
  lOut->Write();  
  NEvents->Write();
  SumWeights->Write();
  SumScaleWeights->Write();
  lFile->ls();
  lFile->Close();


  std::cout <<"SelectedEvents : "<< neventstest << std::endl;
}
