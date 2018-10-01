#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TRandom2.h"
#include "TLorentzVector.h"
///#include "wzAnalyzer/factors.h"

void nanoAOD() {

  vector<TString>infilenamev;
  vector<Int_t> infilecatv;
  vector<double> infilexsecv;

  //  TString outputDirectory = "./";
  // char output[200]; 
  //sprintf(output,outputDirectory+"outAnalysis.root");
  TFile* outFile = new TFile("output.root","recreate");
  //  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/NanoTestPost6/180429_134358/0000/tree_1.root");
  //  infilecatv.push_back(6); 
  // infilexsecv.push_back(0.0175);

  
  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/DoubleMuon/merged.root");
  infilecatv.push_back(0); 
  infilexsecv.push_back(1);


  for(unsigned int ifile=0; ifile<1./*infilenamev.size()*/; ifile++){
    cout<<infilenamev[ifile]<<endl;
    TFile *the_input_file = TFile::Open((TString)infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("Events");
    TTree *the_input_runs = (TTree*)the_input_file->FindObjectAny("Runs");

    int numberOfEvents = the_input_tree->GetEntriesFast();
    cout<<"number of Events: "<<numberOfEvents<<endl;


    TH1D* histo_dimuon_mass = new TH1D( "histo_dimuon_mass"  , "histo_dimuon_mass"  , 400, 0, 400); 
    histo_dimuon_mass->Sumw2();
    
    UInt_t  nMuon;
    Float_t Muon_pt[50];
    Float_t Muon_eta[50];
    Float_t Muon_phi[50]; 
    Float_t Muon_mass[50];
    Bool_t  Muon_tightId[50];  
    Float_t Muon_pfRelIso04_all[50]; 
    
    the_input_tree->SetBranchAddress("nMuon",&nMuon); 
    the_input_tree->SetBranchAddress("Muon_pt",Muon_pt); 
    the_input_tree->SetBranchAddress("Muon_eta",Muon_eta); 
    the_input_tree->SetBranchAddress("Muon_phi",Muon_phi); 
    the_input_tree->SetBranchAddress("Muon_mass",Muon_mass); 
    the_input_tree->SetBranchAddress("Muon_tightId",Muon_tightId); 
    the_input_tree->SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all); 
    
    Int_t PV_npvs;
    the_input_tree->SetBranchAddress("PV_npvs",&PV_npvs); 
    
    Bool_t hltIsoMu24;
    the_input_tree->SetBranchAddress("HLT_IsoMu24",&hltIsoMu24); 
    Bool_t hltMu17Mu8;
    the_input_tree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",&hltMu17Mu8);
    
    
    for (int i=0; i < numberOfEvents; ++i) {
      
      the_input_tree->GetEntry(i);
      if(i%10000==0) cout<<"Event: "<<i<<endl;
      
      if(!(hltIsoMu24 || hltMu17Mu8)) continue;
      TLorentzVector mu1(0,0,0,0);
      TLorentzVector mu2(0,0,0,0);
      TLorentzVector dimu(0,0,0,0);

      
      for(unsigned int m=0; m<nMuon; m++) {
	for(unsigned int n=m+1; n<nMuon; n++) {
       	  if(Muon_pt[m]>20. && Muon_pt[n]>10.) {
	    //&& TMath::Abs(Muon_eta[m])<2.4 
	    mu1.SetPtEtaPhiM(Muon_pt[m],Muon_eta[m],Muon_phi[m],Muon_mass[m]);
	    mu2.SetPtEtaPhiM(Muon_pt[n],Muon_eta[n],Muon_phi[n],Muon_mass[n]);
	    dimu = (mu1 + mu2) ;
	    }
	}
      }
      
      histo_dimuon_mass->Fill(dimu.M());
      
      
    }

  outFile->cd();
  histo_dimuon_mass->Write();


  }

}
