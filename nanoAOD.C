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

void nanoAOD() {
  
  vector<TString>infilenamev;
  vector<Int_t> infilecatv;
  vector<double> infilexsecv;
  
  TFile* outFile = new TFile("output.root","recreate");
  
  TH1D* massforLimit_CatA[300];
  TH1D* massforLimit_CatB[300];
  TH1D* massforLimit_CatC[300];
  float m=10.;
  for(int j=0; j<300.; j++){
    m = m+(m*0.01); 
    massforLimit_CatA[j] = new TH1D(Form("massforLimit_CatA%d",j),Form("massforLimit_CatA%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatA[j]->Sumw2();
    massforLimit_CatB[j] = new TH1D(Form("massforLimit_CatB%d",j),Form("massforLimit_CatB%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatB[j]->Sumw2();
    massforLimit_CatC[j] = new TH1D(Form("massforLimit_CatC%d",j),Form("massforLimit_CatC%d",j),100,m-(m*0.01*10.),m+(m*0.01*10.));  massforLimit_CatC[j]->Sumw2();
    //cout<<m<<"  "<<j<<endl;
  }
  
  TH1D* histo_dimuon_mass = new TH1D( "histo_dimuon_mass"  , "histo_dimuon_mass"  , 400, 0, 400); 
  histo_dimuon_mass->Sumw2();
  
  TH1D* histo_dimuon_pt = new TH1D( "histo_dimuon_pt"  , "histo_dimuon_pt"  , 400, 0, 400); 
  histo_dimuon_pt->Sumw2();
  
  TH1D* histo_dimuon_eta = new TH1D( "histo_dimuon_eta"  , "histo_dimuon_eta", 600, -3, 3); 
  histo_dimuon_eta->Sumw2();

  //run B  
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017B-31Mar2018-v1/181004_141919/0000/tree_1.root");  
  //  infilecatv.push_back(0);  infilexsecv.push_back(1);
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017B-31Mar2018-v1/181004_141919/0000/tree_2.root");  

  //run C 
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017C-31Mar2018-v1/181004_141948/0000/tree_1.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017C-31Mar2018-v1/181004_141948/0000/tree_2.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017C-31Mar2018-v1/181004_141948/0000/tree_3.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017C-31Mar2018-v1/181004_141948/0000/tree_4.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017C-31Mar2018-v1/181004_141948/0000/tree_5.root");
  
  //run D
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017D-31Mar2018-v1/181004_142018/0000/tree_1.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017D-31Mar2018-v1/181004_142018/0000/tree_2.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017D-31Mar2018-v1/181004_142018/0000/tree_3.root");

  //run E
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_1.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_2.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_3.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_4.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_5.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_6.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_7.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_8.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_9.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_10.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_11.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_12.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_13.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_14.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask_Retry_RunE/DoubleMuon/Run2017E-31Mar2018-v1/181005_082406/0000/tree_15.root");

  //run F
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_1.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_2.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_3.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_4.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_5.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_6.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_7.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_8.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_9.root");
  infilenamev.push_back("/eos/cms/store/group/phys_exotica/darkPhoton/swagata/Oct4_without_lumiMask/DoubleMuon/Run2017F-31Mar2018-v1/181004_142119/0000/tree_10.root");


  int num=infilenamev.size();

  std::cout << "Number of files " << num << std::endl;

  for(unsigned int ifile=0; ifile<num; ifile++){
    cout<<infilenamev[ifile]<<endl;
    TFile *the_input_file = TFile::Open((TString)infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("Events");
    TTree *the_input_runs = (TTree*)the_input_file->FindObjectAny("Runs");

    int numberOfEvents = the_input_tree->GetEntriesFast();
    cout<<"number of Events: "<<numberOfEvents<<endl;

    UInt_t  nMuon;
    Float_t Muon_pt[50];
    Float_t Muon_eta[50];
    Float_t Muon_phi[50]; 
    Float_t Muon_mass[50];
    Bool_t  Muon_tightId[50];  
    Float_t Muon_pfRelIso04_all[50]; 
    Int_t Muon_charge[50];
    
    the_input_tree->SetBranchAddress("nMuon",&nMuon); 
    the_input_tree->SetBranchAddress("Muon_pt",Muon_pt); 
    the_input_tree->SetBranchAddress("Muon_eta",Muon_eta); 
    the_input_tree->SetBranchAddress("Muon_phi",Muon_phi); 
    the_input_tree->SetBranchAddress("Muon_mass",Muon_mass); 
    the_input_tree->SetBranchAddress("Muon_tightId",Muon_tightId); 
    the_input_tree->SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all); 
    the_input_tree->SetBranchAddress("Muon_charge",Muon_charge);
    
    Int_t PV_npvs;
    Int_t PV_npvsGood;
    the_input_tree->SetBranchAddress("PV_npvs",&PV_npvs); 
    the_input_tree->SetBranchAddress("PV_npvsGood",&PV_npvsGood);
    
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

      int pair_found=0;
      float maxEta=99;
	  
      for(unsigned int m1=0; m1<nMuon; m1++) {
	for(unsigned int n=m1+1; n<nMuon; n++) {
	  //	  std::cout << "mu pt m, n " << Muon_pt[m] << "   " << Muon_pt[n] << std::endl;
	  if (pair_found>0) continue;
	  
	  if (PV_npvsGood>0) {
	    if(Muon_pt[m1]>19. && Muon_pt[n]>10.) {
	      if ( TMath::Abs(Muon_eta[m1])<2.4  &&  TMath::Abs(Muon_eta[n])<2.4 ) {  
		if (Muon_tightId[m1] && Muon_tightId[n] ) { 
		  if (Muon_pfRelIso04_all[m1]<0.15 && Muon_pfRelIso04_all[n]<0.15) {
		    if ( Muon_charge[m1]*Muon_charge[n] < 0 ) {
		      mu1.SetPtEtaPhiM(Muon_pt[m1],Muon_eta[m1],Muon_phi[m1],Muon_mass[m1]);
		      mu2.SetPtEtaPhiM(Muon_pt[n],Muon_eta[n],Muon_phi[n],Muon_mass[n]);
		      dimu = (mu1 + mu2) ;
		      maxEta=TMath::Max(abs(Muon_eta[m1]),abs(Muon_eta[n]));
		      pair_found+=1;
		      // std::cout << "muon pair found, mass " << dimu.M() << std::endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      histo_dimuon_mass->Fill(dimu.M());
      
      float slidePt1 = dimu.M()/3.;
      float slidePt2 = dimu.M()/4.;
      
      float ma=10.;
      for(int j=0; j<300.; j++){
	ma = ma+(ma*0.01); 
	if( (dimu.M() >= ma-(ma*0.01*10.)) &&  (dimu.M() < ma+(ma*0.01*10.)) ) {
	  if(mu1.Pt()>slidePt1 && mu2.Pt()>slidePt2 && maxEta<0.9 ){ massforLimit_CatA[j]->Fill(dimu.M()); }
	  if(mu1.Pt()>slidePt1 && mu2.Pt()>slidePt2 && maxEta>=0.9 && maxEta<1.2){ massforLimit_CatB[j]->Fill(dimu.M()); }
	  if(mu1.Pt()>slidePt1 && mu2.Pt()>slidePt2 && maxEta>=1.2 && maxEta<2.4){ massforLimit_CatC[j]->Fill(dimu.M()); }
	  
	}
      }
      
      histo_dimuon_pt->Fill(dimu.Pt());
      histo_dimuon_eta->Fill(dimu.Eta());
    }
    
  }
  
  outFile->cd();
  histo_dimuon_mass->Write();
  histo_dimuon_pt->Write();
  histo_dimuon_eta->Write();
  
  
  for(int j=0; j<300.;j++){
    massforLimit_CatA[j]->Write();
    massforLimit_CatB[j]->Write();
    massforLimit_CatC[j]->Write();
  }
  
  
  
}
