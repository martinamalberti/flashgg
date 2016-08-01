#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "./EGEnergyCorrectorSemiParamFromSimpleTrees.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>


using namespace std;

/*struct photonVariables{

  int numberOfClusters;
  float rawEnergy;
  float full5x5_r9;
  float etaWidth;
  float phiWidth;
  float hadronicOverEm;
  float scEta;
  float scPhi;
  float seedEta;
  float seedPhi;
  float seedEnergy;
  float e3x3;
  float e5x5;
  float sigmaIetaIeta;
  float sigmaIphiIphi;
  float sigmaIetaIphi;
  float maxEnergyXtal;
  float e2nd;
  float eTop;
  float eBottom;
  float eLeft;
  float eRight;
  float e2x5Max;
  float e2x5Left;
  float e2x5Right;
  float e2x5Top;
  float e2x5Bottom;
  float preshowerEnergy;
  float preshowerEnergyPlane1;
  float preshowerEnergyPlane2;

  bool isEB;
  
  float cryPhi;
  float cryEta;
  int ieta;
  int iphi;
  
  }*/





int main(int argc, char** argv)
{
  
  TFile *file = TFile::Open(argv[1]);
  cout << "Analyzing " << argv[1] <<endl;

  TTree* tree = (TTree*)file->Get("tagsDumper/trees/test_13TeV_UntaggedTag");

  // Declaration of leaf types
  // Declaration of leaf types                                                                                                                                                      
  Int_t           candidate_id;
  Float_t         weight;
  Float_t         mass;
  Float_t         pt;
  Float_t         diphotonMVA;
  Float_t         pho1_et;
  Float_t         pho1_energy;
  Float_t         pho1_eTrue;
  Float_t         pho1_eta;
  Float_t         pho1_phi;
  Float_t         pho1_idmva;
  Float_t         pho2_et;
  Float_t         pho2_energy;
  Float_t         pho2_eTrue;
  Float_t         pho2_eta;
  Float_t         pho2_phi;
  Float_t         pho2_idmva;

  Float_t         pho1_isEB;
  Float_t         pho1_numberOfClusters;
  //Bool_t          pho1_missing_clusters;
  Float_t         pho1_rawEnergy;
  Float_t         pho1_scEnergy;
  Float_t         pho1_full5x5_r9;
  Float_t         pho1_r9;
  Float_t         pho1_etaWidth;
  Float_t         pho1_phiWidth;
  Float_t         pho1_hadronicOverEm;
  Float_t         pho1_scEta;
  Float_t         pho1_scPhi;
  Float_t         pho1_seedEta;
  Float_t         pho1_seedPhi;
  Float_t         pho1_seedEnergy;
  Float_t         pho1_e3x3;
  Float_t         pho1_e5x5;
  Float_t         pho1_sigmaIetaIeta;
  Float_t         pho1_sigmaIphiIphi;
  Float_t         pho1_sigmaIetaIphi;
  Float_t         pho1_maxEnergyXtal;
  Float_t         pho1_e2nd;
  Float_t         pho1_eTop;
  Float_t         pho1_eBottom;
  Float_t         pho1_eLeft;
  Float_t         pho1_eRight;
  Float_t         pho1_e2x5Max;
  Float_t         pho1_e2x5Left;
  Float_t         pho1_e2x5Right;
  Float_t         pho1_e2x5Top;
  Float_t         pho1_e2x5Bottom;
  Float_t         pho1_preshowerEnergy;
  Float_t         pho1_preshowerEnergyPlane1;
  Float_t         pho1_preshowerEnergyPlane2;
  Float_t         pho1_cryEta;
  Float_t         pho1_cryPhi;
  Float_t         pho1_ieta;
  Float_t         pho1_iphi;
 
  Float_t          pho2_isEB;
  Float_t         pho2_numberOfClusters;
  //  Bool_t          pho2_missing_clusters;
  Float_t         pho2_rawEnergy;
  Float_t         pho2_scEnergy;
  Float_t         pho2_full5x5_r9;
  Float_t         pho2_r9;
  Float_t         pho2_etaWidth;
  Float_t         pho2_phiWidth;
  Float_t         pho2_hadronicOverEm;
  Float_t         pho2_scEta;
  Float_t         pho2_scPhi;
  Float_t         pho2_seedEta;
  Float_t         pho2_seedPhi;
  Float_t         pho2_seedEnergy;
  Float_t         pho2_e3x3;
  Float_t         pho2_e5x5;
  Float_t         pho2_sigmaIetaIeta;
  Float_t         pho2_sigmaIphiIphi;
  Float_t         pho2_sigmaIetaIphi;
  Float_t         pho2_maxEnergyXtal;
  Float_t         pho2_e2nd;
  Float_t         pho2_eTop;
  Float_t         pho2_eBottom;
  Float_t         pho2_eLeft;
  Float_t         pho2_eRight;
  Float_t         pho2_e2x5Max;
  Float_t         pho2_e2x5Left;
  Float_t         pho2_e2x5Right;
  Float_t         pho2_e2x5Top;
  Float_t         pho2_e2x5Bottom;
  Float_t         pho2_preshowerEnergy;
  Float_t         pho2_preshowerEnergyPlane1;
  Float_t         pho2_preshowerEnergyPlane2;
  Float_t         pho2_cryEta;
  Float_t         pho2_cryPhi;
  Float_t         pho2_ieta;
  Float_t         pho2_iphi;
  Float_t         rho;
  Int_t           nvtx;


  tree->SetBranchAddress("candidate_id", &candidate_id);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("diphotonMVA", &diphotonMVA);
  tree->SetBranchAddress("pho1_et", &pho1_et);
  tree->SetBranchAddress("pho1_energy", &pho1_energy);
  tree->SetBranchAddress("pho1_eTrue", &pho1_eTrue);
  tree->SetBranchAddress("pho1_eta", &pho1_eta);
  tree->SetBranchAddress("pho1_phi", &pho1_phi);
  tree->SetBranchAddress("pho1_idmva", &pho1_idmva);
  tree->SetBranchAddress("pho2_et", &pho2_et);
  tree->SetBranchAddress("pho2_energy", &pho2_energy);
  tree->SetBranchAddress("pho2_eTrue", &pho2_eTrue);
  tree->SetBranchAddress("pho2_eta", &pho2_eta);
  tree->SetBranchAddress("pho2_phi", &pho2_phi);
  tree->SetBranchAddress("pho2_idmva", &pho2_idmva);

  tree->SetBranchAddress("pho1_isEB", &pho1_isEB);
  tree->SetBranchAddress("pho1_numberOfClusters", &pho1_numberOfClusters);
  //tree->SetBranchAddress("pho1_missing_clusters", &pho1_missing_clusters);
  tree->SetBranchAddress("pho1_rawEnergy", &pho1_rawEnergy);
  tree->SetBranchAddress("pho1_scEnergy", &pho1_scEnergy);
  tree->SetBranchAddress("pho1_full5x5_r9", &pho1_full5x5_r9);
  tree->SetBranchAddress("pho1_r9", &pho1_r9);
  tree->SetBranchAddress("pho1_etaWidth", &pho1_etaWidth);
  tree->SetBranchAddress("pho1_phiWidth", &pho1_phiWidth);
  tree->SetBranchAddress("pho1_hadronicOverEm", &pho1_hadronicOverEm);
  tree->SetBranchAddress("pho1_scEta", &pho1_scEta);
  tree->SetBranchAddress("pho1_scPhi", &pho1_scPhi);
  tree->SetBranchAddress("pho1_seedEta", &pho1_seedEta);
  tree->SetBranchAddress("pho1_seedPhi", &pho1_seedPhi);
  tree->SetBranchAddress("pho1_seedEnergy", &pho1_seedEnergy);
  tree->SetBranchAddress("pho1_e3x3", &pho1_e3x3);
  tree->SetBranchAddress("pho1_e5x5", &pho1_e5x5);
  tree->SetBranchAddress("pho1_sigmaIetaIeta", &pho1_sigmaIetaIeta);
  tree->SetBranchAddress("pho1_sigmaIphiIphi", &pho1_sigmaIphiIphi);
  tree->SetBranchAddress("pho1_sigmaIetaIphi", &pho1_sigmaIetaIphi);
  tree->SetBranchAddress("pho1_maxEnergyXtal", &pho1_maxEnergyXtal);
  tree->SetBranchAddress("pho1_e2nd", &pho1_e2nd);
  tree->SetBranchAddress("pho1_eTop", &pho1_eTop);
  tree->SetBranchAddress("pho1_eBottom", &pho1_eBottom);
  tree->SetBranchAddress("pho1_eLeft", &pho1_eLeft);
  tree->SetBranchAddress("pho1_eRight", &pho1_eRight);
  tree->SetBranchAddress("pho1_e2x5Max", &pho1_e2x5Max);
  tree->SetBranchAddress("pho1_e2x5Left", &pho1_e2x5Left);
  tree->SetBranchAddress("pho1_e2x5Right", &pho1_e2x5Right);
  tree->SetBranchAddress("pho1_e2x5Top", &pho1_e2x5Top);
  tree->SetBranchAddress("pho1_e2x5Bottom", &pho1_e2x5Bottom);
  tree->SetBranchAddress("pho1_preshowerEnergy", &pho1_preshowerEnergy);
  tree->SetBranchAddress("pho1_preshowerEnergyPlane1", &pho1_preshowerEnergyPlane1);
  tree->SetBranchAddress("pho1_preshowerEnergyPlane2", &pho1_preshowerEnergyPlane2);
  tree->SetBranchAddress("pho1_cryEta", &pho1_cryEta);
  tree->SetBranchAddress("pho1_cryPhi", &pho1_cryPhi);
  tree->SetBranchAddress("pho1_ieta", &pho1_ieta);
  tree->SetBranchAddress("pho1_iphi", &pho1_iphi);

  tree->SetBranchAddress("pho2_isEB", &pho2_isEB);
  tree->SetBranchAddress("pho2_numberOfClusters", &pho2_numberOfClusters);
  //  tree->SetBranchAddress("pho2_missing_clusters", &pho2_missing_clusters);
  tree->SetBranchAddress("pho2_rawEnergy", &pho2_rawEnergy);
  tree->SetBranchAddress("pho2_scEnergy", &pho2_scEnergy);
  tree->SetBranchAddress("pho2_full5x5_r9", &pho2_full5x5_r9);
  tree->SetBranchAddress("pho2_r9", &pho2_r9);
  tree->SetBranchAddress("pho2_etaWidth", &pho2_etaWidth);
  tree->SetBranchAddress("pho2_phiWidth", &pho2_phiWidth);
  tree->SetBranchAddress("pho2_hadronicOverEm", &pho2_hadronicOverEm);
  tree->SetBranchAddress("pho2_scEta", &pho2_scEta);
  tree->SetBranchAddress("pho2_scPhi", &pho2_scPhi);
  tree->SetBranchAddress("pho2_seedEta", &pho2_seedEta);
  tree->SetBranchAddress("pho2_seedPhi", &pho2_seedPhi);
  tree->SetBranchAddress("pho2_seedEnergy", &pho2_seedEnergy);
  tree->SetBranchAddress("pho2_e3x3", &pho2_e3x3);
  tree->SetBranchAddress("pho2_e5x5", &pho2_e5x5);
  tree->SetBranchAddress("pho2_sigmaIetaIeta", &pho2_sigmaIetaIeta);
  tree->SetBranchAddress("pho2_sigmaIphiIphi", &pho2_sigmaIphiIphi);
  tree->SetBranchAddress("pho2_sigmaIetaIphi", &pho2_sigmaIetaIphi);
  tree->SetBranchAddress("pho2_maxEnergyXtal", &pho2_maxEnergyXtal);
  tree->SetBranchAddress("pho2_e2nd", &pho2_e2nd);
  tree->SetBranchAddress("pho2_eTop", &pho2_eTop);
  tree->SetBranchAddress("pho2_eBottom", &pho2_eBottom);
  tree->SetBranchAddress("pho2_eLeft", &pho2_eLeft);
  tree->SetBranchAddress("pho2_eRight", &pho2_eRight);
  tree->SetBranchAddress("pho2_e2x5Max", &pho2_e2x5Max);
  tree->SetBranchAddress("pho2_e2x5Left", &pho2_e2x5Left);
  tree->SetBranchAddress("pho2_e2x5Right", &pho2_e2x5Right);
  tree->SetBranchAddress("pho2_e2x5Top", &pho2_e2x5Top);
  tree->SetBranchAddress("pho2_e2x5Bottom", &pho2_e2x5Bottom);
  tree->SetBranchAddress("pho2_preshowerEnergy", &pho2_preshowerEnergy);
  tree->SetBranchAddress("pho2_preshowerEnergyPlane1", &pho2_preshowerEnergyPlane1);
  tree->SetBranchAddress("pho2_preshowerEnergyPlane2", &pho2_preshowerEnergyPlane2);
  tree->SetBranchAddress("pho2_cryEta", &pho2_cryEta);
  tree->SetBranchAddress("pho2_cryPhi", &pho2_cryPhi);
  tree->SetBranchAddress("pho2_ieta", &pho2_ieta);
  tree->SetBranchAddress("pho2_iphi", &pho2_iphi);

  tree->SetBranchAddress("rho", &rho);
  tree->SetBranchAddress("nvtx", &nvtx);


  // book histograms
  TH1F *h_EBEB_highR9 = new TH1F("h_EBEB_highR9","h_EBEB_highR9",160,70,110);
  TH1F *h_EBEB_lowR9  = new TH1F("h_EBEB_lowR9","h_EBEB_lowR9",160,70,110);

  TH1F *h_EEEE_highR9 = new TH1F("h_EEEE_highR9","h_EEEE_highR9",160,70,110);
  TH1F *h_EEEE_lowR9  = new TH1F("h_EEEE_lowR9","h_EEEE_lowR9",160,70,110);


  TH1F *hcorr_EBEB_highR9 = new TH1F("hcorr_EBEB_highR9","hcorr_EBEB_highR9",160,70,110);
  TH1F *hcorr_EBEB_lowR9  = new TH1F("hcorr_EBEB_lowR9","hcorr_EBEB_lowR9",160,70,110);

  TH1F *hcorr_EEEE_highR9 = new TH1F("hcorr_EEEE_highR9","hcorr_EEEE_highR9",160,70,110);
  TH1F *hcorr_EEEE_lowR9  = new TH1F("hcorr_EEEE_lowR9","hcorr_EEEE_lowR9",160,70,110);
  
  cout << "Number of events to be analyzed : " << tree->GetEntries() <<endl;
  cout <<endl;

  float minEt1   = 40.;
  float minEt2   = 30.;
  float minMass = 70.;
  float r9cut = 0.94;

  photonVariables pho1Vars, pho2Vars;
  EGEnergyCorrectorSemiParam _cor;
  _cor.setGBRForestFromTooFile("cms_wereg_ph_bx25_db.root"   , 1);
  
  for (int i=0; i<tree->GetEntries(); i++) {

    tree -> GetEntry(i);
    if (i%1000000==0) cout << "Analyzing event "<< i <<endl;

    // take only the first diphoton pair (highest SumpT) if there is more than one diphoton candidate
    if (candidate_id !=0 ) continue; 

    // selections                  
    if ( pho1_et < minEt1  ) continue;                                                         
    if ( pho2_et < minEt2  ) continue;                                                         
    if ( mass    < minMass ) continue;

    // weights
    float w = weight;

    
    //set photon variables
    pho1Vars.numberOfClusters = pho1_numberOfClusters;
    //pho1Vars.missing_clusters = pho1_missing_clusters;
    pho1Vars.missing_clusters = 0;
    pho1Vars.rawEnergy = pho1_rawEnergy;
    pho1Vars.full5x5_r9 = pho1_full5x5_r9;
    pho1Vars.r9 = pho1_r9;
    pho1Vars.etaWidth =  pho1_etaWidth;
    pho1Vars.phiWidth =  pho1_phiWidth;
    pho1Vars.hadronicOverEm =  pho1_hadronicOverEm;
    pho1Vars.scEta =  pho1_scEta;
    pho1Vars.scPhi =  pho1_scPhi;
    pho1Vars.seedEnergy =  pho1_seedEnergy;
    pho1Vars.e3x3 =  pho1_e3x3;
    pho1Vars.e5x5 =  pho1_e5x5;
    pho1Vars.sigmaIetaIeta =  pho1_sigmaIetaIeta;
    pho1Vars.sigmaIphiIphi =  pho1_sigmaIphiIphi;
    pho1Vars.sigmaIetaIphi =  pho1_sigmaIetaIphi;
    pho1Vars.maxEnergyXtal = pho1_maxEnergyXtal;
    pho1Vars.e2nd = pho1_e2nd ;
    pho1Vars.eTop = pho1_eTop ;
    pho1Vars.eBottom = pho1_eBottom ;
    pho1Vars.eLeft = pho1_eLeft ;
    pho1Vars.eRight = pho1_eRight ;
    pho1Vars.e2x5Max = pho1_e2x5Max ;
    pho1Vars.e2x5Left = pho1_e2x5Left ;
    pho1Vars.e2x5Right = pho1_e2x5Right ;
    pho1Vars.e2x5Top = pho1_e2x5Top ;
    pho1Vars.e2x5Bottom = pho1_e2x5Bottom ;
    pho1Vars.isEB = bool(pho1_isEB) ;
    pho1Vars.cryPhi = pho1_cryPhi ;
    pho1Vars.cryEta = pho1_cryEta ;
    pho1Vars.ieta = int(pho1_ieta);
    pho1Vars.iphi = int(pho1_iphi);
    pho1Vars.preshowerEnergy = pho1_preshowerEnergy ;
    pho1Vars.preshowerEnergyPlane1 = pho1_preshowerEnergyPlane1 ;
    pho1Vars.preshowerEnergyPlane2 = pho1_preshowerEnergyPlane2 ;

    pho2Vars.numberOfClusters = pho2_numberOfClusters;
    //pho2Vars.missing_clusters = pho2_missing_clusters;
    pho2Vars.missing_clusters = 0;
    pho2Vars.rawEnergy = pho2_rawEnergy;
    pho2Vars.full5x5_r9 = pho2_full5x5_r9;
    pho2Vars.r9 = pho2_r9;
    pho2Vars.etaWidth =  pho2_etaWidth;
    pho2Vars.phiWidth =  pho2_phiWidth;
    pho2Vars.hadronicOverEm =  pho2_hadronicOverEm;
    pho2Vars.scEta =  pho2_scEta;
    pho2Vars.scPhi =  pho2_scPhi;
    pho2Vars.seedEnergy =  pho2_seedEnergy;
    pho2Vars.e3x3 =  pho2_e3x3;
    pho2Vars.e5x5 =  pho2_e5x5;
    pho2Vars.sigmaIetaIeta =  pho2_sigmaIetaIeta;
    pho2Vars.sigmaIphiIphi =  pho2_sigmaIphiIphi;
    pho2Vars.sigmaIetaIphi =  pho2_sigmaIetaIphi;
    pho2Vars.maxEnergyXtal = pho2_maxEnergyXtal;
    pho2Vars.e2nd = pho2_e2nd ;
    pho2Vars.eTop = pho2_eTop ;
    pho2Vars.eBottom = pho2_eBottom ;
    pho2Vars.eLeft = pho2_eLeft ;
    pho2Vars.eRight = pho2_eRight ;
    pho2Vars.e2x5Max = pho2_e2x5Max ;
    pho2Vars.e2x5Left = pho2_e2x5Left ;
    pho2Vars.e2x5Right = pho2_e2x5Right ;
    pho2Vars.e2x5Top = pho2_e2x5Top ;
    pho2Vars.e2x5Bottom = pho2_e2x5Bottom ;
    pho2Vars.isEB = bool(pho2_isEB) ;
    pho2Vars.cryPhi = pho2_cryPhi ;
    pho2Vars.cryEta = pho2_cryEta ;
    pho2Vars.ieta = int(pho2_ieta);
    pho2Vars.iphi = int(pho2_iphi);
    pho2Vars.preshowerEnergy = pho2_preshowerEnergy ;
    pho2Vars.preshowerEnergyPlane1 = pho2_preshowerEnergyPlane1 ;
    pho2Vars.preshowerEnergyPlane2 = pho2_preshowerEnergyPlane2 ;
    
    // recompute energies with regression

    float pho1Ecorr = _cor.CorrectedEnergyAndResolution( pho1Vars , rho, nvtx).first;     
    float pho2Ecorr = _cor.CorrectedEnergyAndResolution( pho2Vars , rho, nvtx).first;     
    float massCorr = mass * sqrt(pho1Ecorr*pho2Ecorr)/sqrt(pho1_energy*pho2_energy);
    

    cout << "pho scEnergy from microAOD =" << pho1_scEnergy << "      after regression on the fly = " << pho1Ecorr<< "   ratio = "<< pho1Ecorr/pho1_scEnergy<< endl; 
    //cout << "pho energy   from microAODs = " << pho2_energy << "      after regression on the fly = " << pho2Ecorr<<endl; 


    // both electrons in the same cat.
    if (fabs(pho1_eta)<1.5 && fabs(pho2_eta)<1.5 ) {
      if (pho1_full5x5_r9 > r9cut && pho2_full5x5_r9 > r9cut) {
	h_EBEB_highR9->Fill(mass,w);
	hcorr_EBEB_highR9->Fill(massCorr,w);
      }
      if (pho1_full5x5_r9 < r9cut && pho2_full5x5_r9 < r9cut) {
	h_EBEB_lowR9->Fill(mass,w);
	hcorr_EBEB_lowR9->Fill(massCorr,w);
      }
    }
    if (fabs(pho1_eta)>1.5 && fabs(pho2_eta)>1.5 ) {
      if (pho1_full5x5_r9 > r9cut && pho2_full5x5_r9 > r9cut) {
	h_EEEE_highR9->Fill(mass,w);
	hcorr_EEEE_highR9->Fill(massCorr,w);
      }
      if (pho1_full5x5_r9 < r9cut && pho2_full5x5_r9 < r9cut) {
	h_EEEE_lowR9->Fill(mass,w);
	hcorr_EEEE_lowR9->Fill(massCorr,w);
      }
    }

  }// end loop over events
  
  
  TFile *fout = new TFile(argv[2],"recreate");


  h_EBEB_highR9->Write();
  h_EBEB_lowR9->Write();
  h_EEEE_highR9->Write();
  h_EEEE_lowR9->Write();

  hcorr_EBEB_highR9->Write();
  hcorr_EBEB_lowR9->Write();
  hcorr_EEEE_highR9->Write();
  hcorr_EEEE_lowR9->Write();

  cout << "Closing file..." << endl;
  fout->Close();

  
}
