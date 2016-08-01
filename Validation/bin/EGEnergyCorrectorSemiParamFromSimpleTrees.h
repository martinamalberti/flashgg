//--------------------------------------------------------------------------------------------------
// $Id $
//
// SCEnergyCorrectorSemiParm
//
// Helper Class for applying regression-based energy corrections with optimized BDT implementation
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef EGAMMATOOLS_SCEnergyCorrectorSemiParm_H
#define EGAMMATOOLS_SCEnergyCorrectorSemiParm_H
    
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloTopology/interface/CaloTopology.h" 
//#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include <string>

struct photonVariables{

  int numberOfClusters;
  bool missing_clusters;
  float rawEnergy;
  float full5x5_r9;
  float r9;
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

};



class EGEnergyCorrectorSemiParam {
 public:
  EGEnergyCorrectorSemiParam();
  ~EGEnergyCorrectorSemiParam();

  void setGBRForestFromTooFile( const std::string rootfilename, int config = 0 );
  std::pair<double,double> CorrectedEnergyAndResolution(struct photonVariables phoVars, double rho, int nVtx);

 
  
 protected:    
  const GBRForestD *foresteb_;
  const GBRForestD *forestee_;
  const GBRForestD *forestsigmaeb_;
  const GBRForestD *forestsigmaee_;

  int _forestConfig;
};




#include <vdt/vdtMath.h>

#include "TFile.h"
#include <iostream>
using namespace reco;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParam::EGEnergyCorrectorSemiParam() :
  foresteb_(0),
  forestee_(0),
  forestsigmaeb_(0),
  forestsigmaee_(0),
  _forestConfig(0)
{}

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParam::~EGEnergyCorrectorSemiParam()
{
  if( foresteb_ != 0 ) delete foresteb_;
  if( forestee_ != 0 ) delete forestee_;
  if( forestsigmaeb_ != 0 ) delete forestsigmaeb_;
  if( forestsigmaee_ != 0 ) delete forestsigmaee_;


}


void EGEnergyCorrectorSemiParam::setGBRForestFromTooFile( const std::string rootfilename, int config ) {
  edm::FileInPath infile(rootfilename.c_str());
  
  //load forests from file
  TFile *fgbr = TFile::Open(infile.fullPath().c_str(),"READ");    
  fgbr->GetObject("EBCorrection", foresteb_);
  fgbr->GetObject("EECorrection", forestee_);
  //  fgbr->GetObject("EBUncertainty", forestsigmaeb_);
  //  fgbr->GetObject("EEUncertainty", forestsigmaee_);
  fgbr->Close();   

  _forestConfig = config;


}


 



std::pair<double,double> EGEnergyCorrectorSemiParam::CorrectedEnergyAndResolution(struct photonVariables phoVars, 
										   double rho, int nVtx) {

  std::array<float, 39> eval;

  const int numberOfClusters =  phoVars.numberOfClusters;
  const bool missing_clusters =  phoVars.missing_clusters;

  if( missing_clusters ) return  std::pair<double,double>(-999,-999); // do not apply corrections in case of missing info (slimmed MiniAOD electrons)

  const double raw_energy = phoVars.rawEnergy;

  int ivar = 0;

  // SET INPUTS
  eval[ivar++]  = raw_energy;
  if( _forestConfig == 0 ) eval[ivar++]  = phoVars.scEta;
  if( _forestConfig == 0 ) eval[ivar++]  = phoVars.scPhi;
  eval[ivar++]  = phoVars.r9;
  eval[ivar++]  = phoVars.etaWidth;
  eval[ivar++]  = phoVars.phiWidth;
  eval[ivar++]  = std::max(0,numberOfClusters - 1);
  if( _forestConfig != 3 ) eval[ivar++]  = phoVars.hadronicOverEm;
  eval[ivar++] = rho;
  eval[ivar++] = nVtx;
  eval[ivar++] = phoVars.seedEta-phoVars.scEta;
  eval[ivar++] = reco::deltaPhi(phoVars.seedPhi,phoVars.scPhi);
  eval[ivar++] = phoVars.seedEnergy/raw_energy;
  eval[ivar++] = phoVars.e3x3/phoVars.e5x5;
  eval[ivar++] = phoVars.sigmaIetaIeta;
  eval[ivar++] = phoVars.sigmaIphiIphi;
  eval[ivar++] = phoVars.sigmaIetaIphi/(phoVars.sigmaIphiIphi*phoVars.sigmaIetaIeta);
  eval[ivar++] = phoVars.maxEnergyXtal/phoVars.e5x5;
  eval[ivar++] = phoVars.e2nd/phoVars.e5x5;
  eval[ivar++] = phoVars.eTop/phoVars.e5x5;
  eval[ivar++] = phoVars.eBottom/phoVars.e5x5;
  eval[ivar++] = phoVars.eLeft/phoVars.e5x5;
  eval[ivar++] = phoVars.eRight/phoVars.e5x5;
  eval[ivar++] = phoVars.e2x5Max/phoVars.e5x5;
  eval[ivar++] = phoVars.e2x5Left/phoVars.e5x5;
  eval[ivar++] = phoVars.e2x5Right/phoVars.e5x5;
  eval[ivar++] = phoVars.e2x5Top/phoVars.e5x5;
  eval[ivar++] = phoVars.e2x5Bottom/phoVars.e5x5;

  const bool iseb = phoVars.isEB;
  if (iseb) {
    eval[ivar++] = phoVars.e5x5/phoVars.seedEnergy;
    eval[ivar++] = phoVars.ieta;
    eval[ivar++] = phoVars.iphi;

    eval[ivar++] = (phoVars.ieta-1*abs(phoVars.ieta)/phoVars.ieta)%5;
    eval[ivar++] = (phoVars.iphi-1)%2;
    eval[ivar++] = (abs(phoVars.ieta)<=25)*((phoVars.ieta-1*abs(phoVars.ieta)/phoVars.ieta)%25) + (abs(phoVars.ieta)>25)*((phoVars.ieta-26*abs(phoVars.ieta)/phoVars.ieta)%20);
    eval[ivar++] = (phoVars.iphi-1)%20;
    if( _forestConfig == 0 ) eval[ivar++] = phoVars.cryPhi;
    if( _forestConfig == 0 ) eval[ivar++] = phoVars.cryEta;
    eval[ivar++] = phoVars.ieta;
    eval[ivar++] = phoVars.iphi;
    
  } else {
    if( _forestConfig != 2 && _forestConfig != 3 ) {
      eval[ivar++] = phoVars.preshowerEnergy/raw_energy;
      eval[ivar++] = phoVars.preshowerEnergyPlane1/ raw_energy;
      eval[ivar++] = phoVars.preshowerEnergyPlane2/ raw_energy;
    }
    //    EEDetId eeseedid(theseed->seed());
    //    eval[ivar++] = eeseedid.ix();
    //    eval[ivar++] = eeseedid.iy();
    
    if( _forestConfig == 0 ) eval[ivar++] = phoVars.cryEta;
    if( _forestConfig == 0 ) eval[ivar++] = phoVars.cryPhi;
    eval[ivar++] = phoVars.ieta;
    eval[ivar++] = phoVars.iphi;
   
  }
  //  std::string ecalStr[2] = {"ee", "eb"};
  //  std::cout << "  config: " << _forestConfig << " #vars: " << ivar << " Ecal: " << ecalStr[int(iseb)]  << std::endl;

  //magic numbers for MINUIT-like transformation of BDT output onto limited range
  //(These should be stored inside the conditions object in the future as well)
  const double meanlimlow  = 0.2;
  const double meanlimhigh = 2.0;
  const double meanoffset  = meanlimlow + 0.5*(meanlimhigh-meanlimlow);
  const double meanscale   = 0.5*(meanlimhigh-meanlimlow);

  const double sigmalimlow  = 0.0002;
  const double sigmalimhigh = 0.5;
  const double sigmaoffset  = sigmalimlow + 0.5*(sigmalimhigh-sigmalimlow);
  const double sigmascale   = 0.5*(sigmalimhigh-sigmalimlow);

  //these are the actual BDT responses
  const GBRForestD *cor_ph_forest = 0;
  //  const GBRForestD *sig_ph_forest = 0;
  if( iseb ) {
    cor_ph_forest = foresteb_;
    //    sig_ph_forest = forestsigmaeb_;
  } else {
    cor_ph_forest = forestee_;
    //    sig_ph_forest = forestsigmaee_;
  }
  //these are the actual BDT responses
  double rawcor = cor_ph_forest->GetResponse(eval.data());
  double rawsig = -1;//sig_ph_forest->GetResponse(eval.data());
  //apply transformation to limited output range (matching the training)
  double cor = meanoffset  + meanscale *vdt::fast_sin(rawcor);
  double sig = sigmaoffset + sigmascale*vdt::fast_sin(rawsig);

  double      ecor = cor*eval[0];
  if( !iseb && 
      _forestConfig != 2 && 
      _forestConfig != 3 )  ecor = cor*(eval[0]+phoVars.preshowerEnergy);

  double sigmacor = sig*ecor;
  
  return std::pair<double,double>(ecor,sigmacor);

}


#endif
