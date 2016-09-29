#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "PhysicsTools/UtilAlgos/interface/BasicAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"

#include "TTree.h"

#include <vector>
#include <string>

using namespace std;
using namespace edm;
using namespace flashgg;


// **********************************************************************

// define the structures used to create tree branches and fill the trees

struct eventInfo {

    float weight;
    float puweight;

    int run;
    int event;
    int lumi;
    
    int nvtx;
    int npu;

    int passHLT; 
    
    float pho1_e;
    float pho1_pt;
    float pho1_eta;
    float pho1_phi;
    float pho1_idmva;
    int pho1_genMatchType;

    float pho2_e;
    float pho2_pt;
    float pho2_eta;
    float pho2_phi;
    float pho2_idmva;
    int pho2_genMatchType;

    float dipho_m;
    float dipho_pt;
    float dipho_mva;
    int dipho_vtxind;
    
    vector<float> jet_e;
    vector<float> jet_pt;
    vector<float> jet_eta;
    vector<float> jet_phi;
    vector<float> jet_pujetid;
    vector<float> jet_bdiscriminant;
    vector<int>   jet_hadronFlavour;
    vector<int>   jet_partonFlavour;
    vector<int>   jet_isMatchedToGen;

    vector<float> ele_e;
    vector<float> ele_pt;
    vector<float> ele_eta;
    vector<float> ele_phi;
    vector<float> ele_idmva;
    vector<float> ele_drGsfToPho1;
    vector<float> ele_drGsfToPho2;
    vector<float> ele_iso;
    vector<float> ele_dz;
    vector<float> ele_d0;
    vector<int> ele_nMissingHits;
    vector<int> ele_passCutBasedIdLoose;
    vector<int> ele_isMatchedToGen;
    vector<int> ele_charge;
    
    vector<float> mu_pt;
    vector<float> mu_eta;
    vector<float> mu_phi;
    vector<float> mu_iso;
    vector<bool> mu_isTight;
    vector<bool> mu_isMedium;
    vector<bool> mu_isLoose;
    vector<int> mu_isMatchedToGen;
    vector<int> mu_charge;

    float met;
    float metx;
    float mety;
    float metphi;
    float metSumEt;

    float uncorr_met;
    float uncorr_metx;
    float uncorr_mety;
    float uncorr_metphi;
    float uncorr_metSumEt;

};
// ******************************************************************************************

// ******************************************************************************************
float passCutBasedElectronIdLoose(edm::Ptr<flashgg::Electron> electron, const std::vector<edm::Ptr<reco::Vertex> > &pvPointers){

  bool pass = false;

  float elfull5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();
  float eldEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  float eldPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
  float elhOverE = electron->hcalOverEcal();
  //float elRelIsoEA = electron->standardHggIso()/electron->pt(); // for 76X microAODs
  float elRelIsoEA = electron->standardHggIso();// for 80X microAODs it is laready stored as relative isolation

  float elooEmooP =-999 ; 

  if( electron->ecalEnergy() == 0 ){
    elooEmooP = 1e30;
  }else if( !std::isfinite(electron->ecalEnergy())){    
    elooEmooP = 1e30;
  }else{
    elooEmooP = fabs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy() );
  }

  double vtx_dz = 1000000.;
  unsigned int min_dz_vtx = -1;
        
  for( unsigned int vtxi = 0; vtxi < pvPointers.size(); vtxi++ ) {            
    Ptr<reco::Vertex> vtx = pvPointers[vtxi];            
    if( vtx_dz > fabs(electron->gsfTrack()->dz( vtx->position() )) ) {                
      vtx_dz = fabs( electron->gsfTrack()->dz( vtx->position() ) );
      min_dz_vtx = vtxi;
    }
  }
        
  Ptr<reco::Vertex> best_vtx_elec = pvPointers[min_dz_vtx];
  float elDxy = fabs( electron->gsfTrack()->dxy( best_vtx_elec->position()) ) ;
  float elDz = fabs( electron->gsfTrack()->dz( best_vtx_elec->position())) ;

  int elMissedHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);

  bool isEB = (fabs(electron->superCluster()->eta())<1.479);

  if (isEB){
    if ( elfull5x5_sigmaIetaIeta < 0.0103
	 && fabs(eldEtaIn) < 0.0105
	 && fabs(eldPhiIn) < 0.115
	 && elhOverE < 0.104
	 && elRelIsoEA < 0.0893
	 && elooEmooP < 0.102
	 && fabs(elDxy) < 0.0261
	 && fabs(elDz) < 0.41
	 && elMissedHits <=2 ) 
      pass =  true;
  }
  else {
    if ( elfull5x5_sigmaIetaIeta < 0.0301
	 && fabs(eldEtaIn) < 0.00814
	 && fabs(eldPhiIn) < 0.182
	 && elhOverE < 0.0897
	 && elRelIsoEA < 0.121
	 && elooEmooP < 0.126
	 && fabs(elDxy) < 0.118
	 && fabs(elDz) < 0.822
	 && elMissedHits <=1 )
      pass = true; 
  }

  return pass;

}



// ******************************************************************************************
float electronIsolation(edm::Ptr<flashgg::Electron> electron, double rho){
    // -- compute combined relative isolation: IsoCh + max( 0.0, IsoNh + IsoPh - PU ) )/pT, PU = rho * Aeff 
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    // effetive areas:  https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
    float Aeff = 0;
    float eta = fabs(electron->eta());
    if( eta <  1.0 )                  { Aeff = 0.1752; }
    if( eta >= 1.0   && eta < 1.479 ) { Aeff = 0.1862; }
    if( eta >= 1.479 && eta < 2.0 )   { Aeff = 0.1411; }
    if( eta >= 2.0   && eta < 2.2 )   { Aeff = 0.1534; }
    if( eta >= 2.2   && eta < 2.3 )   { Aeff = 0.1903; }
    if( eta >= 2.3   && eta < 2.4 )   { Aeff = 0.2243; }
    if( eta >= 2.4 )                  { Aeff = 0.2687; }

    //float iso = electron->chargedHadronIso() + std::max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. );  //???? 
    reco::GsfElectron::PflowIsolationVariables pfIso = electron->pfIsolationVariables();
    float iso = pfIso.sumChargedHadronPt + std::max( pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * Aeff, 0. );

    //cout << electron->chargedHadronIso() << "  " <<  pfIso.sumChargedHadronPt << "   pt = " << electron->pt() << endl; 
    //cout << electron->neutralHadronIso() << "  " << pfIso.sumNeutralHadronEt << endl;
    //cout << electron->photonIso() << "  " << pfIso.sumPhotonEt <<endl;
    //cout << electron->chargedHadronIso() + std::max( electron->neutralHadronIso() + electron->photonIso() - rho * Aeff, 0. ) << "   "<< iso<< endl;

    return (iso/ electron->pt());
    
}
// ******************************************************************************************
int electronMatchingToGen(edm::Ptr<flashgg::Electron> electron,  Handle<View<reco::GenParticle> > genParticles){

    int mcmatch = 0;
    for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
        Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
        if ( fabs(gen->pdgId()) != 11 ) continue;
        if ( !(gen->isPromptFinalState())) continue;
        float dR = deltaR( electron->eta(), electron->phi(), gen->eta(), gen->phi() );
        if (dR < 0.1){ //??? 0.1 ok???
            mcmatch = 1;
        }
    }
    return (mcmatch);
}
// ******************************************************************************************


// ******************************************************************************************
int muonMatchingToGen(edm::Ptr<flashgg::Muon> muon, Handle<View<reco::GenParticle> > genParticles){

    int mcmatch = 0;
    for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
        Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
        if ( fabs(gen->pdgId()) != 13 ) continue;
        if ( !(gen)->isPromptFinalState()) continue;
        float dR = deltaR( muon->eta(), muon->phi(), gen->eta(), gen->phi() );
        //cout << " *** Found muon:  ***"<<endl;
        //cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
        //cout << "dR = " << dR <<endl;
        if (dR < 0.1){ //??? 0.1 ok???
            mcmatch = 1;
        }
    }
    return (mcmatch);
}
// ******************************************************************************************


// ******************************************************************************************
float getPhotonEffectiveArea(float eta){
    float effectiveArea = 0;

    if (fabs(eta) < 1. ) effectiveArea = 0.0725;
    if (fabs(eta) > 1.000 && fabs(eta) < 1.479) effectiveArea = 0.0604;
    if (fabs(eta) > 1.479 && fabs(eta) < 2.000) effectiveArea = 0.0320;
    if (fabs(eta) > 2.000 && fabs(eta) < 2.000) effectiveArea = 0.0512;
    if (fabs(eta) > 2.200 && fabs(eta) < 2.300) effectiveArea = 0.0766;
    if (fabs(eta) > 2.300 && fabs(eta) < 2.400) effectiveArea = 0.0949;
    if (fabs(eta) > 2.400) effectiveArea = 0.1160; 

    return effectiveArea;
}

// ******************************************************************************************


// ******************************************************************************************
bool passDiphotonPreselection(edm::Ptr<flashgg::DiPhotonCandidate> dipho, double rho){
    
    if (dipho->leadingPhoton()->pt() < 30.) return false;
    if (dipho->subLeadingPhoton()->pt() < 20.) return false;
    if (fabs(dipho->leadingPhoton()-> superCluster()->eta()) > 2.5) return false;
    if (fabs(dipho->subLeadingPhoton()-> superCluster()->eta()) > 2.5) return false;
    if (fabs(dipho->leadingPhoton()-> superCluster()->eta()) > 1.4442 && fabs(dipho->leadingPhoton()-> superCluster()->eta()) < 1.566) return false; // EB-EE gap veto
    if (fabs(dipho->subLeadingPhoton()-> superCluster()->eta()) > 1.4442 && fabs(dipho->subLeadingPhoton()-> superCluster()->eta()) < 1.566) return false; // EB-EE gap veto
    if (dipho->mass() < 95.) return false;
    if (dipho->leadingPhoton()->hadronicOverEm() > 0.08 ) return false;
    if (dipho->subLeadingPhoton()->hadronicOverEm() > 0.08 ) return false;
    if (!dipho->leadingPhoton()->passElectronVeto()) return false;
    if (!dipho->subLeadingPhoton()->passElectronVeto()) return false;

    // if at least one photon is low R9 (cat1,cat3) --> additional cuts on the low R9 leg  
    if ( (dipho->leadingPhoton()->isEB() && dipho->leadingPhoton()->full5x5_r9()<0.85) || 
         (dipho->subLeadingPhoton()->isEB() && dipho->subLeadingPhoton()->full5x5_r9()<0.85) ||
         (dipho->leadingPhoton()->isEE() && dipho->leadingPhoton()->full5x5_r9()<0.90) ||
         (dipho->subLeadingPhoton()->isEE() && dipho->subLeadingPhoton()->full5x5_r9()<0.90) )
        {
            if ( dipho->leadingPhoton()->isEB() && dipho->leadingPhoton()->full5x5_r9()<0.85 ){
                if (dipho->leadingPhoton()->full5x5_r9() < 0.5) return false;
                if (dipho->leadingPhoton()->full5x5_sigmaIetaIeta() > 0.015) return false;
                float pfPhoIso = dipho->leadingPhoton()->pfPhoIso03() - rho*getPhotonEffectiveArea(dipho->leadingPhoton()-> superCluster()->eta());
                if (pfPhoIso > 4.  ) return false;
                if (dipho->leadingPhoton()->trkSumPtHollowConeDR03() > 6. ) return false;
            }

            if (dipho->subLeadingPhoton()->isEB() && dipho->subLeadingPhoton()->full5x5_r9()<0.85){
                if (dipho->subLeadingPhoton()->full5x5_r9() < 0.5) return false;
                if (dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta() > 0.015) return false;
                float pfPhoIso = dipho->subLeadingPhoton()->pfPhoIso03() - rho*getPhotonEffectiveArea(dipho->subLeadingPhoton()-> superCluster()->eta());
                if (pfPhoIso > 4.  ) return false;
                if (dipho->subLeadingPhoton()->trkSumPtHollowConeDR03() > 6. ) return false;
            }

            if (dipho->leadingPhoton()->isEE() && dipho->leadingPhoton()->full5x5_r9()<0.90){
                if (dipho->leadingPhoton()->full5x5_r9() < 0.8) return false;
                if (dipho->leadingPhoton()->full5x5_sigmaIetaIeta() > 0.035) return false;
                float pfPhoIso = dipho->leadingPhoton()->pfPhoIso03() - rho*getPhotonEffectiveArea(dipho->leadingPhoton()-> superCluster()->eta());
                if (pfPhoIso > 4.  ) return false;
                if (dipho->leadingPhoton()->trkSumPtHollowConeDR03() > 6. ) return false;
            }

            if (dipho->subLeadingPhoton()->isEE() && dipho->subLeadingPhoton()->full5x5_r9()<0.90){
                if (dipho->subLeadingPhoton()->full5x5_r9() < 0.8) return false;
                if (dipho->subLeadingPhoton()->full5x5_sigmaIetaIeta() > 0.035) return false;
                float pfPhoIso = dipho->subLeadingPhoton()->pfPhoIso03() - rho*getPhotonEffectiveArea(dipho->subLeadingPhoton()-> superCluster()->eta());
                if (pfPhoIso > 4.  ) return false;
                if (dipho->subLeadingPhoton()->trkSumPtHollowConeDR03() > 6. ) return false;
            }
        }
    
    return true;
}
// ******************************************************************************************

class SimpleTreeMaker : public edm::BasicAnalyzer
{
 public:
  //  SimpleTreeMaker( const edm::ParameterSet & iConfig, TFileDirectory& fs);
  SimpleTreeMaker( const edm::ParameterSet & iConfig, TFileDirectory& fs, edm::ConsumesCollector && cc);
  virtual ~SimpleTreeMaker();
  void beginJob();
  void analyze( const edm::EventBase& event );
  void endJob();
 
 private:
  void initEventStructure();

  TTree *eventTree;
  eventInfo evInfo;
  int ngen;
  int npre;
  int nfullpre;

  EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
  EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  EDGetTokenT<edm::View<PileupSummaryInfo> >  PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
  EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
  std::vector<edm::InputTag> inputTagJets_;
  EDGetTokenT<View<reco::GenJet> > genJetToken_;
  EDGetTokenT<View<Electron> > electronToken_;
  EDGetTokenT<View<Muon> > muonToken_;
  EDGetTokenT<View<pat::MET> > METToken_;
  EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<double> rhoToken_;
  std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > tokenJets_;

  typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

  double electronPtThreshold_;
  double muonPtThreshold_;
  double jetPtThreshold_;
  string bTag_;

  double lumiWeight_;

  bool isControlSample_;


  GlobalVariablesDumper *globalVarsDumper_;
};
// ******************************************************************************************                                                                                                                     


