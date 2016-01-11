#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
    vector<float> ele_iso;
    vector<float> ele_dz;
    vector<float> ele_d0;
    vector<int> ele_isMatchedToGen;
    
    vector<float> mu_pt;
    vector<float> mu_eta;
    vector<float> mu_phi;
    vector<float> mu_iso;
    vector<bool> mu_isTight;
    vector<bool> mu_isMedium;
    vector<bool> mu_isLoose;
    vector<int> mu_isMatchedToGen;


    float met;
    float metx;
    float mety;
    float metphi;

};
// ******************************************************************************************


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


// ******************************************************************************************
class tthOptimizationTreeMaker : public edm::EDAnalyzer
{
public:
    explicit tthOptimizationTreeMaker( const edm::ParameterSet & );
    ~tthOptimizationTreeMaker();
private:
    edm::Service<TFileService> fs_;
    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;
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
    
    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;
    
    edm::InputTag rhoFixedGrid_;
    
    double electronPtThreshold_;
    double muonPtThreshold_;
    double jetPtThreshold_;
    string bTag_;

    double lumiWeight_;

    GlobalVariablesDumper *globalVarsDumper_;
};
// ******************************************************************************************


// ******************************************************************************************
// constructors and destructor
//
tthOptimizationTreeMaker::tthOptimizationTreeMaker( const edm::ParameterSet &iConfig ):
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag>( "genParticleTag" ) ) ),
    genInfoToken_(consumes<GenEventInfoProduct>( iConfig.getParameter<InputTag> ( "generatorInfo" ) ) ),
    PileUpToken_(consumes<View<PileupSummaryInfo> >( iConfig.getParameter<InputTag> ( "PileUpTag" ) ) ), 
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    genJetToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
    electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
    muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) ),
    METToken_( consumes<View<pat::MET> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
    triggerBitsToken_( consumes<edm::TriggerResults>( iConfig.getParameter<InputTag>( "triggerBits" ) ) )
{
    jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", 20. );
    bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
    electronPtThreshold_ = iConfig.getUntrackedParameter<double>( "electronPtThreshold", 20. );
    muonPtThreshold_ = iConfig.getUntrackedParameter<double>( "muonPtThreshold", 20. );
    lumiWeight_ = iConfig.getUntrackedParameter<double>( "lumiWeight", 1000. ); //pb
    rhoFixedGrid_  = iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" );
    globalVarsDumper_ = new GlobalVariablesDumper( iConfig.getParameter<edm::ParameterSet>( "globalVariables" ) );
}

tthOptimizationTreeMaker::~tthOptimizationTreeMaker()
{
}
// ******************************************************************************************


// ******************************************************************************************
// analyzer
//
void tthOptimizationTreeMaker::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{
    
    // access edm objects
    Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken( triggerBitsToken_, triggerBits );
    
    Handle<double> rhoHandle;
    iEvent.getByLabel( rhoFixedGrid_, rhoHandle );
    double rho = *( rhoHandle.product() );
    
    Handle<View<reco::Vertex> > vertices;
    iEvent.getByToken( vertexToken_, vertices );

    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );
    
    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
    iEvent.getByToken( mvaResultToken_, mvaResults );
    
    JetCollectionVector Jets( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByLabel( inputTagJets_[j], Jets[j] );
    }

    Handle<View<flashgg::Electron> > electrons;
    iEvent.getByToken( electronToken_, electrons );

    Handle<View<flashgg::Muon> > muons;
    iEvent.getByToken( muonToken_, muons );

    Handle<View<pat::MET> > METs;
    iEvent.getByToken( METToken_, METs );
    if( METs->size() != 1 )
        { std::cout << "WARNING number of MET is not equal to 1" << std::endl; }
    Ptr<pat::MET> theMET = METs->ptrAt( 0 );

    // only if MC
    Handle<GenEventInfoProduct> genInfo;
    Handle<View<reco::GenJet> > genJets;
    Handle<View< PileupSummaryInfo> > PileupInfos;
    Handle<View<reco::GenParticle> > genParticles;
    if ( !iEvent.isRealData() ) {
        iEvent.getByToken( genInfoToken_, genInfo );
        iEvent.getByToken( genJetToken_, genJets );
        iEvent.getByToken( PileUpToken_, PileupInfos );
        iEvent.getByToken( genParticleToken_, genParticles );
    }
    
    //for( unsigned int i = 0 ; i < genParticles->size(); i++ ) {
    //    Ptr<reco::GenParticle> gen = genParticles->ptrAt(i);
    //    cout << " pdgId = "<< gen->pdgId()<< " prompt final state = "<< gen->isPromptFinalState() << "  status = " << gen->status() << "   isPrompt = " << gen->statusFlags().isPrompt() <<endl;
    //}    
    
    // -- initialize tree
    initEventStructure();
    
    // -- check if event passes HLT: "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v1"  
    const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );
    vector<std::string> const &names = triggerNames.triggerNames();  
    for( unsigned index = 0; index < triggerNames.size(); ++index ) {
        if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95") ) {
            //cout << TString::Format((triggerNames.triggerName( index )).c_str()) << " " << triggerBits->accept( index ) << endl;
            evInfo.passHLT =  triggerBits->accept( index );
        }
    }
       
    // -- pre-select best di-photon pair
    //    * pt cut, diphoton preselection, id mva cut on leading and subleading photons 
    //    * if more then one di-photon candidate, take the one with highest sumpt = pt_lead+pt_sublead (DiPhotonCandidates are ordered by decreasing sumpt) 
    //    * di-pho mva cut not applied, needs optimization
    int bestIndex = -1;
    for ( unsigned int idipho = 0; idipho < diphotons->size(); idipho++){
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( idipho );        
        //if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()!=1  || dipho->subLeadingPhoton()->genMatchType()!=1 )) continue;
        //ngen++;
        //if (! passDiphotonPreselection(dipho, rho)) continue; 
        //npre++;
        if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()==1 && dipho->subLeadingPhoton()->genMatchType()==1 )){
            cout << dipho->leadingPhoton()->pt() << "  " << dipho->subLeadingPhoton()-> pt() <<endl;
            ngen++;
            npre++;
        }
        // - pt threshold
        if (dipho->leadingPhoton()->pt() < dipho->mass()/3. ) continue;
        if (dipho->subLeadingPhoton()->pt() < dipho->mass()/4. ) continue;
        // - photon id mva cut (~99% efficient on signal photons after preselection)
        if (dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < -0.9 ) continue;
        if (dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() ) < -0.9 ) continue;
        bestIndex = idipho;
        if (! iEvent.isRealData() && (dipho->leadingPhoton()->genMatchType()==1 && dipho->subLeadingPhoton()->genMatchType()==1 )){
            nfullpre++;
        }
        break;
    }

    
    
    // -- analyze event if there is at least one good di-photon candidate
    
    if ( bestIndex > -1) {
        
        // -- event weight (Lumi x cross section x gen weight)
        float w = 1.;
        if( ! iEvent.isRealData() ) {
            w = lumiWeight_;
            if( genInfo.isValid() ) {
                const auto &weights = genInfo->weights();
                if( ! weights.empty() ) {
                    w *= weights[0];
                }
            }
        }
        evInfo.weight = w;

        // -- pileup weights
        globalVarsDumper_->fill( iEvent );
        evInfo.run = globalVarsDumper_->cache().run;
        evInfo.lumi = globalVarsDumper_->cache().lumi;
        evInfo.event = globalVarsDumper_->cache().event;

        if( globalVarsDumper_->puReWeight() ) {
            evInfo.puweight = globalVarsDumper_->cache().puweight;
        }


        // -- number of pileup events
        float pu = 0.; 
        if( ! iEvent.isRealData() ) {
            for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
                Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
                if( pu_bunchcrossing == 0 ) {
                    pu = PileupInfos->ptrAt( PVI )->getPU_NumInteractions();
                }
            }
        }
        evInfo.npu = pu;
        
        // -- number of reco vertices
        evInfo.nvtx = vertices->size() ;

        // -- photons
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( bestIndex );
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( bestIndex );
        
        evInfo.pho1_e  = dipho->leadingPhoton()->energy();
        evInfo.pho1_pt  = dipho->leadingPhoton()->pt();
        evInfo.pho1_eta = dipho->leadingPhoton()->eta();
        evInfo.pho1_phi = dipho->leadingPhoton()->phi();
        evInfo.pho1_idmva = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        evInfo.pho1_genMatchType = dipho->leadingPhoton()->genMatchType();
        
        evInfo.pho2_e  = dipho->subLeadingPhoton()->energy();
        evInfo.pho2_pt  = dipho->subLeadingPhoton()->pt();
        evInfo.pho2_eta = dipho->subLeadingPhoton()->eta();
        evInfo.pho2_phi = dipho->subLeadingPhoton()->phi();
        evInfo.pho2_idmva = dipho->subLeadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        evInfo.pho2_genMatchType = dipho->subLeadingPhoton()->genMatchType();
                
        evInfo.dipho_pt  = dipho->pt();
        evInfo.dipho_m   = dipho->mass();
        evInfo.dipho_mva = mvares->result ;
        
        // -- jets
        // take the jets corresponding to the diphoton candidate
        unsigned int jetCollectionIndex = diphotons->ptrAt( bestIndex )->jetCollectionIndex();
        
        int njets = 0;
        
        for( UInt_t ijet = 0; ijet < Jets[jetCollectionIndex]->size() ; ijet++ ) {
            
            Ptr<flashgg::Jet> jet  = Jets[jetCollectionIndex]->ptrAt( ijet );
            
            float dRJetPhoLead = deltaR(jet->eta(), jet->phi(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi());
            float dRJetPhoSubLead = deltaR(jet->eta(), jet->phi(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi());

            if( dRJetPhoLead < 0.4 || dRJetPhoSubLead < 0.4 ) { continue; }
            if( !jet->passesJetID( flashgg::Loose ) ) continue;// pass jet id (reject surios detector noise)
            if( !jet->passesPuJetId(diphotons->ptrAt( bestIndex ))){ continue;} // pass PU jet id (always = 1, not implemented yet)
            if( jet->pt() < jetPtThreshold_ ) { continue; }
            
            njets++;
            
            // matching to gen jets
            int isMatchedToGen = 0; 
            if( ! iEvent.isRealData() ) {
                for( unsigned int jg = 0 ; jg < genJets->size() ; jg++ ) {
                    float dr = deltaR(jet->eta(), jet->phi(), genJets->ptrAt( jg )->eta() , genJets->ptrAt( jg )->phi() );
                    if (dr > 0.4) continue;
                    isMatchedToGen = 1;
                }
            }
            
            evInfo.jet_e.push_back(jet->energy());
            evInfo.jet_pt.push_back(jet->pt());
            evInfo.jet_eta.push_back(jet->eta());
            evInfo.jet_phi.push_back(jet->phi());
            evInfo.jet_bdiscriminant.push_back(jet->bDiscriminator( bTag_ ));
            evInfo.jet_partonFlavour.push_back(jet->partonFlavour());
            evInfo.jet_hadronFlavour.push_back(jet->hadronFlavour());
            evInfo.jet_isMatchedToGen.push_back(isMatchedToGen);
        }
        
        // -- leptons (e, mu)

        // -- electrons
        for (UInt_t iele = 0 ; iele < electrons->size(); iele++){
            edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( iele );
            if (fabs(electron->eta()) > 2.4) { continue; }
            if (electron->pt()  < electronPtThreshold_) { continue; }
            if( electron->hasMatchedConversion() ) { continue; } // remove conversions
            // missing hits: from cut-based selection: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            if( electron->isEB() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 2 ) { continue; } 
            if( electron->isEE() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 1 ) { continue; }

            // add dr photon - gsf track
            
            Ptr<reco::Vertex> ele_vtx = chooseElectronVertex( electron,  vertices->ptrs() );
            float d0 = electron->gsfTrack()->dxy( ele_vtx->position() );
            float dz = electron->gsfTrack()->dz( ele_vtx->position() );
            float isol = electronIsolation(electron, rho); 
            int mcMatch = -1;
            if( ! iEvent.isRealData() )
                mcMatch = electronMatchingToGen(electron, genParticles); 

            evInfo.ele_e.push_back(electron->energy());
            evInfo.ele_pt.push_back(electron->pt());
            evInfo.ele_eta.push_back(electron->eta());
            evInfo.ele_phi.push_back(electron->phi());
            evInfo.ele_idmva.push_back(electron->nonTrigMVA());
            evInfo.ele_iso.push_back(isol);
            evInfo.ele_d0.push_back(d0);
            evInfo.ele_dz.push_back(dz);
            evInfo.ele_isMatchedToGen.push_back(mcMatch);
        }       


        for (UInt_t imu = 0 ; imu < muons->size(); imu++){
            edm::Ptr<flashgg::Muon> muon = muons->ptrAt( imu );
            if (fabs(muon->eta()) > 2.4) { continue; }
            if (muon->pt()  < muonPtThreshold_) { continue; }
            // muon ID and isolation: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
            float muPFCombRelIso = ( muon->pfIsolationR04().sumChargedHadronPt + max( 0.,muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5 * muon->pfIsolationR04().sumPUPt ) ) / ( muon->pt() );

            int vtxInd = 0;
            double dzmin = 9999;
            for( size_t ivtx = 0 ; ivtx < vertices->size(); ivtx++ ) {
                Ptr<reco::Vertex> vtx = vertices->ptrAt(ivtx);
                if( !muon->innerTrack() ) { continue; }
                if( fabs( muon->innerTrack()->vz() - vtx->position().z() ) < dzmin ) {
                    dzmin = fabs( muon->innerTrack()->vz() - vtx->position().z() );
                    vtxInd = ivtx;
                }
            }
            Ptr<reco::Vertex> muonVtx = vertices->ptrAt(vtxInd);

            int mcMatch =  -1;
            if( ! iEvent.isRealData() ) mcMatch = muonMatchingToGen(muon, genParticles); 

            evInfo.mu_pt.push_back(muon->pt());
            evInfo.mu_eta.push_back(muon->eta());
            evInfo.mu_phi.push_back(muon->phi());
            evInfo.mu_iso.push_back(muPFCombRelIso);
            evInfo.mu_isTight.push_back(muon::isTightMuon( *muon, *muonVtx ));
            evInfo.mu_isMedium.push_back(muon::isMediumMuon( *muon ));
            evInfo.mu_isLoose.push_back(muon::isLooseMuon( *muon ));
            evInfo.mu_isMatchedToGen.push_back(mcMatch); 
        }


        // -- MET 
        evInfo.met = theMET->pt();
        evInfo.metx = theMET->px();
        evInfo.mety = theMET->py();
        evInfo.metphi = theMET->phi();

        
        // --- fill the tree
        //if ( njets > 0. ) // fill only if min number of jets?
        eventTree->Fill();
    }
}
// ******************************************************************************************


// ******************************************************************************************
void
tthOptimizationTreeMaker::beginJob()
{
  ngen = 0;
  npre = 0;
  nfullpre =0 ;

  // per-event tree
  eventTree = fs_->make<TTree>( "event", "event" );

  eventTree->Branch( "run", &evInfo.run, "run/I" );
  eventTree->Branch( "lumi", &evInfo.lumi, "lumi/I" );
  eventTree->Branch( "event", &evInfo.event, "event/I" );

  eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
  eventTree->Branch( "puweight", &evInfo.puweight,"puweight/F");

  eventTree->Branch( "npu", &evInfo.npu, "npu/I" );
  eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );

  eventTree->Branch( "passHLT", &evInfo.passHLT, "passHLT/I" );
  
  eventTree->Branch( "pho1_e", &evInfo.pho1_e, "pho1_e/F" );
  eventTree->Branch( "pho1_pt", &evInfo.pho1_pt, "pho1_pt/F" );
  eventTree->Branch( "pho1_eta", &evInfo.pho1_eta, "pho1_eta/F" );
  eventTree->Branch( "pho1_phi", &evInfo.pho1_phi, "pho1_phi/F" );
  eventTree->Branch( "pho1_idmva", &evInfo.pho1_idmva, "pho1_idmva/F" );
  eventTree->Branch( "pho1_genMatchType", &evInfo.pho1_genMatchType, "pho1_genMatchType/I" );

  eventTree->Branch( "pho2_e", &evInfo.pho2_e, "pho2_e/F" );
  eventTree->Branch( "pho2_pt", &evInfo.pho2_pt, "pho2_pt/F" );
  eventTree->Branch( "pho2_eta", &evInfo.pho2_eta, "pho2_eta/F" );
  eventTree->Branch( "pho2_phi", &evInfo.pho2_phi, "pho2_phi/F" );
  eventTree->Branch( "pho2_idmva", &evInfo.pho2_idmva, "pho2_idmva/F" );
  eventTree->Branch( "pho2_genMatchType", &evInfo.pho2_genMatchType, "pho2_genMatchType/I" );

  eventTree->Branch( "dipho_pt", &evInfo.dipho_pt, "dipho_pt/F" );
  eventTree->Branch( "dipho_m", &evInfo.dipho_m, "dipho_m/F" );
  eventTree->Branch( "dipho_mva", &evInfo.dipho_mva, "dipho_mva/F" );

  eventTree->Branch( "jet_e", &evInfo.jet_e);
  eventTree->Branch( "jet_pt", &evInfo.jet_pt);
  eventTree->Branch( "jet_eta", &evInfo.jet_eta);
  eventTree->Branch( "jet_phi", &evInfo.jet_phi);
  eventTree->Branch( "jet_pujetid", &evInfo.jet_pujetid);
  eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant);
  eventTree->Branch( "jet_partonFlavour", &evInfo.jet_partonFlavour);
  eventTree->Branch( "jet_hadronFlavour", &evInfo.jet_hadronFlavour);
  eventTree->Branch( "jet_isMatchedToGen", &evInfo.jet_isMatchedToGen);

  eventTree->Branch( "ele_e", &evInfo.ele_e);
  eventTree->Branch( "ele_pt", &evInfo.ele_pt);
  eventTree->Branch( "ele_eta", &evInfo.ele_eta);
  eventTree->Branch( "ele_phi", &evInfo.ele_phi);
  eventTree->Branch( "ele_idmva", &evInfo.ele_idmva);
  eventTree->Branch( "ele_iso", &evInfo.ele_iso);
  eventTree->Branch( "ele_dz", &evInfo.ele_dz);
  eventTree->Branch( "ele_d0", &evInfo.ele_d0);
  eventTree->Branch( "ele_isMatchedToGen", &evInfo.ele_isMatchedToGen);

  eventTree->Branch( "mu_pt", &evInfo.mu_pt);
  eventTree->Branch( "mu_eta", &evInfo.mu_eta);
  eventTree->Branch( "mu_phi", &evInfo.mu_phi);
  eventTree->Branch( "mu_iso", &evInfo.mu_iso);
  eventTree->Branch( "mu_isTight", &evInfo.mu_isTight);
  eventTree->Branch( "mu_isMedium", &evInfo.mu_isMedium);
  eventTree->Branch( "mu_isLoose", &evInfo.mu_isLoose);
  eventTree->Branch( "mu_isMatchedToGen", &evInfo.mu_isMatchedToGen);

  eventTree->Branch( "met", &evInfo.met);
  eventTree->Branch( "metx", &evInfo.metx);
  eventTree->Branch( "mety", &evInfo.mety);
  eventTree->Branch( "metphi", &evInfo.metphi);

}
// ******************************************************************************************


// ******************************************************************************************
void
tthOptimizationTreeMaker::endJob()
{
    cout << "Total nuber of events before preselection = "<< ngen << endl;
    cout << "Number of events after preselection       = "<< npre << endl;
    cout << "Number of events after full preselection  = "<< nfullpre << endl;
 
} // end of endJob
// ******************************************************************************************


// ******************************************************************************************
void
tthOptimizationTreeMaker::initEventStructure()
{
    // per-event tree:
    evInfo.weight = -999.;
    evInfo.puweight = -999.;

    evInfo.run = -999;
    evInfo.lumi = -999.;
    evInfo.event = -999.;

    evInfo.npu = -999;
    evInfo.nvtx = -999;
    evInfo.passHLT = -1;
    
    evInfo.pho1_e  = -999.;
    evInfo.pho1_pt  = -999.;
    evInfo.pho1_eta = -999.;
    evInfo.pho1_phi = -999.;
    evInfo.pho1_idmva = -999.;
    evInfo.pho1_genMatchType = -999.;
    
    evInfo.pho2_e  = -999.;
    evInfo.pho2_pt  = -999.;
    evInfo.pho2_eta = -999.;
    evInfo.pho2_phi = -999.;
    evInfo.pho2_idmva = -999.;
    evInfo.pho2_genMatchType = -999.;

    evInfo.dipho_pt   = -999.;
    evInfo.dipho_m    = -999.;
    evInfo.dipho_mva  = -999.;
    
    evInfo.jet_e .clear();
    evInfo.jet_pt .clear();
    evInfo.jet_eta .clear();
    evInfo.jet_phi .clear();
    evInfo.jet_pujetid .clear();
    evInfo.jet_bdiscriminant .clear();
    evInfo.jet_partonFlavour .clear();
    evInfo.jet_hadronFlavour .clear();
    evInfo.jet_isMatchedToGen .clear();
    
    evInfo.ele_e .clear();
    evInfo.ele_pt .clear();
    evInfo.ele_eta .clear();
    evInfo.ele_phi .clear();
    evInfo.ele_idmva .clear();
    evInfo.ele_iso .clear();
    evInfo.ele_dz .clear();
    evInfo.ele_d0 .clear();
    evInfo.ele_isMatchedToGen .clear();

    evInfo.mu_pt .clear();
    evInfo.mu_eta .clear();
    evInfo.mu_phi .clear();
    evInfo.mu_iso .clear();
    evInfo.mu_isTight .clear();
    evInfo.mu_isMedium .clear();
    evInfo.mu_isLoose .clear();
    evInfo.mu_isMatchedToGen .clear();

    evInfo.met = -999;
    evInfo.metx = -999;
    evInfo.mety = -999;
    evInfo.metphi = -999;

}
// ******************************************************************************************



typedef tthOptimizationTreeMaker FlashggtthOptimizationTreeMaker;
DEFINE_FWK_MODULE( FlashggtthOptimizationTreeMaker );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

