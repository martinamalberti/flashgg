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
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

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
    
    int npu;
    int nvtx;
    
    float pho1_pt;
    float pho1_eta;
    float pho1_phi;
    float pho1_idmva;

    float pho2_pt;
    float pho2_eta;
    float pho2_phi;
    float pho2_idmva;

    float dipho_m;
    float dipho_pt;
    float dipho_mva;
    
    vector<float> jet_pt;
    vector<float> jet_eta;
    vector<float> jet_phi;
    vector<float> jet_pujetid;
    vector<float> jet_bdiscriminant;
    vector<int>   jet_isMatchedToGen;

    vector<float> ele_pt;
    vector<float> ele_eta;
    vector<float> ele_phi;
    vector<float> ele_idmva;
    vector<float> ele_iso;
    vector<float> ele_dz;
    vector<float> ele_d0;
    
    vector<float> mu_pt;
    vector<float> mu_eta;
    vector<float> mu_phi;
    vector<float> mu_iso;
    vector<bool> mu_isTight;
    vector<bool> mu_isMedium;
    vector<bool> mu_isLoose;



};
// **********************************************************************


// **********************************************************************
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
    
    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
    std::vector<edm::InputTag> inputTagJets_;
    EDGetTokenT<View<reco::GenJet> > genJetToken_;
    EDGetTokenT<View<Electron> > electronToken_;
    EDGetTokenT<View<Muon> > muonToken_;

    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

    double electronPtThreshold_;
    double muonPtThreshold_;
    double jetPtThreshold_;
    string bTag_;

};
// ******************************************************************************************


// ******************************************************************************************
// constructors and destructor
//
tthOptimizationTreeMaker::tthOptimizationTreeMaker( const edm::ParameterSet &iConfig ):
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( iConfig.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
    genJetToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "GenJetTag" ) ) ),
    electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
    muonToken_( consumes<View<flashgg::Muon> >( iConfig.getParameter<InputTag>( "MuonTag" ) ) )
{
    jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", 20. );
    bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
    electronPtThreshold_ = iConfig.getUntrackedParameter<double>( "electronPtThreshold", 20. );
    muonPtThreshold_ = iConfig.getUntrackedParameter<double>( "muonPtThreshold", 20. );
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
    
    Handle<View<reco::GenJet> > genJets;
    iEvent.getByToken( genJetToken_, genJets );

    Handle<View<flashgg::Electron> > electrons;
    iEvent.getByToken( electronToken_, electrons );

    Handle<View<flashgg::Muon> > muons;
    iEvent.getByToken( muonToken_, muons );

    // initialize tree
    initEventStructure();
    
    // fill variables for photons, jets, etcc...
    
    
    if (diphotons->size() > 0) {
        // event weight
        // boh???
        
        
        // number of vertices
        evInfo.nvtx = vertices->size() ;

        // -- photons
        // take the diphoton candidate with index = 0 (--> highest pt1+pt2)... to be checked..
        int candIndex = 0;
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( candIndex );
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( candIndex );
        
        evInfo.pho1_pt  = dipho->leadingPhoton()->pt();
        evInfo.pho1_eta = dipho->leadingPhoton()->eta();
        evInfo.pho1_phi = dipho->leadingPhoton()->phi();
        evInfo.pho1_idmva = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        
        evInfo.pho2_pt  = dipho->subLeadingPhoton()->pt();
        evInfo.pho2_eta = dipho->subLeadingPhoton()->eta();
        evInfo.pho2_phi = dipho->subLeadingPhoton()->phi();
        evInfo.pho2_idmva = dipho->leadingPhoton()->phoIdMvaDWrtVtx( dipho->vtx() );
        
        evInfo.dipho_pt  = dipho->pt();
        evInfo.dipho_m   = dipho->mass();
        evInfo.dipho_mva = mvares->result ;
        
        // -- jets
        // take the jets corresponding to the diphoton candidate
        unsigned int jetCollectionIndex = diphotons->ptrAt( candIndex )->jetCollectionIndex();
        
        int njets = 0;
        
        for( UInt_t ijet = 0; ijet < Jets[jetCollectionIndex]->size() ; ijet++ ) {
            
            Ptr<flashgg::Jet> jet  = Jets[jetCollectionIndex]->ptrAt( ijet );
            
            float dRJetPhoLead = deltaR(jet->eta(), jet->phi(), dipho->leadingPhoton()->eta(), dipho->leadingPhoton()->phi());
            float dRJetPhoSubLead = deltaR(jet->eta(), jet->phi(), dipho->subLeadingPhoton()->eta(), dipho->subLeadingPhoton()->phi());

            if( dRJetPhoLead < 0.5 || dRJetPhoSubLead < 0.5 ) { continue; } // ?? can change to 0.4???
            if( !jet->passesJetID( flashgg::Loose ) ) continue;// pass jet id (reject surios detector noise)
            if( !jet->passesPuJetId(diphotons->ptrAt( candIndex ))){ continue;} // pass PU jet id
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
            
            evInfo.jet_pt.push_back(jet->pt());
            evInfo.jet_eta.push_back(jet->eta());
            evInfo.jet_phi.push_back(jet->phi());
            evInfo.jet_bdiscriminant.push_back(jet->bDiscriminator( bTag_ ));
            evInfo.jet_isMatchedToGen.push_back(isMatchedToGen);
        }
        
        // -- leptons (e, mu)
        for (UInt_t iele = 0 ; iele < electrons->size(); iele++){
            edm::Ptr<flashgg::Electron> electron = electrons->ptrAt( iele );
            if (fabs(electron->eta()) > 2.4) { continue; }
            if (electron->pt()  < electronPtThreshold_) { continue; }
            if( electron->hasMatchedConversion() ) { continue; } // remove conversions
            // missing hits: from cut-based selection: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            if( electron->isEB() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 2 ) { continue; } 
            if( electron->isEE() && electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ) > 1 ) { continue; }

            Ptr<reco::Vertex> ele_vtx = chooseElectronVertex( electron,  vertices->ptrs() );
            float d0 = electron->gsfTrack()->dxy( ele_vtx->position() );
            float dz = electron->gsfTrack()->dz( ele_vtx->position() );
            float isol = -99.; /// fix isolation computation

            evInfo.ele_pt.push_back(electron->pt());
            evInfo.ele_eta.push_back(electron->eta());
            evInfo.ele_phi.push_back(electron->phi());
            evInfo.ele_idmva.push_back(electron->nonTrigMVA());
            evInfo.ele_iso.push_back(isol);
            evInfo.ele_d0.push_back(d0);
            evInfo.ele_dz.push_back(dz);
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

            evInfo.mu_pt.push_back(muon->pt());
            evInfo.mu_eta.push_back(muon->eta());
            evInfo.mu_phi.push_back(muon->phi());
            evInfo.mu_iso.push_back(muPFCombRelIso);
            evInfo.mu_isTight.push_back(muon::isTightMuon( *muon, *muonVtx ));
            evInfo.mu_isMedium.push_back(muon::isMediumMuon( *muon ));
            evInfo.mu_isLoose.push_back(muon::isLooseMuon( *muon ));
        }
        
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
  // per-event tree
  eventTree = fs_->make<TTree>( "event", "event" );
  eventTree->Branch( "weight", &evInfo.weight, "weight/F" );
  eventTree->Branch( "npu", &evInfo.npu, "npu/I" );
  eventTree->Branch( "nvtx", &evInfo.nvtx, "nvtx/I" );

  eventTree->Branch( "pho1_pt", &evInfo.pho1_pt, "pho1_pt/F" );
  eventTree->Branch( "pho1_eta", &evInfo.pho1_eta, "pho1_eta/F" );
  eventTree->Branch( "pho1_phi", &evInfo.pho1_phi, "pho1_phi/F" );
  eventTree->Branch( "pho1_idmva", &evInfo.pho1_idmva, "pho1_idmva/F" );

  eventTree->Branch( "pho2_pt", &evInfo.pho2_pt, "pho2_pt/F" );
  eventTree->Branch( "pho2_eta", &evInfo.pho2_eta, "pho2_eta/F" );
  eventTree->Branch( "pho2_phi", &evInfo.pho2_phi, "pho2_phi/F" );
  eventTree->Branch( "pho2_idmva", &evInfo.pho1_idmva, "pho2_idmva/F" );

  eventTree->Branch( "dipho_pt", &evInfo.dipho_pt, "dipho_pt/F" );
  eventTree->Branch( "dipho_m", &evInfo.dipho_m, "dipho_m/F" );
  eventTree->Branch( "dipho_mva", &evInfo.dipho_mva, "dipho_mva/F" );

  eventTree->Branch( "jet_pt", &evInfo.jet_pt);
  eventTree->Branch( "jet_eta", &evInfo.jet_eta);
  eventTree->Branch( "jet_phi", &evInfo.jet_phi);
  eventTree->Branch( "jet_pujetid", &evInfo.jet_pujetid);
  eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant);
  eventTree->Branch( "jet_isMatchedToGen", &evInfo.jet_isMatchedToGen);

  eventTree->Branch( "ele_pt", &evInfo.ele_pt);
  eventTree->Branch( "ele_eta", &evInfo.ele_eta);
  eventTree->Branch( "ele_phi", &evInfo.ele_phi);
  eventTree->Branch( "ele_idmva", &evInfo.ele_idmva);
  eventTree->Branch( "ele_iso", &evInfo.ele_iso);
  eventTree->Branch( "ele_dz", &evInfo.ele_dz);
  eventTree->Branch( "ele_d0", &evInfo.ele_d0);

  eventTree->Branch( "mu_pt", &evInfo.mu_pt);
  eventTree->Branch( "mu_eta", &evInfo.mu_eta);
  eventTree->Branch( "mu_phi", &evInfo.mu_phi);
  eventTree->Branch( "mu_iso", &evInfo.mu_iso);
  eventTree->Branch( "mu_isTight", &evInfo.mu_isTight);
  eventTree->Branch( "mu_isMedium", &evInfo.mu_isMedium);
  eventTree->Branch( "mu_isLoose", &evInfo.mu_isLoose);

}
// ******************************************************************************************


// ******************************************************************************************
void
tthOptimizationTreeMaker::endJob()
{

} // end of endJob
// ******************************************************************************************


// ******************************************************************************************
void
tthOptimizationTreeMaker::initEventStructure()
{
    // per-event tree:
    evInfo.weight = -999.;
    evInfo.npu = -999;
    evInfo.nvtx = -999;
    
    evInfo.pho1_pt  = -999.;
    evInfo.pho1_eta = -999.;
    evInfo.pho1_phi = -999.;
    evInfo.pho1_idmva = -999.;
    
    evInfo.pho2_pt  = -999.;
    evInfo.pho2_eta = -999.;
    evInfo.pho2_phi = -999.;
    evInfo.pho2_idmva = -999.;
    
    evInfo.dipho_pt   = -999.;
    evInfo.dipho_m    = -999.;
    evInfo.dipho_mva  = -999.;
    
    evInfo.jet_pt .clear();
    evInfo.jet_eta .clear();
    evInfo.jet_phi .clear();
    evInfo.jet_pujetid .clear();
    evInfo.jet_bdiscriminant .clear();
    evInfo.jet_isMatchedToGen .clear();
    
    evInfo.ele_pt .clear();
    evInfo.ele_eta .clear();
    evInfo.ele_phi .clear();
    evInfo.ele_idmva .clear();
    evInfo.ele_iso .clear();
    evInfo.ele_dz .clear();
    evInfo.ele_d0 .clear();

    evInfo.mu_pt .clear();
    evInfo.mu_eta .clear();
    evInfo.mu_phi .clear();
    evInfo.mu_iso .clear();
    evInfo.mu_isTight .clear();
    evInfo.mu_isMedium .clear();
    evInfo.mu_isLoose .clear();

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

