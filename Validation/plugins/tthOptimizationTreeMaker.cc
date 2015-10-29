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

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TTree.h"

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

  float pho2_pt;
  float pho2_eta;
  float pho2_phi;

  float dipho_m;
  float dipho_pt;
  float dipho_mva;

  vector<float> jet_pt;
  vector<float> jet_eta;
  vector<float> jet_phi;
  vector<float> jet_pujetid;
  vector<float> jet_bdiscriminant;

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

    typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

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
    inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) )
{
    jetPtThreshold_ = iConfig.getUntrackedParameter<double>( "jetPtThreshold", 20. );
    bTag_ = iConfig.getUntrackedParameter<string> ( "bTag", "pfCombinedInclusiveSecondaryVertexV2BJetTags" );
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
    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<View<flashgg::DiPhotonCandidate> > diphotons;
    iEvent.getByToken( diphotonToken_, diphotons );
    
    Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;
    iEvent.getByToken( mvaResultToken_, mvaResults );
    
    JetCollectionVector Jets( inputTagJets_.size() );
    for( size_t j = 0; j < inputTagJets_.size(); ++j ) {
        iEvent.getByLabel( inputTagJets_[j], Jets[j] );
    }
    
    // initialize tree
    initEventStructure();
    
    // fill variables for photons, jets, etcc...
    
    
    if (diphotons->size() > 0) {
        // event weight
        // boh???
        
        // number of vertices
        evInfo.nvtx = primaryVertices->size() ;

        // -- photons
        // take the diphoton candidate with index = 0 (--> highest pt1+pt2)... to be checked..
        int candIndex = 0;
        edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( candIndex );
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( candIndex );
        
        evInfo.pho1_pt  = dipho->leadingPhoton()->pt();
        evInfo.pho1_eta = dipho->leadingPhoton()->eta();
        evInfo.pho1_phi = dipho->leadingPhoton()->phi();
        
        evInfo.pho2_pt  = dipho->subLeadingPhoton()->pt();
        evInfo.pho2_eta = dipho->subLeadingPhoton()->eta();
        evInfo.pho2_phi = dipho->subLeadingPhoton()->phi();
        
        evInfo.dipho_pt  = dipho->pt();
        evInfo.dipho_m   = dipho->mass();
        evInfo.dipho_mva = mvares->result ;
        
        // -- jets
        // take the jets corresponding to the diphoton candidate
        unsigned int jetCollectionIndex = diphotons->ptrAt( candIndex )->jetCollectionIndex();
        
        int njets = 0;
        
        for( UInt_t ijet = 0; ijet < Jets[jetCollectionIndex]->size() ; ijet++ ) {
            
            Ptr<flashgg::Jet> jet  = Jets[jetCollectionIndex]->ptrAt( ijet );
            
            float dEtaLead = jet->eta() - dipho->leadingPhoton()->eta();
            float dEtaSublead = jet->eta() - dipho->subLeadingPhoton()->eta();
            
            float dPhiLead = deltaPhi( jet->phi(), dipho->leadingPhoton()->phi() );
            float dPhiSublead = deltaPhi( jet->phi(), dipho->subLeadingPhoton()->phi() );
            
            float dRJetPhoLead = sqrt( dEtaLead * dEtaLead + dPhiLead * dPhiLead );
            float dRJetPhoSubLead = sqrt( dEtaSublead * dEtaSublead + dPhiSublead * dPhiSublead );
            
            if( dRJetPhoLead < 0.5 || dRJetPhoSubLead < 0.5 ) { continue; } // ?? can change to 0.4???
            if( jet->pt() < jetPtThreshold_ ) { continue; }
            
            njets++;
            evInfo.jet_pt.push_back(jet->pt());
            evInfo.jet_eta.push_back(jet->eta());
            evInfo.jet_phi.push_back(jet->phi());
            //evInfo.jet_pujetid.push_back( ?? );
            evInfo.jet_bdiscriminant.push_back(jet->bDiscriminator( bTag_ ));
        }      
        
        // fill the tree
        if ( njets > 0. )
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
  eventTree->Branch( "pho2_pt", &evInfo.pho2_pt, "pho2_pt/F" );
  eventTree->Branch( "pho2_eta", &evInfo.pho2_eta, "pho2_eta/F" );
  eventTree->Branch( "pho2_phi", &evInfo.pho2_phi, "pho2_phi/F" );
  eventTree->Branch( "dipho_pt", &evInfo.dipho_pt, "dipho_pt/F" );
  eventTree->Branch( "dipho_m", &evInfo.dipho_m, "dipho_m/F" );
  eventTree->Branch( "dipho_mva", &evInfo.dipho_mva, "dipho_mva/F" );

  eventTree->Branch( "jet_pt", &evInfo.jet_pt);
  eventTree->Branch( "jet_eta", &evInfo.jet_eta);
  eventTree->Branch( "jet_phi", &evInfo.jet_phi);
  eventTree->Branch( "jet_pujetid", &evInfo.jet_pujetid);
  eventTree->Branch( "jet_bdiscriminant", &evInfo.jet_bdiscriminant);

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
  evInfo.pho2_pt  = -999.;
  evInfo.pho2_eta = -999.;
  evInfo.pho2_phi = -999.;
  
  evInfo.dipho_pt   = -999.;
  evInfo.dipho_m    = -999.;
  evInfo.dipho_mva  = -999.;

  evInfo.jet_pt .clear();
  evInfo.jet_eta .clear();
  evInfo.jet_phi .clear();
  evInfo.jet_pujetid .clear();
  evInfo.jet_bdiscriminant .clear();
};
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

