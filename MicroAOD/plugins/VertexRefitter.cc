// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDProducer.h"
 
///Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "flashgg/DataFormats/interface/Muon.h"

using namespace std;
using namespace edm;


//
// class declaration
//

class VertexRefitter : public edm::EDProducer 
{
public:
  explicit VertexRefitter(const edm::ParameterSet&);
  ~VertexRefitter();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  EDGetTokenT<View<pat::PackedCandidate> > pfcandidateToken_;
  EDGetTokenT<View<pat::Muon> > muonToken_;

  // ----------member data ---------------------------
  
  //double ptMin_;
};


//
// constructors and destructor
//
VertexRefitter::VertexRefitter(const ParameterSet& iConfig):
  pfcandidateToken_( consumes<View<pat::PackedCandidate> >( iConfig.getParameter<InputTag> ( "PFCandidatesTag" ) ) ),
  muonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "muonTag" ) ) )
{
  produces< reco::TrackCollection >();
}


VertexRefitter::~VertexRefitter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
VertexRefitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << " Entering VertexRefitter " << std::endl;

  Handle<View<pat::Muon> >  muons;
  iEvent.getByToken( muonToken_, muons );

  Handle<View<pat::PackedCandidate> > pfCandidates;
  iEvent.getByToken( pfcandidateToken_, pfCandidates );

  /*** get beamspot and TransientTrackBuilder from the event/eventSetup ***/
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
  reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;
  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

  /*** build a new tracklist ***/
  vector<reco::TransientTrack> selectedTracks;

  for( unsigned int i = 0 ; i < pfCandidates->size() ; i++ ) {

    Ptr<pat::PackedCandidate> cand = pfCandidates->ptrAt( i );
    if( cand->charge() == 0 ) { continue; } // skip neutrals     

    /*** use all others for fitting, needs transient tracks ***/
    reco::Track trk = reco::Track(cand->pseudoTrack());
    reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(&trk);
    transientTrack.setBeamSpot(vertexBeamSpot);
    selectedTracks.push_back(transientTrack);
  }

  //  std::vector<TransientVertex> pvs;

  /*** fit the vertex with the selected tracks ***/
  //  if( selectedTracks.size()>1 ){
  // AdaptiveVertexFitter  theFitter;
  // //TransientVertex myVertex = theFitter.vertex(mytracks, vertexBeamSpot);  // if you want the beam constraint
  //TransientVertex myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
  // // now you have a new vertex, can e.g. be compared with the original
  //std::cout << "TEST   " << myVertex->position().z() << " " << myVertex.position().z() << std::endl;
  //}else{
  //std::cout << "not enough tracks left" <<std::endl;
  //}
 
  
  //iEvent.put(AllTracks);  
  //iEvent.put(MuonLessTracks);  
}

// ------------ method called once each job just before starting event loop  ------------
void 
VertexRefitter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexRefitter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexRefitter);
