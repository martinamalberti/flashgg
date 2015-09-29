// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDProducer.h"
 
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "flashgg/DataFormats/interface/Muon.h"

using namespace std;
using namespace edm;



//
// class declaration
//

class NoMuonTrackProducer : public edm::EDProducer 
{
public:
  explicit NoMuonTrackProducer(const edm::ParameterSet&);
  ~NoMuonTrackProducer();
  
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
NoMuonTrackProducer::NoMuonTrackProducer(const ParameterSet& iConfig):
  pfcandidateToken_( consumes<View<pat::PackedCandidate> >( iConfig.getParameter<InputTag> ( "PFCandidatesTag" ) ) ),
  muonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "muonTag" ) ) )
{
  produces< reco::TrackCollection >();
}


NoMuonTrackProducer::~NoMuonTrackProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
NoMuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << " Entering NoMuonTrackProducer " << std::endl;

  Handle<View<pat::Muon> >  muons;
  iEvent.getByToken( muonToken_, muons );

  Handle<View<pat::PackedCandidate> > pfCandidates;
  iEvent.getByToken( pfcandidateToken_, pfCandidates );

  std::auto_ptr<reco::TrackCollection> AllTracks(new reco::TrackCollection);
  std::auto_ptr<reco::TrackCollection> MuonLessTracks(new reco::TrackCollection);

  for( unsigned int i = 0 ; i < pfCandidates->size() ; i++ ) {

    Ptr<pat::PackedCandidate> cand = pfCandidates->ptrAt( i );
    if( cand->charge() == 0 ) { continue; } // skip neutrals

    reco::Track trk = reco::Track(cand->pseudoTrack());

    AllTracks->push_back(cand->pseudoTrack());

    if (!cand->isMuon()) 
      MuonLessTracks->push_back(cand->pseudoTrack());
 
  }
  
  iEvent.put(AllTracks);  
  //iEvent.put(MuonLessTracks);  
}

// ------------ method called once each job just before starting event loop  ------------
void 
NoMuonTrackProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NoMuonTrackProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(NoMuonTrackProducer);
