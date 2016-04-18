#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/WenuCandidate.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Photon.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include <vector>
#include <algorithm>
#include <string>
#include <utility>
#include "TLorentzVector.h"
#include "TMath.h"

using namespace std;
using namespace edm;


namespace flashgg {
    class WenuCandidateProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        WenuCandidateProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<Photon> > photonToken_;
        EDGetTokenT<View<pat::MET> > METToken_;

        // thresholds for selections
        double minElectronPt_;
        double maxElectronEta_;
        double minMet_;
        //double electronWP_;
        
    };

    WenuCandidateProducer::WenuCandidateProducer( const ParameterSet &iConfig ) :
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        photonToken_( consumes<View<flashgg::Photon> >( iConfig.getParameter<InputTag>( "PhotonTag" ) ) ),
        METToken_( consumes<View<pat::MET> >( iConfig.getParameter<InputTag> ( "METTag" ) ) )
    {

        minElectronPt_  = iConfig.getParameter<double>( "minElectronPt");
        maxElectronEta_ = iConfig.getParameter<double>( "maxElectronEta");
        minMet_  = iConfig.getParameter<double>( "minMet");

        produces<vector<WenuCandidate> >();
    }

    
    void WenuCandidateProducer::produce( Event &evt, const EventSetup & )
        
    {

        Handle<View<flashgg::Electron> > electrons;
        evt.getByToken( electronToken_, electrons );

        Handle<View<flashgg::Photon> > photons;
        evt.getByToken( photonToken_, photons );

        Handle<View<pat::MET> > METs;
        evt.getByToken( METToken_, METs );

        if( METs->size() != 1 ) { 
            std::cout << "WARNING number of MET is not equal to 1" << std::endl; 
        }
        Ptr<pat::MET> theMet = METs->ptrAt( 0 );


        std::auto_ptr<vector<WenuCandidate> > WenuCandidates( new vector<WenuCandidate> );
        
        for( unsigned int iele = 0; iele < electrons->size(); iele++) {

            Ptr<flashgg::Electron> theElectron = electrons->ptrAt(iele);

            // --- check if the electron passes pt, eta cuts and electron ID
            if (theElectron->pt() < minElectronPt_) continue;
            if (fabs(theElectron->eta()) > maxElectronEta_) continue;
            //if (!wpXX) continue;

            

            // --- find photons corresponding to the electrons
            int imatch = -1;
            for (unsigned int ipho = 0; ipho < photons->size(); ipho++){
                
                if( &( *photons->ptrAt( ipho )->superCluster() ) == &( *theElectron ->superCluster() ) ) {
                    imatch = ipho;
                    break;
                }
            }

            
            if (imatch != -1){
                Ptr<flashgg::Photon> thePhoton = photons->ptrAt(imatch);
                WenuCandidate WenuCand(theElectron, thePhoton, theMet);
                WenuCand.setElectron( theElectron  );
                WenuCand.setPhoton( thePhoton );
                WenuCand.setMet( theMet );
                WenuCandidates->push_back(WenuCand);
           }

        } // endl loop over electrons

        evt.put(WenuCandidates);

    }

}

typedef flashgg::WenuCandidateProducer FlashggWenuCandidateProducer;
DEFINE_FWK_MODULE( FlashggWenuCandidateProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

