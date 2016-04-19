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
#include "DataFormats/VertexReco/interface/Vertex.h"

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
        bool passCutBasedID(edm::Ptr<flashgg::Electron> electron, const std::vector<edm::Ptr<reco::Vertex> > &pvPointers, string wp);

        EDGetTokenT<View<Electron> > electronToken_;
        EDGetTokenT<View<Photon> > photonToken_;
        EDGetTokenT<View<pat::MET> > METToken_;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;

        // thresholds for selections
        double minElectronPt_;
        double maxElectronEta_;
        double minMet_;
        string electronIdWP_;
      
    };

    WenuCandidateProducer::WenuCandidateProducer( const ParameterSet &iConfig ) :
        electronToken_( consumes<View<flashgg::Electron> >( iConfig.getParameter<InputTag>( "ElectronTag" ) ) ),
        photonToken_( consumes<View<flashgg::Photon> >( iConfig.getParameter<InputTag>( "PhotonTag" ) ) ),
        METToken_( consumes<View<pat::MET> >( iConfig.getParameter<InputTag> ( "METTag" ) ) ),
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) )
    {

        minElectronPt_  = iConfig.getParameter<double>( "minElectronPt");
        maxElectronEta_ = iConfig.getParameter<double>( "maxElectronEta");
        minMet_  = iConfig.getParameter<double>( "minMet");
        electronIdWP_  = iConfig.getParameter<string>( "electronIdWP");

        produces<vector<WenuCandidate> >();
    }

    
    bool WenuCandidateProducer::passCutBasedID(edm::Ptr<flashgg::Electron> electron, const std::vector<edm::Ptr<reco::Vertex> > &pvPointers, string wp)
    {

        bool pass = false;

        float full5x5_sigmaIetaIeta = electron->full5x5_sigmaIetaIeta();
        float dEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
        float dPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
        float hOverE = electron->hcalOverEcal();
        float relIsoEA = electron->standardHggIso()/electron->pt();

        float ooEmooP =-999 ; 

        if( electron->ecalEnergy() == 0 ){
            ooEmooP = 1e30;
        }else if( !std::isfinite(electron->ecalEnergy())){    
            ooEmooP = 1e30;
        }else{
            ooEmooP = fabs(1.0/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy() );
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
        float dxy = fabs( electron->gsfTrack()->dxy( best_vtx_elec->position()) ) ;
        float dz = fabs( electron->gsfTrack()->dz( best_vtx_elec->position())) ;

        int missedHits = electron->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS);

        bool isEB = (fabs(electron->superCluster()->eta())<1.479);

        if (isEB){
            if ( wp == "loose" 
                 && full5x5_sigmaIetaIeta < 0.0103
                 && fabs(dEtaIn) < 0.0105
                 && fabs(dPhiIn) < 0.115
                 && hOverE < 0.104
                 && relIsoEA < 0.0893
                 && ooEmooP < 0.102
                 && fabs(dxy) < 0.0261
                 && fabs(dz) < 0.41
                 && missedHits <=2 ) 
                pass =  true;

            if ( wp == "medium" 
                 && full5x5_sigmaIetaIeta < 0.0101
                 && fabs(dEtaIn) < 0.0103
                 && fabs(dPhiIn) < 0.0336
                 && hOverE < 0.0876
                 && relIsoEA < 0.0766
                 && ooEmooP < 0.0174
                 && fabs(dxy) < 0.0118
                 && fabs(dz) < 0.373
                 && missedHits <=2 ) 
                pass =  true;

            if ( wp == "tight" 
                 && full5x5_sigmaIetaIeta < 0.0101
                 && fabs(dEtaIn) < 0.00926
                 && fabs(dPhiIn) < 0.0336
                 && hOverE < 0.0597
                 && relIsoEA < 0.0354
                 && ooEmooP < 0.012
                 && fabs(dxy) < 0.0111
                 && fabs(dz) < 0.0466
                 && missedHits <=2 ) 
                pass =  true;
        }
        else {
            if ( wp == "loose"
                 && full5x5_sigmaIetaIeta < 0.0301
                 && fabs(dEtaIn) < 0.00814
                 && fabs(dPhiIn) < 0.182
                 && hOverE < 0.0897
                 && relIsoEA < 0.121
                 && ooEmooP < 0.126
                 && fabs(dxy) < 0.118
                 && fabs(dz) < 0.822
                 && missedHits <=1 )
                pass = true; 

            if ( wp == "medium"
                 && full5x5_sigmaIetaIeta < 0.0283
                 && fabs(dEtaIn) < 0.00733
                 && fabs(dPhiIn) < 0.114
                 && hOverE < 0.0678
                 && relIsoEA < 0.0678
                 && ooEmooP < 0.0898
                 && fabs(dxy) < 0.0739
                 && fabs(dz) < 0.602
                 && missedHits <=1 )
                pass = true; 
            
            if ( wp == "tight"
                 && full5x5_sigmaIetaIeta < 0.0279
                 && fabs(dEtaIn) < 0.00724
                 && fabs(dPhiIn) < 0.0918
                 && hOverE < 0.0615
                 && relIsoEA < 0.0646
                 && ooEmooP < 0.00999
                 && fabs(dxy) < 0.0351
                 && fabs(dz) < 0.417
                 && missedHits <=1 )
                pass = true; 

        }

        return pass;

    
        
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


        Handle<View<reco::Vertex> > vertices;
        evt.getByToken( vertexToken_, vertices );

        std::auto_ptr<vector<WenuCandidate> > WenuCandidates( new vector<WenuCandidate> );
        
        for( unsigned int iele = 0; iele < electrons->size(); iele++) {

            Ptr<flashgg::Electron> theElectron = electrons->ptrAt(iele);

            // --- check if the electron passes pt, eta cuts and electron ID
            if (theElectron->pt() < minElectronPt_) continue;
            if (fabs(theElectron->eta()) > maxElectronEta_) continue;
            //if (!wpXX) continue;
            if( theElectron->hasMatchedConversion() ) continue; 
            if ( !passCutBasedID(theElectron, vertices->ptrs(), electronIdWP_) ) continue;

            // --- find photons corresponding to the electrons
            int imatch = -1;
            for (unsigned int ipho = 0; ipho < photons->size(); ipho++){
                
                if( &( *photons->ptrAt( ipho )->superCluster() ) == &( *theElectron ->superCluster() ) ) {
                    imatch = ipho;
                    break;
                }
            }


            // min met cut
            if (thMet->pt() < minMet_) continue;
            
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

