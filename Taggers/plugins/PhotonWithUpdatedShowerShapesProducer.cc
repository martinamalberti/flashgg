#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"

#include "TFile.h"
#include "TGraph.h"

namespace flashgg {

    class PhotonWithUpdatedShowerShapesProducer : public edm::EDProducer
    {
    public:
        PhotonWithUpdatedShowerShapesProducer( const edm::ParameterSet & );
        void produce( edm::Event &, const edm::EventSetup & ) override;

    private:
        edm::EDGetTokenT<edm::View<flashgg::Photon> > token_;
        edm::FileInPath correctionFile_;
        bool correctInputs_;
        bool debug_;
        //        std::vector<TGraph*> corrections_;
        std::vector<std::unique_ptr<TGraph> > corrections_;
    };

    PhotonWithUpdatedShowerShapesProducer::PhotonWithUpdatedShowerShapesProducer( const edm::ParameterSet &ps ) :
        token_(consumes<edm::View<flashgg::Photon> >(ps.getParameter<edm::InputTag>("src"))),
        debug_( ps.getParameter<bool>( "Debug" ) )
    {
        correctInputs_ = ps.existsAs<edm::FileInPath>("correctionFile") ? true: false;
        if (correctInputs_) {
            correctionFile_ = ps.getParameter<edm::FileInPath>( "correctionFile" );
            TFile* f = TFile::Open(correctionFile_.fullPath().c_str());
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EB"))->Clone() );
            f->Close();
        }

        produces<std::vector<flashgg::Photon> >();
    }

    void PhotonWithUpdatedShowerShapesProducer::produce( edm::Event &evt, const edm::EventSetup & )
    {
        edm::Handle<edm::View<flashgg::Photon> > objects;
        evt.getByToken( token_, objects );

        auto_ptr<std::vector<flashgg::Photon> > out_obj( new std::vector<flashgg::Photon>() );

        for (const auto & obj : *objects) {
            flashgg::Photon *new_obj = obj.clone();
            //            new_obj->makePhotonsPersistent();
            if (not evt.isRealData() and correctInputs_) { 
                if (new_obj->isEB()) {
                    if (this->debug_) {
                        std::cout << "Original shower shapes values " << std::endl;
                        std::cout << new_obj->full5x5_r9() << std::endl;
                        std::cout << new_obj->s4() << std::endl;
                    }
                    reco::Photon::ShowerShape newShowerShapes = new_obj->full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->full5x5_r9())*new_obj->superCluster()->rawEnergy();
                    new_obj->full5x5_setShowerShapeVariables(newShowerShapes);
                    new_obj->setS4(corrections_[2]->Eval(new_obj->s4()));

                    if (this->debug_) {
                        std::cout << "Corrected shower shapes values " << std::endl;
                        std::cout << new_obj->full5x5_r9() << std::endl;
                        std::cout << new_obj->s4() << std::endl;
                    }
                }
            }

            out_obj->push_back(*new_obj);
            delete new_obj;
        }
        evt.put(out_obj);
    }
}

typedef flashgg::PhotonWithUpdatedShowerShapesProducer FlashggPhotonWithUpdatedShowerShapesProducer;
DEFINE_FWK_MODULE( FlashggPhotonWithUpdatedShowerShapesProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
