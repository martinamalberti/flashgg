#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TFile.h"
#include "TGraph.h"

namespace flashgg {

  class PhotonWithUpdatedIdMVAProducer : public edm::EDProducer
  {
  public:
    PhotonWithUpdatedIdMVAProducer( const edm::ParameterSet & );
    void produce( edm::Event &, const edm::EventSetup & ) override;

  private:
    edm::EDGetTokenT<edm::View<flashgg::Photon> > token_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
    PhotonIdUtils phoTools_;
    edm::FileInPath phoIdMVAweightfileEB_, phoIdMVAweightfileEE_, correctionFile_;
    bool correctInputs_;
    bool debug_;
    //        std::vector<TGraph*> corrections_;                                                                                                                                               
    std::vector<std::unique_ptr<TGraph> > corrections_;

  };

  PhotonWithUpdatedIdMVAProducer::PhotonWithUpdatedIdMVAProducer( const edm::ParameterSet &ps ) :
    token_(consumes<edm::View<flashgg::Photon> >(ps.getParameter<edm::InputTag>("src"))),
    rhoToken_( consumes<double>( ps.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
    vertexToken_( consumes<edm::View<reco::Vertex> >( ps.getParameter<edm::InputTag> ( "VertexTag" ) ) ),
    debug_( ps.getParameter<bool>( "Debug" ) )
  {
    phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB" );
    phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE" );
    phoTools_.setupMVA( phoIdMVAweightfileEB_.fullPath(), phoIdMVAweightfileEE_.fullPath() );

    correctInputs_ = ps.existsAs<edm::FileInPath>("correctionFile") ? true: false;
    if (correctInputs_) {
      correctionFile_ = ps.getParameter<edm::FileInPath>( "correctionFile" );
      TFile* f = TFile::Open(correctionFile_.fullPath().c_str());
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EB"))->Clone() );
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEB"))->Clone() );
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EB"))->Clone() );
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EE"))->Clone() );
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEE"))->Clone() );
      corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EE"))->Clone() );
      f->Close();
    }
    
    produces<std::vector<flashgg::Photon> >();
  }

  void PhotonWithUpdatedIdMVAProducer::produce( edm::Event &evt, const edm::EventSetup & )
  {
      std::cout << " ciao " <<std::endl;
      
    edm::Handle<edm::View<flashgg::Photon> > objects;
    evt.getByToken( token_, objects );

    edm::Handle<double> rhoHandle;
    evt.getByToken( rhoToken_, rhoHandle );
    const double rhoFixedGrd = *( rhoHandle.product() );

    edm::Handle<edm::View<reco::Vertex> > primaryVertices;
    evt.getByToken( vertexToken_, primaryVertices );
    const std::vector<edm::Ptr<reco::Vertex>> &pvPointers = primaryVertices->ptrs();
    edm::Ptr<reco::Vertex> pvx = pvPointers[0]; //selected vertex 0 for Zmumugamma


    auto_ptr<std::vector<flashgg::Photon> > out_obj( new std::vector<flashgg::Photon>() );

    for (const auto & obj : *objects) {
        flashgg::Photon *new_obj = obj.clone();
        //            new_obj->makePhotonsPersistent(); ???
        double correctedEtaWidth = 0;
        if (not evt.isRealData() and correctInputs_) { 
            if (new_obj->isEB()) {
                std::cout << "Correction in EB" << std::endl;
                if (this->debug_) {
                    std::cout << new_obj->full5x5_r9() << std::endl;
                    std::cout << new_obj->r9() << std::endl;
                }
                reco::Photon::ShowerShape newShowerShapes = new_obj->full5x5_showerShapeVariables();
                newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->full5x5_r9())*new_obj->superCluster()->rawEnergy();
                new_obj->full5x5_setShowerShapeVariables(newShowerShapes);
                correctedEtaWidth = corrections_[1]->Eval(new_obj->superCluster()->etaWidth());
                new_obj->getSuperCluster()->setEtaWidth(correctedEtaWidth);
                new_obj->setS4(corrections_[2]->Eval(new_obj->s4()));
                
                if (this->debug_) {
                    std::cout << new_obj->full5x5_r9() << std::endl;
                    std::cout << new_obj->r9() << std::endl;
                }
            }
            
            if (new_obj->isEE()) {
                std::cout << "Correction in EE" << std::endl;
                if (this->debug_) {
                    std::cout << new_obj->full5x5_r9() << std::endl;
                    std::cout << new_obj->r9() << std::endl;
                }
                reco::Photon::ShowerShape newShowerShapes = new_obj->full5x5_showerShapeVariables();
                newShowerShapes.e3x3 = corrections_[3]->Eval(new_obj->full5x5_r9())*new_obj->superCluster()->rawEnergy();
                new_obj->full5x5_setShowerShapeVariables(newShowerShapes);
                correctedEtaWidth = corrections_[4]->Eval(new_obj->superCluster()->etaWidth());
                //std::cout<<new_obj->superCluster()->etaWidth() << "  " << correctedEtaWidth <<std::endl;
                new_obj->getSuperCluster()->setEtaWidth(correctedEtaWidth);
                new_obj->setS4(corrections_[5]->Eval(new_obj->s4()));
                
                if (this->debug_) {
                    std::cout << new_obj->full5x5_r9() << std::endl;
                    std::cout << new_obj->r9() << std::endl;
                }
            }
        }

        std::cout << "Recomputing photon ID MVA ..." <<std::endl;

        if (this->debug_) {
            std::cout << " Output Photon lead IDMVA: " << new_obj->phoIdMvaDWrtVtx(pvx)   << std::endl;
        }
        float newleadmva = phoTools_.computeMVAWrtVtx( *new_obj, pvx, rhoFixedGrd, correctedEtaWidth );
        new_obj->setPhoIdMvaWrtVtx( pvx, newleadmva);
        if (this->debug_) {
            std::cout << " Output Photon lead IDMVA new: " << new_obj->phoIdMvaDWrtVtx(pvx)   << std::endl;
        }
        
        out_obj->push_back(*new_obj);
        delete new_obj;
    }
    evt.put(out_obj);
  }
}

typedef flashgg::PhotonWithUpdatedIdMVAProducer FlashggPhotonWithUpdatedIdMVAProducer;
DEFINE_FWK_MODULE( FlashggPhotonWithUpdatedIdMVAProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
