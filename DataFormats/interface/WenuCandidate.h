#ifndef FLASHgg_WenuCandidate_h
#define FLASHgg_WenuCandidate_h

#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

namespace flashgg {

  class WenuCandidate: public reco::CompositeCandidate
  {
  public:
    WenuCandidate();
    WenuCandidate( edm::Ptr<Electron> electron, edm::Ptr<Photon> photon, edm::Ptr<pat::MET> met  );

    ~WenuCandidate();

    WenuCandidate *clone() const { return ( new WenuCandidate( *this ) ); }

    const edm::Ptr<flashgg::Electron>  electron() const {return electron_;}
    const edm::Ptr<flashgg::Photon> photon() const {return photon_;}
    const edm::Ptr<pat::MET> met() const {return met_;}

    void setElectron( edm::Ptr<Electron> electron ) {electron_ = electron;}
    void setPhoton( edm::Ptr<Photon> photon ) {photon_ = photon;}
    void setMet( edm::Ptr<pat::MET> met ) {met_ = met;}



  private:
      edm::Ptr<Electron>  electron_;
      edm::Ptr<Photon> photon_;
      edm::Ptr<pat::MET> met_;
  };
}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
