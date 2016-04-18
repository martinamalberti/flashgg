#include "flashgg/DataFormats/interface/WenuCandidate.h"

using namespace flashgg;

WenuCandidate::WenuCandidate() {}

WenuCandidate::~WenuCandidate() {}

WenuCandidate::WenuCandidate( edm::Ptr<Electron> electron, edm::Ptr<Photon> photon, edm::Ptr<pat::MET> met  )
{
  electron_ = electron;
  photon_ = photon;
  met_ = met;
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
