#include "FWCore/Framework/interface/MakerMacros.h"
#include "flashgg/Taggers/interface/WenuDumper.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"

typedef edm::AnalyzerWrapper<flashgg::WenuDumper> WenuDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedWenuDumper> CutBasedWenuDumper;

DEFINE_FWK_MODULE( WenuDumper );
DEFINE_FWK_MODULE( CutBasedWenuDumper );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

