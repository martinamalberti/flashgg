#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
#include "flashgg/Validation/interface/SimpleTreeMaker.h"

typedef edm::AnalyzerWrapper<SimpleTreeMaker> EDSimpleTreeMaker;
DEFINE_FWK_MODULE(EDSimpleTreeMaker);
