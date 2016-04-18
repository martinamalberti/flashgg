#ifndef flashgg_WenuDumper_h
#define flashgg_WenuDumper_h

#include "flashgg/DataFormats/interface/WenuCandidate.h"
#include "flashgg/Taggers/interface/CollectionDumper.h"


namespace flashgg {

    typedef CollectionDumper<std::vector<WenuCandidate> > WenuDumper;
    typedef CollectionDumper<std::vector<WenuCandidate>,
            WenuCandidate,
            CutBasedClassifier<WenuCandidate> > CutBasedWenuDumper;

}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
