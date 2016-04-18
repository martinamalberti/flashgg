#include "flashgg/Taggers/interface/WenuDumper.h"

#include "flashgg/Taggers/interface/PluggableAnalyzer.h"

namespace flashgg {

    namespace fwlite {
        PLUGGABLE_ANALYZER( WenuDumper );
        PLUGGABLE_ANALYZER( CutBasedWenuDumper );
    }
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
