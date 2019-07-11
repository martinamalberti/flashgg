#ifndef FLASHgg_DiElectronCandidate_h
#define FLASHgg_DiElectronCandidate_h

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

namespace flashgg {
    class DiElectronCandidate : public reco::CompositeCandidate
    {
    public:
        DiElectronCandidate();
        DiElectronCandidate( edm::Ptr<pat::Electron>, edm::Ptr<pat::Electron> );
        DiElectronCandidate( const pat::Electron &, const pat::Electron & );
        ~DiElectronCandidate();

        const pat::Electron *leadingElectron() const;
        const pat::Electron *subleadingElectron() const;

    private:


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
