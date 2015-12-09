import FWCore.ParameterSet.Config as cms

scaleBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","r9"),
    bins = cms.VPSet(
            # scale for prompt 2015, values OK, errors as in Run1 (since Run2 stat. only)
            # see Fasanella et al., ECAL DPG 03/12/2015, https://cern.ch/go/zFH8
            # w.r.t. the numbers in slide 3, the scale is computed as 1. / scale(MC) - 1.
                     cms.PSet( lowBounds = cms.vdouble(0.000,0.940), upBounds = cms.vdouble(1.000,999.000), values = cms.vdouble( 0.00118139404497 ), uncertainties = cms.vdouble( 0.00050 ) ),
                     cms.PSet( lowBounds = cms.vdouble(0.000,-999.000), upBounds = cms.vdouble(1.000,0.940), values = cms.vdouble( 0.00458088885317 ), uncertainties = cms.vdouble( 0.00050 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.000,0.940), upBounds = cms.vdouble(1.500,999.000), values = cms.vdouble( -0.00645802285147 ), uncertainties = cms.vdouble( 0.00060 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.000,-999.000), upBounds = cms.vdouble(1.500,0.940), values = cms.vdouble( 0.00339146314543 ), uncertainties = cms.vdouble( 0.00120 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.500,-999.000), upBounds = cms.vdouble(2.000,0.940), values = cms.vdouble( 0.0138594588018 ), uncertainties = cms.vdouble( 0.00200 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.500,0.940), upBounds = cms.vdouble(2.000,999.000), values = cms.vdouble( 0.00466162996303 ), uncertainties = cms.vdouble( 0.00300 ) ),
                     cms.PSet( lowBounds = cms.vdouble(2.000,-999.000), upBounds = cms.vdouble(3.000,0.940), values = cms.vdouble( 0.021878416906 ), uncertainties = cms.vdouble( 0.00200 ) ),
                     cms.PSet( lowBounds = cms.vdouble(2.000,0.940), upBounds = cms.vdouble(3.000,999.000), values = cms.vdouble( 0.014538334331 ), uncertainties = cms.vdouble( 0.00300 ) ),
                    ))

flashggDiPhotonScale = cms.EDProducer('FlashggDiPhotonSystematicProducer',
                                      src = cms.InputTag("flashggFinalEGamma","finalDiPhotons"),
                                      SystMethods2D = cms.VPSet(),
                                      # the number of syst methods matches the number of nuisance parameters
                                      # assumed for a given systematic uncertainty and is NOT required
                                      # to match 1-to-1 the number of bins above.
                                      SystMethods = cms.VPSet(
        cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
                  MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                  Label = cms.string("MCScaleHighR9EB"),
                  NSigmas = cms.vint32(),
                  OverallRange = cms.string("r9>0.94&&abs(eta)<1.5"),
                  BinList = scaleBins,
                  Debug = cms.untracked.bool(False)
                  ),
        cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
                  MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                  Label = cms.string("MCScaleLowR9EB"),
                  NSigmas = cms.vint32(),
                  OverallRange = cms.string("r9<0.94&&abs(eta)<1.5"),
                  BinList = scaleBins,
                  Debug = cms.untracked.bool(False)
                  ),
        cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
                  MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                  Label = cms.string("MCScaleHighR9EE"),
                  NSigmas = cms.vint32(),
                  OverallRange = cms.string("r9>0.94&&abs(eta)>=1.5"),
                  BinList = scaleBins,
                  Debug = cms.untracked.bool(False)
                  ),
        cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
                  MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                  Label = cms.string("MCScaleLowR9EE"),
                  NSigmas = cms.vint32(),
                  OverallRange = cms.string("r9<0.94&&abs(eta)>=1.5"),
                  BinList = scaleBins,
                  Debug = cms.untracked.bool(False)
                  ),
        
        )
                                      )
