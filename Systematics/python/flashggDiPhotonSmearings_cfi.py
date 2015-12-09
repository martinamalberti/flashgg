import FWCore.ParameterSet.Config as cms

smearBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","r9"),
    bins = cms.VPSet(
            # smearings for prompt 2015, values OK, errors as statistical of Run2 (since Run1 different parametrization)
            # see Fasanella et al., ECAL DPG 03/12/2015, https://cern.ch/go/zFH8
                     cms.PSet( lowBounds = cms.vdouble(0.000,-999.000), upBounds = cms.vdouble(1.000,0.940), values = cms.vdouble( 0.013654 ), uncertainties = cms.vdouble( 0.00024652 ) ),
                     cms.PSet( lowBounds = cms.vdouble(0.000,0.940), upBounds = cms.vdouble(1.000,999.000), values = cms.vdouble( 0.014142 ), uncertainties = cms.vdouble( 0.00018319 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.000,-999.000), upBounds = cms.vdouble(1.500,0.940), values = cms.vdouble( 0.020859 ), uncertainties = cms.vdouble( 0.00032435 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.000,0.940), upBounds = cms.vdouble(1.500,999.000), values = cms.vdouble( 0.017120 ), uncertainties = cms.vdouble( 0.00098551 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.500,-999.000), upBounds = cms.vdouble(2.000,0.940), values = cms.vdouble( 0.028083 ), uncertainties = cms.vdouble( 0.00045816 ) ),
                     cms.PSet( lowBounds = cms.vdouble(1.500,0.940), upBounds = cms.vdouble(2.000,999.000), values = cms.vdouble( 0.027289 ), uncertainties = cms.vdouble( 0.00076386 ) ),
                     cms.PSet( lowBounds = cms.vdouble(2.000,-999.000), upBounds = cms.vdouble(3.000,0.940), values = cms.vdouble( 0.031793 ), uncertainties = cms.vdouble( 0.00061517 ) ),
                     cms.PSet( lowBounds = cms.vdouble(2.000,0.940), upBounds = cms.vdouble(3.000,999.000), values = cms.vdouble( 0.030831 ), uncertainties = cms.vdouble( 0.00042066 ) ),
                    ))



flashggDiPhotonSmearings = cms.EDProducer('FlashggDiPhotonSystematicProducer',
                                            src = cms.InputTag("flashggFinalEGamma","finalDiPhotons"),
                SystMethods2D = cms.VPSet(),
                # the number of syst methods matches the number of nuisance parameters
                # assumed for a given systematic uncertainty and is NOT required
                # to match 1-to-1 the number of bins above.
                SystMethods = cms.VPSet(
                  cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
                             MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                             Label = cms.string("MCSmearHighR9EE"),
                             NSigmas = cms.vint32(),
                             OverallRange = cms.string("r9>0.94&&abs(eta)>=1.5"),
                                BinList = smearBins,
                             Debug = cms.untracked.bool(False),
                             ExaggerateShiftUp = cms.untracked.bool(False),
                             ),
                   cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
                             MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                             Label = cms.string("MCSmearLowR9EE"),
                             NSigmas = cms.vint32(),
                             OverallRange = cms.string("r9<0.94&&abs(eta)>=1.5"),
                             BinList = smearBins,
                             Debug = cms.untracked.bool(False),
                             ExaggerateShiftUp = cms.untracked.bool(False),
                             ),
                   cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
                             MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                             Label = cms.string("MCSmearHighR9EB"),
                             NSigmas = cms.vint32(),
                             OverallRange = cms.string("r9>0.94&&abs(eta)<1.5"),
                             BinList = smearBins,
                             Debug = cms.untracked.bool(False),
                             ExaggerateShiftUp = cms.untracked.bool(False),
                             ),
                   cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
                             MethodName = cms.string("FlashggDiPhotonFromPhoton"),
                             Label = cms.string("MCSmearLowR9EB"),
                             NSigmas = cms.vint32(),
                             OverallRange = cms.string("r9<=0.94&&abs(eta)<1.5"),
                             BinList = smearBins,
                             Debug = cms.untracked.bool(False),
                             ExaggerateShiftUp = cms.untracked.bool(False),
                             )
                   )
)
