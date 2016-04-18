import FWCore.ParameterSet.Config as cms


flashggWenu = cms.EDProducer('FlashggWenuCandidateProducer',
                             PhotonTag=cms.InputTag('flashggPhotons'), 
                             ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                             METTag=cms.InputTag('slimmedMETs'),   
                             minElectronPt = cms.double(30.),
                             maxElectronEta = cms.double(2.5),
                             minMet = cms.double(25.)
                             )

