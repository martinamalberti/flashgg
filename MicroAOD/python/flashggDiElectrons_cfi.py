import FWCore.ParameterSet.Config as cms

flashggDiElectrons = cms.EDProducer('FlashggDiElectronProducer',
                                  ElectronTag=cms.InputTag('flashggRandomizedElectrons'),
                                  ##Parameters                                                
                                  minElectronPt=cms.double(7.),
                                  maxElectronEta=cms.double(2.4)
                                  )
