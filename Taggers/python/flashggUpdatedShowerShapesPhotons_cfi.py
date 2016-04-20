import FWCore.ParameterSet.Config as cms

flashggUpdatedShowerShapesPhotons = cms.EDProducer("FlashggPhotonWithUpdatedShowerShapesProducer",
                                              src                      = cms.InputTag("flashggRandomizedPhotons"),
                                              correctionFile           = cms.FileInPath("flashgg/MicroAOD/data/transformation_76X_v2.root"),
                                              Debug                    = cms.bool(False)
                                              )
