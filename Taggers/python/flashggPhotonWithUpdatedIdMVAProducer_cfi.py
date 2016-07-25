import FWCore.ParameterSet.Config as cms

flashggUpdatedIdMVAPhotons = cms.EDProducer("FlashggPhotonWithUpdatedIdMVAProducer",
                                              src                      = cms.InputTag("flashggRandomizedPhotons"),
                                              rhoFixedGridCollection   = cms.InputTag('fixedGridRhoAll'),
                                              VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),     
                                              photonIdMVAweightfile_EB = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_80X_barrel_ICHEPvtx.xml"),
                                              photonIdMVAweightfile_EE = cms.FileInPath("flashgg/MicroAOD/data/MVAweights_80X_endcap_ICHEPvtx.xml"),
                                              # commenting out this parameter will disable all corrections performed by this module
                                              correctionFile           = cms.FileInPath("flashgg/MicroAOD/data/transformation_80X_v2.root"),
                                              Debug                    = cms.bool(True)
                                              )
