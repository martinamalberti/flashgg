import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggTkVtxMap_cfi import flashggVertexMapUnique,flashggVertexMapNonUnique
from flashgg.MicroAOD.flashggPhotons_cfi import flashggPhotons
from flashgg.MicroAOD.flashggPhotonJet_cfi import flashggPhotonJet
from flashgg.MicroAOD.flashggMicroAODGenSequence_cff import *

eventCount = cms.EDProducer("EventCountProducer")
weightsCount = cms.EDProducer("WeightsCountProducer",
                              generator=cms.InputTag("generator"),
                              pileupInfo=cms.InputTag("addPileupInfo"),
                              doObsPileup=cms.untracked.bool(True),
                              minObsPileup=cms.double(-0.5),
                              maxObsPileup=cms.double(100.5),
                              nbinsObsPileup=cms.int32(101),
                              )

flashggMicroAODPhotonJetValidationSequence = cms.Sequence(eventCount+weightsCount
                                                          +flashggVertexMapUnique+flashggVertexMapNonUnique
                                                          +flashggMicroAODGenSequence
                                                          +flashggPhotons
                                                          #+flashggFinalEGamma
                                                          +flashggPhotonJet
                                                          )

