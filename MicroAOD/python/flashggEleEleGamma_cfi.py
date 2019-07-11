import FWCore.ParameterSet.Config as cms

flashggEleEleGamma = cms.EDProducer('FlashggEleEleGammaProducer',
                                  DiElectronTag=cms.InputTag('flashggDiElectrons'),
                                  PhotonTag=cms.InputTag('flashggPhotons'),
                                  VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                  ##Parameters                                                
                                  minPhotonPT=cms.double(10.)
                                  )
