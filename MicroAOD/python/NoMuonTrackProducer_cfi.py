import FWCore.ParameterSet.Config as cms

NoMuonTrackProducer = cms.EDProducer(
    "NoMuonTrackProducer",
    PFCandidatesTag=cms.InputTag('packedPFCandidates'),
    muonTag = cms.InputTag('slimmedMuons'), 
    )
