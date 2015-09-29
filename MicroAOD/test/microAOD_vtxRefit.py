import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#process.GlobalTag.globaltag = 'auto:run2_mc'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 20 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

# Fix because auto:run2_mc points to MCRUN2_74_V9::All
current_gt = process.GlobalTag.globaltag.value()
if current_gt.count("::All"):
    new_gt = current_gt.replace("::All","")
    print 'Removing "::All" from GlobalTag by hand for condDBv2: was %s, now %s' % (current_gt,new_gt)
    process.GlobalTag.globaltag = new_gt

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
        "/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/B291D129-7802-E511-8893-B8CA3A70A410.root"))
        #"/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/00D4DF0D-7735-E511-BD81-00259073E4F0.root"))
        


process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now

process.load("flashgg/MicroAOD/flashggMicroAODSequence_cff")

from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDefaultOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile.root'),
                               outputCommands = microAODDefaultOutputCommand,
                               )


# re-make vertices from charged packed candidates
process.load("flashgg/MicroAOD/NoMuonTrackProducer_cfi")

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.offlinePrimaryVerticesNoMu=offlinePrimaryVertices.clone()
process.offlinePrimaryVerticesNoMu.TrackLabel = cms.InputTag("NoMuonTrackProducer")

process.offlinePrimaryVerticesNoMu.TkFilterParameters = cms.PSet(
    algorithm=cms.string('filter'),
    maxNormalizedChi2 = cms.double(20000.0),
    minPixelLayersWithHits=cms.int32(1),
    minSiliconLayersWithHits = cms.int32(1),
    maxD0Significance = cms.double(5000.0), 
    minPt = cms.double(0.0),
    trackQuality = cms.string("any"))

process.offlinePrimaryVerticesNoMu.verbose = cms.untracked.bool(True)
process.noMuonVertexReco = cms.Sequence(process.NoMuonTrackProducer*process.offlinePrimaryVerticesNoMu)
#process.noMuonVertexReco = cms.Sequence(process.NoMuonTrackProducer)

# All jets are now handled in MicroAODCustomize.py
# Switch from PFCHS to PUPPI with puppi=1 argument (both if puppi=2)

process.p = cms.Path(process.noMuonVertexReco*process.flashggMicroAODSequence)
#process.p = cms.Path(process.flashggPDFWeightObject*process.flashggMicroAODSequence)
process.e = cms.EndPath(process.out)


from flashgg.MicroAOD.MicroAODCustomize import customize
customize(process)
