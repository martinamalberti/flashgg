import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = ["HLT_Photon50_v*"]

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#process.GlobalTag.globaltag = 'auto:run2_mc'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

# Fix because auto:run2_mc points to MCRUN2_74_V9::All
current_gt = process.GlobalTag.globaltag.value()
if current_gt.count("::All"):
    new_gt = current_gt.replace("::All","")
    print 'Removing "::All" from GlobalTag by hand for condDBv2: was %s, now %s' % (current_gt,new_gt)
    process.GlobalTag.globaltag = new_gt

# Spring15
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
        "/store/data/Run2015D/SinglePhoton/MINIAOD/05Oct2015-v1/10000/C8770AB6-B76F-E511-8D38-0026189437F0.root"
        #"/store/mc/RunIISpring15DR74/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v3/00000/00BB8CCC-9208-E511-8A90-0025905B857C.root"
        ))
process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now

process.load("flashgg/MicroAOD/flashggMicroAODPhotonJetValidationSequence_cff")

from flashgg.MicroAOD.flashggMicroAODPhotonJetValidationOutputCommands_cff import microAODPhotonJetValidationOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myMicroAODOutputFile.root'),
                               outputCommands = microAODPhotonJetValidationOutputCommand
                               )

process.p = cms.Path(process.hltHighLevel+process.flashggMicroAODPhotonJetValidationSequence)
process.e = cms.EndPath(process.out)

#from flashgg.MicroAOD.MicroAODCustomize import customize
#customize(process)
