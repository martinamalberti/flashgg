import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAODAndTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1000 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Spring14miniaod/GluGluToHToGG_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/24926621-F11C-E411-AB9A-02163E008D0B.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Spring14miniaod/TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/049C0F9C-E61E-E411-9388-D8D385AE8466.root"))                                                                                                                            
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_153138/0000/myMicroAODOutputFile_2.root"
        "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_152108/0000/myMicroAODOutputFile_1.root"
        ))

process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

process.treeMaker = cms.EDAnalyzer('FlashggtthOptimizationTreeMaker',
                                   lumiWeight=cms.untracked.double(1000.),
                                   generatorInfo = cms.InputTag('generator'),  
                                   PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                   MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                   inputTagJets= UnpackedJetCollectionVInputTag,
                                   GenJetTag=cms.InputTag( "slimmedGenJets"),
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   jetPtThreshold = cms.untracked.double(20.),
                                   bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                  )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mytree.root")
                                   )

process.p = cms.Path(process.flashggTagSequence*process.treeMaker)
