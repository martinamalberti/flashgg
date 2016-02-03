import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' # keep updated for JEC

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
	#filtered
	#"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-TTHFiltered-1_2_0-25ns/1_2_0/ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-TTHFiltered-1_2_0-25ns-1_2_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160113_173025/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReReco74X-1_1_0-25ns/1_1_0/DoubleEG/RunIISpring15-ReReco74X-1_1_0-25ns-1_1_0-v0-Run2015D-04Dec2015-v2/160112_095813/0000/myMicroAODOutputFile_1.root"
        "/store/group/phys_higgs/cmshgg/sani/flashgg/GJetsHTReMiniAOD_v2/129-g1b8ab40/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GJetsHTReMiniAOD_v2-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1-a5c635afc9005d9960f109ef79f2aebc/160126_101144/0000/myMicroAODOutputFile_568.root",
"/store/group/phys_higgs/cmshgg/sani/flashgg/GJetsHTReMiniAOD_v2/129-g1b8ab40/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GJetsHTReMiniAOD_v2-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1-a5c635afc9005d9960f109ef79f2aebc/160126_101144/0000/myMicroAODOutputFile_569.root",
"/store/group/phys_higgs/cmshgg/sani/flashgg/GJetsHTReMiniAOD_v2/129-g1b8ab40/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GJetsHTReMiniAOD_v2-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1-a5c635afc9005d9960f109ef79f2aebc/160126_101144/0000/myMicroAODOutputFile_570.root",
"/store/group/phys_higgs/cmshgg/sani/flashgg/GJetsHTReMiniAOD_v2/129-g1b8ab40/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GJetsHTReMiniAOD_v2-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1-a5c635afc9005d9960f109ef79f2aebc/160126_101144/0000/myMicroAODOutputFile_571.root",
        ))

# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

# flashgg tag sequence (for dipho MVA) and jet collections
process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

from flashgg.Taggers.globalVariables_cff import globalVariables

process.analysisTree = cms.EDAnalyzer('EDSimpleTreeMaker',
                                      lumiWeight=cms.untracked.double(1.),
                                      generatorInfo = cms.InputTag('generator'),  
                                      genParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),
                                      PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
                                      rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      #DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                      DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                      MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                      inputTagJets= UnpackedJetCollectionVInputTag,
                                      GenJetTag=cms.InputTag( "slimmedGenJets"),
                                      ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                      MuonTag=cms.InputTag('flashggSelectedMuons'),
                                      METTag=cms.InputTag('slimmedMETs'),
                                      jetPtThreshold = cms.untracked.double(20.),
				      isControlSample = cms.untracked.bool(False),
                                      bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                      triggerBits = cms.InputTag('TriggerResults::HLT'),
                                      globalVariables = globalVariables
                                      )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mytree_unfiltered.root")
                                   )


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v1") )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

process.dataRequirements = cms.Sequence()
if customize.processType == "data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter


process.p = cms.Path(process.dataRequirements*
                     process.flashggTagSequence
                     *process.analysisTree)



# set default options if needed                                                                                                                                                                        
#from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.6e+3)
# call the customization
customize(process)

