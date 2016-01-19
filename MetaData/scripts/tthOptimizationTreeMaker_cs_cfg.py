import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'POSTLS170_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Spring14miniaod/GluGluToHToGG_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/24926621-F11C-E411-AB9A-02163E008D0B.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Spring14miniaod/TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/049C0F9C-E61E-E411-9388-D8D385AE8466.root"))                                                                                                                            
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_153138/0000/myMicroAODOutputFile_2.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_152108/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/DoubleEG/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-Run2015D-05Oct2015-v1/151021_151712/0000/myMicroAODOutputFile_1.root"
        
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-1_1_0-25ns/1_1_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-1_1_0-25ns-1_1_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160105_224456/0000/myMicroAODOutputFile_2.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-Prompt-1_1_0-25ns/1_1_0/DoubleEG/RunIISpring15-Prompt-1_1_0-25ns-1_1_0-v0-Run2015D-PromptReco-v4/160105_225454/0000/myMicroAODOutputFile_824.root"
	 "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-Prompt-1_1_0-25ns/1_1_0/DoubleEG/RunIISpring15-Prompt-1_1_0-25ns-1_1_0-v0-Run2015D-PromptReco-v4/160105_225454/0000/myMicroAODOutputFile_824.root"
		))


# flashgg tag sequence (for dipho MVA)
process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

from flashgg.Taggers.globalVariables_cff import globalVariables

process.analysisTree = cms.EDAnalyzer('FlashggtthOptimizationTreeMaker',
                                      lumiWeight=cms.untracked.double(1000.),
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
				      isControlSample =cms.untracked.bool(True),
                                      bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                      triggerBits = cms.InputTag('TriggerResults::HLT'),
                                      globalVariables = globalVariables
                                      )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mytree.root")
                                   )

# import flashgg customization to check if we have data or MC
from flashgg.MetaData.JobConfig import customize
customize.parse()

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v1") )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

customize.processType = "data"

process.dataRequirements = cms.Sequence()
if customize.processType == "data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter


from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
process.flashggPreselectedDiPhotons.cut = cms.string(
		"(leadingPhoton.pt >30.0 && subLeadingPhoton.pt > 20.0)"
		" && (abs(leadingPhoton.superCluster.eta) < 2.5 && abs(subLeadingPhoton.superCluster.eta) < 2.5)"
		" && (abs(leadingPhoton.superCluster.eta) < 1.4442 || abs(leadingPhoton.superCluster.eta) > 1.566)"
		" && (abs(subLeadingPhoton.superCluster.eta) < 1.4442 || abs(subLeadingPhoton.superCluster.eta) > 1.566)" 
		" && (((leadingPhoton.full5x5_r9>0.8||leadingPhoton.egChargedHadronIso<20||leadingPhoton.egChargedHadronIso/leadingPhoton.pt<0.3)"
		" && ( leadingPhoton.hadronicOverEm < 0.08 ) && (leadPhotonId > -0.9 && subLeadPhotonId < -0.9))"
		" || ((subLeadingPhoton.full5x5_r9>0.8||subLeadingPhoton.egChargedHadronIso<20||subLeadingPhoton.egChargedHadronIso/subLeadingPhoton.pt<0.3)"
		" && (subLeadingPhoton.hadronicOverEm < 0.08) &&  (leadPhotonId < -0.9 && subLeadPhotonId > -0.9)))"
        )

process.p = cms.Path(process.dataRequirements*process.flashggTagSequence*process.analysisTree)

# set default options if needed                                                                                                                                                                        
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.63e+3)
# call the customization
customize(process)

