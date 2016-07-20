import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v13')
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_v11')
else:
    raise Exception,"The default setup does not support releases other than 76X and 80X"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


## input file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        #data

        # mc
        "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160707_152047/0000/myMicroAODOutputFile_1.root"
        ))


# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

## import systs. customize
from flashgg.Systematics.SystematicsCustomize import *

# load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")


##customize.processId = 'Data'


# apply scale and smearing corrections
useEGMTools(process)

## if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processId == 'Data':
    print 'data' 
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    customizePhotonSystematicsForMC(process)
    ##syst (1D) 
    vpset   = process.flashggDiPhotonSystematics.SystMethods
    newvpset = cms.VPSet()
    for pset in vpset:
        pset.NSigmas = cms.vint32() # no up/down syst shifts
        pset.ApplyCentralValue = cms.bool(False) # no central value
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
        newvpset+= [pset]
    process.flashggDiPhotonSystematics.SystMethods = newvpset        
    ##syst (2D) : smearings with EGMTool
    vpset2D   = process.flashggDiPhotonSystematics.SystMethods2D
    newvpset2D = cms.VPSet()
    for pset in vpset2D:
        pset.NSigmas = cms.PSet( firstVar = cms.vint32(), secondVar = cms.vint32() ) # only central value, no up/down syst shifts (2D case)
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
            newvpset2D+= [pset]
    process.flashggDiPhotonSystematics.SystMethods2D = newvpset2D       

print 'syst 1D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])
print 'syst 2D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods2D])


# flashgg tag sequence (for dipho MVA) and jet collections
process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
process.flashggTagSequence.remove(process.flashggUpdatedIdMVADiPhotons) # Needs to be run before systematics
massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggUpdatedIdMVADiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))


## global variables to dump
from flashgg.Taggers.globalVariables_cff import globalVariables

## analyzer
process.analysisTree = cms.EDAnalyzer('EDSimpleTreeMaker',
                                      lumiWeight=cms.untracked.double(1000.),
                                      generatorInfo = cms.InputTag('generator'),  
                                      genParticleTag = cms.InputTag( "flashggPrunedGenParticles" ),
                                      PileUpTag = cms.InputTag('slimmedAddPileupInfo'),
                                      rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),
                                      VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                      MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                      inputTagJets= UnpackedJetCollectionVInputTag,
                                      GenJetTag=cms.InputTag( "slimmedGenJets"),
                                      ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                      MuonTag=cms.InputTag('flashggSelectedMuons'),
                                      METTag=cms.InputTag('slimmedMETs'),
                                      jetPtThreshold = cms.untracked.double(20.),
                                      electronPtThreshold = cms.untracked.double(10.),
                                      muonPtThreshold = cms.untracked.double(10.),
				      isControlSample = cms.untracked.bool(False),
                                      bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                      triggerBits = cms.InputTag('TriggerResults::HLT'),
                                      globalVariables = globalVariables
                                      )



process.analysisTree.globalVariables.addTriggerBits = cms.PSet(
    tag = cms.InputTag("TriggerResults::HLT"),
    bits = cms.vstring(
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v1",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v2",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v3",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v4",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v5"
    )
)
if customize.processId == "Data":
    process.analysisTree.triggerBits = cms.InputTag('TriggerResults::HLT') 
    process.analysisTree.globalVariables.addTriggerBits.tag = cms.InputTag("TriggerResults::HLT")
else:
    process.analysisTree.triggerBits = cms.InputTag('TriggerResults::HLT2') 
    process.analysisTree.globalVariables.addTriggerBits.tag = cms.InputTag("TriggerResults::HLT2")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mytree.root")
                                   )


## HLT filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ),
                                         andOr    = cms.bool(True) # True = or between triggers 
                                         )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## EE bas supercluster filter
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

## filters for data
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter



process.p = cms.Path(process.dataRequirements*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggTagSequence*
                     process.analysisTree)


## set default options if needed
customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)

