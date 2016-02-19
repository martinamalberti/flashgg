import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("FlashggAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12' # keep updated for JEC

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
#"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160127_024939/0000/myMicroAODOutputFile_1.root"
"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall15DR76-1_3_0-25ns/1_3_0/DoubleEG/RunIIFall15DR76-1_3_0-25ns-1_3_0-v0-Run2015C_25ns-16Dec2015-v1/160116_105829/0000/myMicroAODOutputFile_1.root"
        ))

# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

# flashgg tag sequence (for dipho MVA) and jet collections
process.load("flashgg/Taggers/flashggTagSequence_cfi")
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag


# load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))

customize.processType = 'data'

# import flashgg customization to check if we have data or MC
# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    newvpset = cms.VPSet()
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if pset.Label.value().count("Scale"):
            pset.NoCentralShift = cms.bool(False) # Turn on central shift for data (it is off for MC)                                                                              
	    pset.NSigmas = cms.vint32() # Do only central value, no syst shifts
	    newvpset += [pset]
	if pset.Label.value().count("SigmaEOverESmear"):
            pset.NSigmas = cms.vint32() # Do only central value, no syst shifts
            newvpset += [pset]	
    process.flashggDiPhotonSystematics.SystMethods = newvpset
else:
    print 'mc'
    newvpset = cms.VPSet()
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if pset.Label.value().count("MCSmear"):
            pset.NSigmas = cms.vint32() # Do only central value, no syst. shifts
            newvpset += [pset]
	if pset.Label.value().count("SigmaEOverESmear"):
            pset.NSigmas = cms.vint32() # Do only central value, no syst shifts
            newvpset += [pset]	
    process.flashggDiPhotonSystematics.SystMethods = newvpset



from flashgg.Taggers.globalVariables_cff import globalVariables

process.analysisTree = cms.EDAnalyzer('EDSimpleTreeMaker',
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
				      isControlSample = cms.untracked.bool(False),
                                      bTag = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                      triggerBits = cms.InputTag('TriggerResults::HLT'),
                                      globalVariables = globalVariables
                                      )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mytree_unfiltered.root")
                                   )


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v1",
                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"),
                                         andOr    = cms.bool(True) # True = or between triggers 
                                         )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

process.dataRequirements = cms.Sequence()
if customize.processType == "data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter
        process.GlobalTag.globaltag = '76X_dataRun2_v15'

process.p = cms.Path(process.dataRequirements*
                     process.flashggDiPhotonSystematics*
                     process.flashggTagSequence
                     *process.analysisTree)



# set default options if needed                                                                                                                                                                        
customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",2.6e+3)
# call the customization
customize(process)

