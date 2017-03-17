#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

## input files
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(      
        #data                    
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016E-23Sep2016-v1/161114_163114/0000/myMicroAODOutputFile_817.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/ReMiniAOD-03Feb2017-2_5_1/2_5_1/SingleElectron/ReMiniAOD-03Feb2017-2_5_1-2_5_1-v0-Run2016D-03Feb2017-v1/170214_121515/0000/myMicroAODOutputFile_550.root"
        # mc
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_0-25ns_Moriond17/2_4_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_30M/RunIISummer16-2_4_0-25ns_Moriond17-2_4_0-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/161225_213425/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/DYJetsToEE_M-50_LTbinned_95To100_5f_LO_13TeV-madgraph_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170118_090503/0000/myMicroAODOutputFile_5.root"
        "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-EGM0_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170121_195541/0000/myMicroAODOutputFile_9.root"
        ))

## output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)


## import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## load module to recompute photon id on-the-fly
process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")
print process.flashggUpdatedIdMVADiPhotons.src

## import flashgg customization to check if we have data, signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

## import systs. customize
from flashgg.Systematics.SystematicsCustomize import *
#customize.processId = 'Data' # for test

## load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
process.flashggDiPhotonSystematics.src = "flashggUpdatedIdMVADiPhotons"

## Or use the official  tool instead  ????????????????
useEGMTools(process)

## if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processId == 'Data':
    print 'data' 
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    customizePhotonSystematicsForMC(process)
#    ##syst (1D) 
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



## preselection Zee
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")
process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"
process.flashggPreselectedDiPhotons.src = "flashggDiPhotonSystematics"

## untagged tag dumper
process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.flashggTagSequence.remove(process.flashggUpdatedIdMVADiPhotons) # Needs to be run before systematics

process.flashggUntagged.Boundaries = cms.vdouble(-2)
process.flashggUntagged.DiPhotonTag    = "flashggPreselectedDiPhotons"
process.flashggDiPhotonMVA.DiPhotonTag = "flashggPreselectedDiPhotons"

#remove un-necessary tags
process.flashggTagSequence.remove(process.flashggVBFTag)
process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
process.flashggTagSequence.remove(process.flashggVHMetTag)
process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
process.flashggTagSequence.remove(process.flashggVHHadronicTag)
process.flashggTagSequence.remove(process.flashggTagSorter)


import flashgg.Taggers.dumperConfigTools as cfgTools
from flashgg.Taggers.tagsDumpers_cfi import createTagDumper
process.diphotonDumper = createTagDumper("UntaggedTag")
process.diphotonDumper.src = "flashggUntagged"
process.diphotonDumper.maxCandPerEvent = -1 # take them all
process.diphotonDumper.dumpTrees = True
process.diphotonDumper.dumpWorkspace = False
process.diphotonDumper.quietRooFit = True
#process.diphotonDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL"
process.diphotonDumper.nameTemplate ="tree_$SQRTS_$LABEL"

process.diphotonDumper.globalVariables.addTriggerBits = cms.PSet(
    tag = cms.InputTag("TriggerResults::HLT"),
    bits = cms.vstring(
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v1",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v2",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v3",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v4",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v5",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v6",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v7",
        "HLT_Ele27_WPTight_Gsf_v1",
        "HLT_Ele27_WPTight_Gsf_v2",
        "HLT_Ele27_WPTight_Gsf_v3",
        "HLT_Ele27_WPTight_Gsf_v4",
        "HLT_Ele27_WPTight_Gsf_v5",
        "HLT_Ele27_WPTight_Gsf_v6",
        "HLT_Ele27_WPTight_Gsf_v7",
        #"HLT_Ele35_WPLoose_Gsf_v1",
        #"HLT_Ele35_WPLoose_Gsf_v2",
        #"HLT_Ele35_WPLoose_Gsf_v3",
        #"HLT_Ele35_WPLoose_Gsf_v4",
        #"HLT_Ele35_WPLoose_Gsf_v5",
        #"HLT_Ele35_WPLoose_Gsf_v6",
        #"HLT_Ele45_WPLoose_Gsf_v1",
        #"HLT_Ele45_WPLoose_Gsf_v2",
        #"HLT_Ele45_WPLoose_Gsf_v3",
        #"HLT_Ele45_WPLoose_Gsf_v4",
        #"HLT_Ele45_WPLoose_Gsf_v5",
        #"HLT_Ele45_WPLoose_Gsf_v6",
    )
)


## define categories and associated objects to dump
cfgTools.addCategory(process.diphotonDumper,
                     "Reject",
                     "0",
                      -1 ## if nSubcat is -1 do not store anythings
                     )

## interestng categories 
cfgTools.addCategories(process.diphotonDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [
                        ("All","1",0)
                       ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=["mass               :=diPhoton.mass",
                                  "pt                 :=diPhoton.pt",
                                  "diphotonMVA        :=diPhotonMVA.result",

                                  "ele1_et            :=diPhoton.leadingPhoton.pt",
                                  "ele1_energy        :=diPhoton.leadingPhoton.energy",
                                  "ele1_rawEnergy     :=diPhoton.leadingPhoton.superCluster.rawEnergy",
                                  "ele1_esEnergy      :=diPhoton.leadingPhoton.superCluster.preshowerEnergy",
                                  "ele1_eTrue         := ?diPhoton.leadingPhoton().hasMatchedGenPhoton()?diPhoton.leadingPhoton().matchedGenPhoton().energy():0",
                                  "ele1_scEta         :=diPhoton.leadingPhoton.superCluster.eta",
                                  "ele1_scPhi         :=diPhoton.leadingPhoton.superCluster.phi",
                                  "ele1_eta           :=diPhoton.leadingPhoton.eta",
                                  "ele1_phi           :=diPhoton.leadingPhoton.phi",
                                  "ele1_r9            :=diPhoton.leadingPhoton.old_r9()",
                                  "ele1_full5x5_r9    :=diPhoton.leadingPhoton.full5x5_r9",
                                  "ele1_etaWidth      :=diPhoton.leadingPhoton.superCluster.etaWidth",
                                  "ele1_s4            :=diPhoton.leadingPhoton.s4",
                                  "ele1_full5x5_sigmaIetaIeta := diPhoton.leadingPhoton.full5x5_sigmaIetaIeta",
                                  "ele1_sigmaEoE      :=diPhoton.leadingPhoton.sigEOverE",
                                  "ele1_pfPhoIso03    :=diPhoton.leadingPhoton.pfPhoIso03",
                                  "ele1_chIsoRV       :=diPhoton.leadingView.pfChIso03WrtChosenVtx()",
                                  "ele1_chIsoWV       :=diPhoton.leadingPhoton.pfChgIsoWrtWorstVtx03",
                                  "ele1_ESoRawE       :=diPhoton.leadingPhoton.superCluster.preshowerEnergy()/diPhoton.leadingPhoton.superCluster.rawEnergy()",
                                  "ele1_idmva         :=diPhoton.leadPhotonId",
                                  "ele1_kHasSwitchToGain6 := diPhoton().leadingPhoton.checkStatusFlag('kHasSwitchToGain6')",
                                  "ele1_kHasSwitchToGain1 := diPhoton().leadingPhoton.checkStatusFlag('kHasSwitchToGain1')",

                                  # additional vars for checks in regression 
                                  "ele1_unsmearedSigmaEoE       :=diPhoton.leadingPhoton.userFloat('unsmearedSigmaEoE')",
                                  "ele1_reco_E                  :=diPhoton.leadingPhoton.userFloat('reco_E')",
                                  "ele1_reco_regr_E             :=diPhoton.leadingPhoton.userFloat('reco_regr_E')",
                                  "ele1_beforeShShTransf_regr_E :=diPhoton.leadingPhoton.userFloat('beforeShShTransf_regr_E')",
                                  "ele1_afterShShTransf_regr_E  :=diPhoton.leadingPhoton.userFloat('afterShShTransf_regr_E')",

                                  "ele2_et            :=diPhoton.subLeadingPhoton.pt",
                                  "ele2_energy        :=diPhoton.subLeadingPhoton.energy",
                                  "ele2_rawEnergy     :=diPhoton.subLeadingPhoton.superCluster.rawEnergy",
                                  "ele2_esEnergy      :=diPhoton.subLeadingPhoton.superCluster.preshowerEnergy",
                                  "ele2_eTrue         := ?diPhoton.subLeadingPhoton().hasMatchedGenPhoton()?diPhoton.subLeadingPhoton().matchedGenPhoton().energy():0",
                                  "ele2_scEta         :=diPhoton.subLeadingPhoton.superCluster.eta",
                                  "ele2_scPhi         :=diPhoton.subLeadingPhoton.superCluster.phi",
                                  "ele2_eta           :=diPhoton.subLeadingPhoton.eta",
                                  "ele2_phi           :=diPhoton.subLeadingPhoton.phi",
                                  "ele2_r9            :=diPhoton.subLeadingPhoton.old_r9()",
                                  "ele2_full5x5_r9    :=diPhoton.subLeadingPhoton.full5x5_r9",
                                  "ele2_etaWidth      :=diPhoton.subLeadingPhoton.superCluster.etaWidth",
                                  "ele2_s4            :=diPhoton.subLeadingPhoton.s4",
                                  "ele2_full5x5_sigmaIetaIeta := diPhoton.subLeadingPhoton.full5x5_sigmaIetaIeta",
                                  "ele2_sigmaEoE      :=diPhoton.subLeadingPhoton.sigEOverE",
                                  "ele2_pfPhoIso03    :=diPhoton.subLeadingPhoton.pfPhoIso03",
                                  "ele2_chIsoRV       :=diPhoton.subLeadingView.pfChIso03WrtChosenVtx()",
                                  "ele2_chIsoWV       :=diPhoton.subLeadingPhoton.pfChgIsoWrtWorstVtx03",
                                  "ele2_ESoRawE       :=diPhoton.subLeadingPhoton.superCluster.preshowerEnergy()/diPhoton.subLeadingPhoton.superCluster.rawEnergy()",
                                  "ele2_idmva         :=diPhoton.subLeadPhotonId",
                                  "ele2_kHasSwitchToGain6 := diPhoton().subLeadingPhoton.checkStatusFlag('kHasSwitchToGain6')",
                                  "ele2_kHasSwitchToGain1 := diPhoton().subLeadingPhoton.checkStatusFlag('kHasSwitchToGain1')",


                                  # additional vars for checks in regression 
                                  "ele2_unsmearedSigmaEoE       :=diPhoton.subLeadingPhoton.userFloat('unsmearedSigmaEoE')",
                                  "ele2_reco_E                  :=diPhoton.subLeadingPhoton.userFloat('reco_E')",
                                  "ele2_reco_regr_E             :=diPhoton.subLeadingPhoton.userFloat('reco_regr_E')",
                                  "ele2_beforeShShTransf_regr_E :=diPhoton.subLeadingPhoton.userFloat('beforeShShTransf_regr_E')",
                                  "ele2_afterShShTransf_regr_E  :=diPhoton.subLeadingPhoton.userFloat('afterShShTransf_regr_E')",
                                  ],
                       histograms=[]
                       )


## HLT + EE bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only) 
process.dataRequirements = cms.Sequence()
if customize.processId == 'Data':
    process.dataRequirements += process.eeBadScFilter
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
            #DoubleEG
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*",
            #SingleEG
            "HLT_Ele27_WPTight_Gsf_v*",
            #"HLT_Ele35_WPLoose_Gsf_v*",
            #"HLT_Ele45_WPLoose_Gsf_v*",
            ) )
else:
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
    #DoubleEG
    #"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*",
    #SingleEG
    #"HLT_Ele23_WPLoose_Gsf_v*" # 7_6_X
    ##"HLT_Ele27_WPLoose_Gsf_v*",
    ##"HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
    ##"HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
    ) )

 

process.p = cms.Path(process.hltHighLevel*
                     process.dataRequirements*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggTagSequence*
                     process.diphotonDumper)
    
#printSystematicInfo(process)

## set default options if needed
customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1e+3)
## call the customization
customize(process)
