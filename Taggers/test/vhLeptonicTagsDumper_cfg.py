#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
import os

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet


process = cms.Process("VHHadronicTagDumper")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

# input files
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
#MC
"/store//group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/161114_134042/0000/myMicroAODOutputFile_1.root"
#data
#"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016E-23Sep2016-v1/161114_163114/0000/myMicroAODOutputFile_817.root"

))


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("vhHadDump.root"),
                                   closeFileFast = cms.untracked.bool(True))

# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

## import systs. customize
from flashgg.Systematics.SystematicsCustomize import *

# load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")

# apply scale and smearing corrections
useEGMTools(process)


## if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processId == 'Data' or customize.processId == 'data':
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


# load tag sequence
process.load("flashgg.Taggers.flashggTagSequence_cfi")
process.flashggTagSequence.remove(process.flashggUpdatedIdMVADiPhotons) # Needs to be run before systematics
massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggUpdatedIdMVADiPhotons"),cms.InputTag("flashggDiPhotonSystematics"))

#remove cut on jets
#process.flashggWHLeptonicTag.jetsNumberThreshold = cms.double(99999.)
#process.flashggVHLeptonicLooseTag.jetsNumberThreshold = cms.double(99999.)
process.flashggZHLeptonicTag.MVAThreshold = cms.double(-1.)
process.flashggWHLeptonicTag.MVAThreshold = cms.double(-1.)
process.flashggVHLeptonicLooseTag.MVAThreshold = cms.double(-1.)

# dumper
from flashgg.Taggers.tagsDumpers_cfi import createTagDumper
import flashgg.Taggers.dumperConfigTools as cfgTools

# set the VH dumpers
process.WHLeptonicTagDumper = createTagDumper("WHLeptonicTag")
process.WHLeptonicTagDumper.dumpTrees     = True
process.WHLeptonicTagDumper.dumpHistos    = True
process.WHLeptonicTagDumper.dumpWorkspace = False

process.VHLeptonicLooseTagDumper = createTagDumper("VHLeptonicLooseTag")
process.VHLeptonicLooseTagDumper.dumpTrees     = True
process.VHLeptonicLooseTagDumper.dumpHistos    = True
process.VHLeptonicLooseTagDumper.dumpWorkspace = False

process.ZHLeptonicTagDumper = createTagDumper("ZHLeptonicTag")
process.ZHLeptonicTagDumper.dumpTrees     = True
process.ZHLeptonicTagDumper.dumpHistos    = True
process.ZHLeptonicTagDumper.dumpWorkspace = False

# get the variable list
diphoton_variables = ["mass            := diPhoton.mass",
                      "diphoton_pt     := diPhoton.pt",
                      "diphoton_mva    := diPhotonMVA.result",
                      "pho1_pt         := diPhoton.leadingPhoton.pt",
                      "pho1_eta        := diPhoton.leadingPhoton.eta",
                      "pho1_phi        := diPhoton.leadingPhoton.phi",
                      "pho1_energy     := diPhoton.leadingPhoton.energy",
                      "pho1_full5x5_r9 := diPhoton.leadingPhoton.full5x5_r9",
                      "pho1_idmva      := diPhoton.leadPhotonId",
                      "pho1_genMatchType:=diPhoton.leadingPhoton.genMatchType",
                      "pho2_pt         := diPhoton.subLeadingPhoton.pt",
                      "pho2_eta        := diPhoton.subLeadingPhoton.eta",
                      "pho2_phi        := diPhoton.subLeadingPhoton.phi",
                      "pho2_energy     := diPhoton.subLeadingPhoton.energy",
                      "pho2_full5x5_r9 := diPhoton.subLeadingPhoton.full5x5_r9",
                      "pho2_idmva      := diPhoton.subLeadPhotonId",
                      "pho2_genMatchType:=diPhoton.subLeadingPhoton.genMatchType",
                      ]

leptons_variables = [ "mu1_pt         :=  ? muons.size()>0 ? muons[0].pt() : -100 ",
                      "mu1_phi        :=  ? muons.size()>0 ? muons[0].phi() : -100 ",
                      "mu1_eta        :=  ? muons.size()>0 ? muons[0].eta() : -100 ",
                      "mu1_energy     :=  ? muons.size()>0 ? muons[0].energy() : -100 ",
                      "mu2_pt         :=  ? muons.size()>1 ? muons[1].pt() : -100 ",
                      "mu2_ph1        :=  ? muons.size()>1 ? muons[1].phi() : -100 ",
                      "mu2_eta        :=  ? muons.size()>1 ? muons[1].eta() : -100 ",
                      "mu2_energy     :=  ? muons.size()>1 ? muons[1].energy() : -100 ",
                      "ele1_pt        :=  ? electrons.size()>0 ? electrons[0].pt() : -100 ",
                      "ele1_phi       :=  ? electrons.size()>0 ? electrons[0].phi() : -100 ",
                      "ele1_eta       :=  ? electrons.size()>0 ? electrons[0].eta() : -100 ",
                      "ele1_energy    :=  ? electrons.size()>0 ? electrons[0].energy() : -100 ",
                      "ele2_pt        :=  ? electrons.size()>1 ? electrons[1].pt() : -100 ",
                      "ele2_phi       :=  ? electrons.size()>1 ? electrons[1].phi() : -100 ",
                      "ele2_eta       :=  ? electrons.size()>1 ? electrons[1].eta() : -100 ",
                      "ele2_energy    :=  ? electrons.size()>1 ? electrons[1].energy() : -100 "
                     ]

jets_variables = ["njets := jets.size()",
                  "jet1_pt     :=  ? jets.size()>0 ? jets[0].pt : -100 ",
                  "jet1_phi    :=  ? jets.size()>0 ? jets[0].phi : -100 ",
                  "jet1_eta    :=  ? jets.size()>0 ? jets[0].eta : -100 ",
                  "jet1_energy :=  ? jets.size()>0 ? jets[0].energy : -100 ",
                  "jet2_pt     :=  ? jets.size()>1 ? jets[1].pt : -100 ",
                  "jet2_phi    :=  ? jets.size()>1 ? jets[1].phi : -100 ",
                  "jet2_eta    :=  ? jets.size()>1 ? jets[1].eta : -100 ",
                  "jet2_energy :=  ? jets.size()>1 ? jets[1].energy : -100 ",
                  "jet3_pt     :=  ? jets.size()>2 ? jets[2].pt : -100 ",
                  "jet3_phi    :=  ? jets.size()>2 ? jets[2].phi : -100 ",
                  "jet3_eta    :=  ? jets.size()>2 ? jets[2].eta : -100 ",
                  "jet3_energy :=  ? jets.size()>2 ? jets[2].energy : -100 "
                  ]

met_variables = ["met_pt  := met.corPt()",
                 "met_phi := met.corPhi()"]


gen_variables = ["hasZ := tagTruth().associatedZ",
                 "hasW := tagTruth().associatedW",
                 ]


all_variables = diphoton_variables + leptons_variables
if ( customize.processId != "Data" and customize.processId != "data"):
    all_variables+=gen_variables
    
cfgTools.addCategories(process.WHLeptonicTagDumper,
                       [
                           ("WHLeptonicTag","diPhoton.leadingPhoton.pt>0",0)
                       ],
                       variables  = all_variables+met_variables+jets_variables,
                       histograms = []
)

cfgTools.addCategories(process.VHLeptonicLooseTagDumper,
                       [
                           ("VHLeptonicLooseTag","diPhoton.leadingPhoton.pt>0",0)
                       ],
                       variables  = all_variables+met_variables+jets_variables,
                       histograms = []
)

cfgTools.addCategories(process.ZHLeptonicTagDumper,
                       [
                           ("ZHLeptonicTag","diPhoton.leadingPhoton.pt>0",0)
                       ],
                       variables  = all_variables,
                       histograms = []
)

process.WHLeptonicTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"
process.VHLeptonicLooseTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"
process.ZHLeptonicTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"





# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*"))
                                                               

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data" or customize.processId == "data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter
    


customize.setDefault("maxEvents" ,1000)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)


process.p1 = cms.Path(
        process.dataRequirements*
        process.flashggUpdatedIdMVADiPhotons*
        process.flashggDiPhotonSystematics*
        process.flashggTagSequence*
        process.WHLeptonicTagDumper*
        process.VHLeptonicLooseTagDumper*
        process.ZHLeptonicTagDumper
)

print process.p1

