#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process = cms.Process("VHHadronicTagDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/160707_151052/0000/myMicroAODOutputFile_1.root"))


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("VHHadTagDump.root"),
                                   closeFileFast = cms.untracked.bool(True))


# met 
process.load("flashgg.MicroAOD.flashggMets_cfi")

# load tag sequence
process.load("flashgg.Taggers.flashggTagSequence_cfi")

#release cuts
process.flashggVHHadronicTag.jetPtThreshold  = cms.double(20.)
process.flashggVHHadronicTag.jetEtaThreshold = cms.double(5.0)
process.flashggVHHadronicTag.dijetMassLowThreshold = cms.double(0.)
process.flashggVHHadronicTag.dijetMassHighThreshold =cms.double(9999.)
process.flashggVHHadronicTag.cosThetaStarThreshold = cms.double(1.)

# use the preselected diphotons
#from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
#massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDiPhotons"),cms.InputTag("flashggPreselectedDiPhotons"))

# dumper
from flashgg.Taggers.tagsDumpers_cfi import createTagDumper
import flashgg.Taggers.dumperConfigTools as cfgTools

# set the VH had dumper
process.vhHadTagDumper = createTagDumper("VHHadronicTag")
process.vhHadTagDumper.dumpTrees     = True
process.vhHadTagDumper.dumpHistos    = True
process.vhHadTagDumper.dumpWorkspace = False


# get the variable list

diphoton_variables = ["mass            := diPhoton.mass",
                      "diphoton_mva    := diPhotonMVA.result",
                      "pho1_pt         := diPhoton.leadingPhoton.pt",
                      "pho1_eta        := diPhoton.leadingPhoton.eta",
                      "pho1_phi        := diPhoton.leadingPhoton.phi",
                      "pho1_energy     := diPhoton.leadingPhoton.energy",
                      "pho1_full5x5_r9 := diPhoton.leadingPhoton.full5x5_r9",
                      "pho1_idmva      := diPhoton.leadPhotonId",
                      "pho2_pt         := diPhoton.subLeadingPhoton.pt",
                      "pho2_eta        := diPhoton.subLeadingPhoton.eta",
                      "pho2_phi        := diPhoton.subLeadingPhoton.phi",
                      "pho2_energy     := diPhoton.subLeadingPhoton.energy",
                      "pho2_full5x5_r9 := diPhoton.subLeadingPhoton.full5x5_r9",
                      "pho2_idmva      := diPhoton.subLeadPhotonId"
                      ]

jet_variables      = ["jet1_pt      := leadingJet.pt()",
                      "jet1_eta     := leadingJet.eta()",
                      "jet1_phi     := leadingJet.phi()",
                      "jet1_energy  := leadingJet.energy()",
                      "jet2_pt      := subLeadingJet.pt()",
                      "jet2_eta     := subLeadingJet.eta()",
                      "jet2_phi     := subLeadingJet.phi()",
                      "jet2_energy  := subLeadingJet.energy()",
                     ]

all_variables = diphoton_variables + jet_variables

cfgTools.addCategories(process.vhHadTagDumper,
                       [
                           ("VHHadronicTag","leadingJet.pt>0",0)
                       ],
                       variables  = all_variables,
                       histograms = []
)

process.vhHadTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents" ,1000)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)

process.p1 = cms.Path(
    process.flashggMets+
    process.flashggTagSequence*
    process.vhHadTagDumper
)

print process.p1
