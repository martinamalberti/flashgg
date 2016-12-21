#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
import os

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
#"/store//group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/161114_134042/0000/myMicroAODOutputFile_1.root"
#data
"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016E-23Sep2016-v1/161114_163114/0000/myMicroAODOutputFile_817.root"

))


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("vhHadDump.root"),
                                   closeFileFast = cms.untracked.bool(True))


# met 
#process.load("flashgg.MicroAOD.flashggMets_cfi")

# load tag sequence
process.load("flashgg.Taggers.flashggTagSequence_cfi")

# ivert photon ID on one photon to get a control sample
process.flashggPreselectedDiPhotons.cut = cms.string(
        "    (leadingPhoton.full5x5_r9>0.8||leadingPhoton.egChargedHadronIso<20||leadingPhoton.egChargedHadronIso/leadingPhoton.pt<0.3)"
        " && (subLeadingPhoton.full5x5_r9>0.8||subLeadingPhoton.egChargedHadronIso<20||subLeadingPhoton.egChargedHadronIso/subLeadingPhoton.pt<0.3)"
        " && (leadingPhoton.hadronicOverEm < 0.08 && subLeadingPhoton.hadronicOverEm < 0.08)"
        " && (leadingPhoton.pt >30.0 && subLeadingPhoton.pt > 20.0)"
        " && (abs(leadingPhoton.superCluster.eta) < 2.5 && abs(subLeadingPhoton.superCluster.eta) < 2.5)"
        " && (abs(leadingPhoton.superCluster.eta) < 1.4442 || abs(leadingPhoton.superCluster.eta) > 1.566)"
        " && (abs(subLeadingPhoton.superCluster.eta) < 1.4442 || abs(subLeadingPhoton.superCluster.eta) > 1.566)"
        " && ( (leadPhotonId > -0.9 && subLeadPhotonId < -0.9) || (leadPhotonId < -0.9 && subLeadPhotonId > -0.9) )" # invert ID on one of the photons
        )

#release cuts
process.flashggVHHadronicTag.jetPtThreshold  = cms.double(20.)
#process.flashggVHHadronicTag.jetEtaThreshold = cms.double(5.0)
process.flashggVHHadronicTag.leadPhoOverMassThreshold = cms.double(0.333)
process.flashggVHHadronicTag.subleadPhoOverMassThreshold = cms.double(0.25)
process.flashggVHHadronicTag.dijetMassLowThreshold = cms.double(0.)
process.flashggVHHadronicTag.dijetMassHighThreshold =cms.double(9999.)
process.flashggVHHadronicTag.cosThetaStarThreshold = cms.double(1.)
process.flashggVHHadronicTag.dRJetToPhoLThreshold = cms.double(0.4)
process.flashggVHHadronicTag.dRJetToPhoSThreshold = cms.double(0.4)
process.flashggVHHadronicTag.phoIdMVAThreshold = cms.double(-1.0)

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

jet_variables      = ["jet1_pt      := leadingJet.pt()",
                      "jet1_eta     := leadingJet.eta()",
                      "jet1_phi     := leadingJet.phi()",
                      "jet1_energy  := leadingJet.energy()",
                      "jet1_qgl     := leadingJet.QGL()",
                      "jet1_bdisc   := leadingJet.bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
                      "jet2_pt      := subLeadingJet.pt()",
                      "jet2_eta     := subLeadingJet.eta()",
                      "jet2_phi     := subLeadingJet.phi()",
                      "jet2_energy  := subLeadingJet.energy()",
                      "jet2_qgl     := subLeadingJet.QGL()",
                      "jet2_bdisc   := subLeadingJet.bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"
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

# Require standard diphoton trigger
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.hltHighLevel
    process.dataRequirements += process.eeBadScFilter
    


customize.setDefault("maxEvents" ,1000)    # max-number of events
customize.setDefault("targetLumi",1e+3) # define integrated lumi
customize(process)


process.p1 = cms.Path(
#    process.flashggMets+
   process.dataRequirements*
   process.flashggTagSequence*
   process.vhHadTagDumper
)

print process.p1
