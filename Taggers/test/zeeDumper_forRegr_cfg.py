#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag,cloneProcessingSnippet

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_v11')
else:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


## apply shower shape corrections
#doUpdatedIdMVADiPhotons = False # set to True for 76X (for 80X shower shape corrections not yet available)                                                                                   

## input files
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        #data                                                                                                                                                                  
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/SingleElectron/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-R\un2016C-PromptReco-v2/160707_145632/0000/myMicroAODOutputFile_1.root"                                                                                                      
        "file:/afs/cern.ch/work/m/malberti/HGG/devel/CMSSW_8_0_8/src/flashgg/myMicroAODOutputFile.root"
        # mc                                                                                                                                                               
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_\MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/160707_143004/0000/myMicroAODOutputFile_1.root"
        ))


## output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)



## load module to recompute photon id on-the-fly
process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")

## preselection Zee
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")
process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"

## load tag sequence, keeping only untagged tags
process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.flashggTagSequence.remove(process.flashggVBFMVA)
process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)
process.flashggTagSequence.remove(process.flashggVBFTag)
process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
process.flashggUntagged.Boundaries = cms.vdouble(-2)
process.flashggTagSorter.TagPriorityRanges = cms.VPSet(cms.PSet(TagName = cms.InputTag('flashggUntagged')))
process.flashggTagSorter.MassCutUpper=cms.double(9999.)
process.flashggTagSorter.MassCutLower=cms.double(60.)
print process.flashggTagSequence


## import flashgg customization
from flashgg.MetaData.JobConfig import customize
customize.parse()


## import systs. customization
from flashgg.Systematics.SystematicsCustomize import *

systlabels = [""]
phosystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []
jetSystematicsInputTags = []
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

## load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
print "input to flashggDiPhotonSystematics = ", process.flashggDiPhotonSystematics.src

# EGM smearer
useEGMTools(process)

# if data, apply only energy scale corrections
if customize.processId == 'Data':
    print 'Data' 
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    customizePhotonSystematicsForMC(process)
    #syst(1D)
    vpset   = process.flashggDiPhotonSystematics.SystMethods
    newvpset = cms.VPSet()
    process.flashggDiPhotonSystematics.SystMethods = newvpset
    #syst (2D) : smearings with EGMTool
    vpset2D   = process.flashggDiPhotonSystematics.SystMethods2D
    newvpset2D = cms.VPSet()
    process.flashggDiPhotonSystematics.SystMethods2D = newvpset2D       

print 'syst 1D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])
print 'syst 2D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods2D])


print 'Cloning Tag seq for each syst...'
cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,jetsystlabels,jetSystematicsInputTags)
#print process.systematicsTagSequences 


#Tag dumper
process.load("flashgg.Taggers.diphotonTagDumper_cfi") ##  import diphotonTagDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools
process.tagsDumper.className = "DiPhotonTagDumper"
process.tagsDumper.src = "flashggSystTagMerger"
process.tagsDumper.dumpTrees = True
process.tagsDumper.dumpWorkspace = False
process.tagsDumper.quietRooFit = True
process.tagsDumper.processId = "test"
process.tagsDumper.nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL")

process.tagsDumper.globalVariables.addTriggerBits = cms.PSet(
    tag = cms.InputTag("TriggerResults::HLT"),
    bits = cms.vstring(
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v1",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v2",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v3",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v4",
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v5",
        "HLT_Ele27_WPTight_Gsf_v1",
        "HLT_Ele27_WPTight_Gsf_v2",
        "HLT_Ele27_WPTight_Gsf_v3",
        "HLT_Ele27_WPTight_Gsf_v4",
        "HLT_Ele27_WPTight_Gsf_v5",
        "HLT_Ele35_WPLoose_Gsf_v1",
        "HLT_Ele35_WPLoose_Gsf_v2",
        "HLT_Ele35_WPLoose_Gsf_v3",
        "HLT_Ele35_WPLoose_Gsf_v4",
        "HLT_Ele35_WPLoose_Gsf_v5"
        #"HLT_Ele27_WPLoose_Gsf_v",
        #"HLT_Ele27_eta2p1_WPLoose_Gsf_v",
        #"HLT_Ele22_eta2p1_WPLoose_Gsf_v"
    )
)

 
## define categories and associated objects to dump
tagList=[
["UntaggedTag",0]
]


definedSysts=set()
process.tagsDumper.classifierCfg.remap=cms.untracked.VPSet()
for tag in tagList: 
    tagName=tag[0]
    tagCats=tag[1]
  # remap return value of class-based classifier
    process.tagsDumper.classifierCfg.remap.append( cms.untracked.PSet( src=cms.untracked.string("flashgg%s"%tagName), dst=cms.untracked.string(tagName) ) )
    for systlabel in systlabels:
        print systlabel
        if not systlabel in definedSysts:
            # the cut corresponding to the systematics can be defined just once
            cutstring = "hasSyst(\"%s\") "%(systlabel)
            print systlabel, cutstring
            definedSysts.add(systlabel)
        else:
            cutstring = None

        # interestng categories 
        cfgTools.addCategory(process.tagsDumper,
                             systlabel, #syst. label
                             classname=tagName,#catname
                             cutbased=cutstring, #cutstring
                             subcats=0,#n subcats                         
                             variables=["mass               :=diPhoton.mass",
                                        "pt                 :=diPhoton.pt",
                                        "diphotonMVA        :=diPhotonMVA.result",
                                        
                                        "pho1_et            :=diPhoton.leadingPhoton.pt",
                                        "pho1_energy        :=diPhoton.leadingPhoton.energy",
                                        "pho1_eTrue         := ?diPhoton.leadingPhoton().hasMatchedGenPhoton()?diPhoton.leadingPhoton().matchedGenPhoton().energy():0",
                                        "pho1_eta           :=diPhoton.leadingPhoton.eta",
                                        "pho1_phi           :=diPhoton.leadingPhoton.phi",
                                        "pho1_idmva         :=diPhoton.leadPhotonId",

                                        "pho2_et            :=diPhoton.subLeadingPhoton.pt",
                                        "pho2_energy        :=diPhoton.subLeadingPhoton.energy",
                                        "pho2_eTrue         := ?diPhoton.subLeadingPhoton().hasMatchedGenPhoton()?diPhoton.subLeadingPhoton().matchedGenPhoton().energy():0",
                                        "pho2_eta           :=diPhoton.subLeadingPhoton.eta",
                                        "pho2_phi           :=diPhoton.subLeadingPhoton.phi",
                                        "pho2_idmva         :=diPhoton.subLeadPhotonId",
                                        
                                        # variables for energy regression
                                        "pho1_isEB             := diPhoton.leadingPhoton.isEB",
                                        "pho1_numberOfClusters := diPhoton.leadingPhoton.superCluster.size()",
                                        #"pho1_missingClusters  := !(diPhoton.leadingPhoton.superCluster.clusters().at(diPhoton.leadingPhoton.superCluster.size()-1).isAvailable())",
                                        #"pho1_missingClusters  := !diPhoton.leadingPhoton.superCluster.clusters().back().isAvailable()",
                                        "pho1_rawEnergy        := diPhoton.leadingPhoton.superCluster.rawEnergy",
                                        "pho1_scEnergy         := diPhoton.leadingPhoton.superCluster.energy",
                                        "pho1_full5x5_r9       := diPhoton.leadingPhoton.full5x5_r9",
                                        "pho1_r9               := diPhoton.leadingPhoton.r9",
                                        "pho1_etaWidth         := diPhoton.leadingPhoton.superCluster.etaWidth",
                                        "pho1_phiWidth         := diPhoton.leadingPhoton.superCluster.phiWidth",
                                        "pho1_hadronicOverEm   := diPhoton.leadingPhoton.hadronicOverEm",
                                        "pho1_scEta            := diPhoton.leadingPhoton.superCluster.eta",
                                        "pho1_scPhi            := diPhoton.leadingPhoton.superCluster.phi",
                                        "pho1_seedEta          := diPhoton.leadingPhoton.superCluster.seed.eta()",
                                        "pho1_seedPhi          := diPhoton.leadingPhoton.superCluster.seed.phi()",
                                        "pho1_seedEnergy       := diPhoton.leadingPhoton.seedEnergy",
                                        "pho1_e3x3             := diPhoton.leadingPhoton.e3x3",
                                        "pho1_e5x5             := diPhoton.leadingPhoton.e5x5",
                                        "pho1_sigmaIetaIeta    := diPhoton.leadingPhoton.showerShapeVariables.sigmaIetaIeta",
                                        "pho1_sigmaIphiIphi    := diPhoton.leadingPhoton.showerShapeVariables.sigmaIphiIphi",
                                        "pho1_sigmaIetaIphi    := diPhoton.leadingPhoton.showerShapeVariables.sigmaIetaIphi",
                                        "pho1_maxEnergyXtal    := diPhoton.leadingPhoton.showerShapeVariables.maxEnergyXtal",
                                        "pho1_e2nd             := diPhoton.leadingPhoton.showerShapeVariables.e2nd",
                                        "pho1_eTop             := diPhoton.leadingPhoton.showerShapeVariables.eTop",
                                        "pho1_eBottom          := diPhoton.leadingPhoton.showerShapeVariables.eBottom",
                                        "pho1_eLeft            := diPhoton.leadingPhoton.showerShapeVariables.eLeft",
                                        "pho1_eRight           := diPhoton.leadingPhoton.showerShapeVariables.eRight",
                                        "pho1_e2x5Max          := diPhoton.leadingPhoton.showerShapeVariables.e2x5Max",
                                        "pho1_e2x5Left         := diPhoton.leadingPhoton.showerShapeVariables.e2x5Left",
                                        "pho1_e2x5Right        := diPhoton.leadingPhoton.showerShapeVariables.e2x5Right",
                                        "pho1_e2x5Top          := diPhoton.leadingPhoton.showerShapeVariables.e2x5Top",
                                        "pho1_e2x5Bottom       := diPhoton.leadingPhoton.showerShapeVariables.e2x5Bottom",
                                        "pho1_preshowerEnergy  := diPhoton.leadingPhoton.superCluster.preshowerEnergy",
                                        "pho1_preshowerEnergyPlane1  := diPhoton.leadingPhoton.superCluster.preshowerEnergyPlane1",
                                        "pho1_preshowerEnergyPlane2  := diPhoton.leadingPhoton.superCluster.preshowerEnergyPlane2",
                                        "pho1_cryEta           := diPhoton.leadingPhoton.cryEta",
                                        "pho1_cryPhi           := diPhoton.leadingPhoton.cryPhi",
                                        "pho1_ieta             := diPhoton.leadingPhoton.iEta",
                                        "pho1_iphi             := diPhoton.leadingPhoton.iPhi",


                                        "pho2_isEB             := diPhoton.subLeadingPhoton.isEB",
                                        "pho2_numberOfClusters := diPhoton.subLeadingPhoton.superCluster.size()",
                                        #"pho2_missingClusters  := !diPhoton.subLeadingPhoton.superCluster.clusters()[diPhoton.subLeadingPhoton.superCluster.size()-1].isAvailable()",
                                        "pho2_rawEnergy        := diPhoton.subLeadingPhoton.superCluster.rawEnergy",
                                        "pho2_scEnergy         := diPhoton.subLeadingPhoton.superCluster.energy",
                                        "pho2_full5x5_r9       := diPhoton.subLeadingPhoton.full5x5_r9",
                                        "pho2_r9               := diPhoton.subLeadingPhoton.r9",
                                        "pho2_etaWidth         := diPhoton.subLeadingPhoton.superCluster.etaWidth",
                                        "pho2_phiWidth         := diPhoton.subLeadingPhoton.superCluster.phiWidth",
                                        "pho2_hadronicOverEm   := diPhoton.subLeadingPhoton.hadronicOverEm",
                                        "pho2_scEta            := diPhoton.subLeadingPhoton.superCluster.eta",
                                        "pho2_scPhi            := diPhoton.subLeadingPhoton.superCluster.phi",
                                        "pho2_seedEta          := diPhoton.subLeadingPhoton.superCluster.seed.eta()",
                                        "pho2_seedPhi          := diPhoton.subLeadingPhoton.superCluster.seed.phi()",
                                        "pho2_seedEnergy       := diPhoton.subLeadingPhoton.seedEnergy",
                                        "pho2_e3x3             := diPhoton.subLeadingPhoton.e3x3",
                                        "pho2_e5x5             := diPhoton.subLeadingPhoton.e5x5",
                                        "pho2_sigmaIetaIeta    := diPhoton.subLeadingPhoton.showerShapeVariables.sigmaIetaIeta",
                                        "pho2_sigmaIphiIphi    := diPhoton.subLeadingPhoton.showerShapeVariables.sigmaIphiIphi",
                                        "pho2_sigmaIetaIphi    := diPhoton.subLeadingPhoton.showerShapeVariables.sigmaIetaIphi",
                                        "pho2_maxEnergyXtal    := diPhoton.subLeadingPhoton.showerShapeVariables.maxEnergyXtal",
                                        "pho2_e2nd             := diPhoton.subLeadingPhoton.showerShapeVariables.e2nd",
                                        "pho2_eTop             := diPhoton.subLeadingPhoton.showerShapeVariables.eTop",
                                        "pho2_eBottom          := diPhoton.subLeadingPhoton.showerShapeVariables.eBottom",
                                        "pho2_eLeft            := diPhoton.subLeadingPhoton.showerShapeVariables.eLeft",
                                        "pho2_eRight           := diPhoton.subLeadingPhoton.showerShapeVariables.eRight",
                                        "pho2_e2x5Max          := diPhoton.subLeadingPhoton.showerShapeVariables.e2x5Max",
                                        "pho2_e2x5Left         := diPhoton.subLeadingPhoton.showerShapeVariables.e2x5Left",
                                        "pho2_e2x5Right        := diPhoton.subLeadingPhoton.showerShapeVariables.e2x5Right",
                                        "pho2_e2x5Top          := diPhoton.subLeadingPhoton.showerShapeVariables.e2x5Top",
                                        "pho2_e2x5Bottom       := diPhoton.subLeadingPhoton.showerShapeVariables.e2x5Bottom",
                                        "pho2_preshowerEnergy  := diPhoton.subLeadingPhoton.superCluster.preshowerEnergy",
                                        "pho2_preshowerEnergyPlane1 := diPhoton.subLeadingPhoton.superCluster.preshowerEnergyPlane1",
                                        "pho2_preshowerEnergyPlane2  := diPhoton.subLeadingPhoton.superCluster.preshowerEnergyPlane2",
                                        "pho2_cryEta           := diPhoton.subLeadingPhoton.cryEta",
                                        "pho2_cryPhi           := diPhoton.subLeadingPhoton.cryPhi",
                                        "pho2_ieta             := diPhoton.subLeadingPhoton.iEta",
                                        "pho2_iphi             := diPhoton.subLeadingPhoton.iPhi",

                                        ],
                             histograms=[]
                             )


# HLT filter + EE bad supercluster filter on data
# import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.eeBadScFilter
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
            ) )
else:
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
    ) )

 
process.p = cms.Path(process.hltHighLevel*
                     process.dataRequirements*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     (process.flashggTagSequence*process.systematicsTagSequences)*
                     process.flashggSystTagMerger*
                     process.tagsDumper
                     )


#printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",1e+3)
# call the customization
customize(process)
