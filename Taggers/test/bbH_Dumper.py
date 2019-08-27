import sys, os

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariablesHTXS,systematicVariablesHTXS
from flashgg.MetaData.MetaConditionsReader import *

# SYSTEMATICS SECTION
process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )


file_names = [
#"/store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/DoubleEG/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-Run2016B-17Jul2018_ver2-v1/190605_220256/0002/myMicroAODOutputFile_2363.root"
#"/store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/ttHJetToGG_M124_13TeV_amcatnloFXFX_madspin_pythia8/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/190715_071956/0000/myMicroAODOutputFile_4.root"
"file:/tmp/myMicroAODOutputFile_4.root"
]


process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring( file_names ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000)


systlabels = [""]
phosystlabels = []
metsystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []

from flashgg.MetaData.MetaConditionsReader import *

# set default options if needed
metaConditions = MetaConditionsReader("MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1.json")

ISDATA = False
if "DoubleEG" in file_names[0] or "EGamma" in file_names[0]:
  ISDATA = True

ISSIGNAL = False
signal_strings = ["ttHJetToGG", "ttHToGG", "THQ", "THW", "VBF", "GluGluHToGG", "VHToGG"]
if any([x in file_names[0] for x in signal_strings]):
  ISSIGNAL = True


### Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
if ISDATA:
    process.GlobalTag.globaltag = str(metaConditions['globalTags']['data'])
else:
    process.GlobalTag.globaltag = str(metaConditions['globalTags']['MC'])

from flashgg.Systematics.SystematicsCustomize import *
class customizer:
  def __init__(self, **kwargs):
    self.metaConditions = kwargs.get("metaConditions")
    self.processType = kwargs.get("processType")

customize = customizer(metaConditions = metaConditions, processType = "data" if ISDATA else "mc")
jetSystematicsInputTags = createStandardSystematicsProducers(process , customize)
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

# needed for 0th vertex from microAOD
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")
process.flashggDiPhotons.whichVertex = cms.uint32(0)
process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)

# Remove unneeded tags
process.flashggTagSequence.remove(process.flashggTTHDiLeptonTag)
process.flashggTagSequence.remove(process.flashggVBFTag)
process.flashggTagSequence.remove(process.flashggVHMetTag)
process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
process.flashggTagSequence.remove(process.flashggVHHadronicTag)
process.flashggTagSequence.remove(process.flashggUntagged)
process.flashggTagSequence.remove(process.flashggVBFMVA)
process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)
process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
process.flashggTagSequence.remove(process.flashggTHQLeptonicTag)

process.flashggTagSorter.TagPriorityRanges = cms.VPSet(   
        cms.PSet(TagName = cms.InputTag('flashggTTHLeptonicTag'))
)

# release selections
process.flashggTTHLeptonicTag.MVAThreshold = cms.vdouble(-999)
process.flashggTTHLeptonicTag.jetsNumberThreshold = cms.double(0.)
process.flashggTTHLeptonicTag.bjetsNumberThreshold = cms.double(0.)
process.flashggTTHLeptonicTag.jetPtThreshold = cms.double(20.)
process.flashggTTHLeptonicTag.MinNLep = cms.int32(0)
process.flashggTTHLeptonicTag.DiLeptonJetThreshold = cms.double(0)
process.flashggTTHLeptonicTag.DiLeptonbJetThreshold = cms.double(0)


# Or use the official tool instead
useEGMTools(process)

if ISSIGNAL:
  #print "Signal MC, so adding systematics and dZ"
  #customizeSystematicsForSignal(process)
  customizeSystematicsForBackground(process)
elif ISDATA:
  print "Data, so turn off all shifts and systematics, with some exceptions"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForData(process)
else:
  print "Background MC, so store mgg and central only"
  variablesToUse = minimalNonSignalVariables
  customizeSystematicsForBackground(process)

printSystematicInfo(process)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
hlt_paths = []
for dset in metaConditions["TriggerPaths"]:
  hlt_paths.extend(metaConditions["TriggerPaths"][dset])
for i in range(len(hlt_paths)):
  hlt_paths[i] = hlt_paths[i].encode("ascii")
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(hlt_paths))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if ISDATA: 
        #process.dataRequirements += process.hltHighLevel # already require triggers in microAOD production, don't need to require again
        process.dataRequirements += process.eeBadScFilter

process.genFilter = cms.Sequence()

# Met Filters
process.load('flashgg/Systematics/flashggMetFilters_cfi')
if ISDATA:
    metFilterSelector = "data" 
else:
    metFilterSelector = "mc"
process.flashggMetFilters.requiredFilterNames = cms.untracked.vstring([filter.encode("ascii") for filter in metaConditions["flashggMetFilters"][metFilterSelector]])

process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("merged_ntuple.root"))


## TAGS DUMPERS ##
from flashgg.Taggers.tagsDumpers_cfi import *

process.tthLeptonicTagDumper = createTagDumper("TTHLeptonicTag")
process.tthLeptonicTagDumper.dumpTrees = True

## define categories and associated objects to dump
import flashgg.Taggers.dumperConfigTools as cfgTools

dipho_variables=["dipho_sumpt      := diPhoton.sumPt",
                 "dipho_cosphi     := abs(cos(diPhoton.leadingPhoton.phi - diPhoton.subLeadingPhoton.phi))",
                 "diphoton_mva    := diPhotonMVA.result",
                 "mass             := diPhoton.mass",
                 "leadPt           := diPhoton.leadingPhoton.pt",
                 "leadEt           := diPhoton.leadingPhoton.et",
		 "leadEnergy       := diPhoton.leadingPhoton.energy",
                 "leadEta          := diPhoton.leadingPhoton.eta",
                 "leadPhi          := diPhoton.leadingPhoton.phi",
                 "lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
                 "lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
                 "lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
                 "lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
                 "leadR9           := diPhoton.leadingPhoton.full5x5_r9",
		 "leadGenMatch     := diPhoton.leadingPhoton.genMatchType",
		 "leadPtGen        := ? diPhoton.leadingPhoton.hasMatchedGenPhoton ? diPhoton.leadingPhoton.matchedGenPhoton.pt : 0",
		 "leadGendeltaR    := ? diPhoton.leadingPhoton.hasMatchedGenPhoton ? sqrt( pow((diPhoton.leadingPhoton.phi -  diPhoton.leadingPhoton.matchedGenPhoton.phi), 2) + pow((diPhoton.leadingPhoton.phi -  diPhoton.leadingPhoton.matchedGenPhoton.phi), 2)) : -999",
                 "subleadPt        := diPhoton.subLeadingPhoton.pt",
                 "subleadEt        := diPhoton.subLeadingPhoton.et",
                 "subleadEnergy    := diPhoton.subLeadingPhoton.energy",
                 "subleadEta       := diPhoton.subLeadingPhoton.eta",
                 "subleadPhi       := diPhoton.subLeadingPhoton.phi",
                 "sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
                 "sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
                 "sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
                 "sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
                 "subleadR9        := diPhoton.subLeadingPhoton.full5x5_r9",
                 "subleadGenMatch  := diPhoton.subLeadingPhoton.genMatchType",
		 "subleadPtGen     := ? diPhoton.subLeadingPhoton.hasMatchedGenPhoton ? diPhoton.subLeadingPhoton.matchedGenPhoton.pt : 0",
		 "subleadGendeltaR := ? diPhoton.subLeadingPhoton.hasMatchedGenPhoton ? sqrt( pow((diPhoton.subLeadingPhoton.phi -  diPhoton.subLeadingPhoton.matchedGenPhoton.phi), 2) + pow((diPhoton.subLeadingPhoton.phi -  diPhoton.subLeadingPhoton.matchedGenPhoton.phi), 2)) : -999",
                 "leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
                 "subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
                 "dipho_rapidity   := diPhoton.rapidity",
                 "vertex_idx       := diPhoton.vertexIndex",
                 "leadPassEVeto    := diPhoton.leadingPhoton.passElectronVeto",
                 "subleadPassEVeto := diPhoton.subLeadingPhoton.passElectronVeto",
                 "leadPixelSeed    := diPhoton.leadingPhoton.hasPixelSeed",
                 "subleadPixelSeed := diPhoton.subLeadingPhoton.hasPixelSeed",
]
## TTH leptonic ##
cfgTools.addCategories(process.tthLeptonicTagDumper,
                       ## categories definition
                       [("all","1",0)
                    ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       variables=dipho_variables+
		       ["n_ele    := electrons.size",
                        "ele1_pt  := ?(electrons.size>0)? electrons.at(0).pt : -1",
                        "ele2_pt  := ?(electrons.size>1)? electrons.at(1).pt : -1",
			"ele1_eta  := ?(electrons.size>0)? electrons.at(0).eta : -999",
                        "ele2_eta  := ?(electrons.size>1)? electrons.at(1).eta : -999",
			"ele1_phi  := ?(electrons.size>0)? electrons.at(0).phi : -999",
                        "ele2_phi  := ?(electrons.size>1)? electrons.at(1).phi : -999",
			"ele1_energy  := ?(electrons.size>0)? electrons.at(0).energy : -999",
                        "ele2_energy  := ?(electrons.size>1)? electrons.at(1).energy : -999",
                        "n_muons  := muons.size",
                        "muon1_pt := ?(muons.size>0)? muons.at(0).pt : -1",
                        "muon2_pt := ?(muons.size>1)? muons.at(1).pt : -1",
			"muon1_eta := ?(muons.size>0)? muons.at(0).eta : -999",
                        "muon2_eta := ?(muons.size>1)? muons.at(1).eta : -999",
			"muon1_phi := ?(muons.size>0)? muons.at(0).phi : -999",
                        "muon2_phi := ?(muons.size>1)? muons.at(1).phi : -999",
			"muon1_energy := ?(muons.size>0)? muons.at(0).energy : -999",
                        "muon2_energy := ?(muons.size>1)? muons.at(1).energy : -999",
			"tthMVA := mvaRes",
                        "jet_pt1                :=  ? jets.size()>0 ? jets[0].pt() : -100 ",
                        "jet_eta1               :=  ? jets.size()>0 ? jets[0].eta() : -100 ",
                        "jet_phi1               :=  ? jets.size()>0 ? jets[0].phi() : -100 ",
                        "jet_bdiscriminant1     :=  ? jets.size()>0 ? jets[0].bDiscriminator('pfDeepCSVJetTags:probb') : -100",
                        "jet_pt2                :=  ? jets.size()>1 ? jets[1].pt() : -100 ",
                        "jet_eta2               :=  ? jets.size()>1 ? jets[1].eta() : -100 ",
                        "jet_phi2               :=  ? jets.size()>1 ? jets[1].phi() : -100 ",
                        "jet_bdiscriminant2     :=  ? jets.size()>1 ? jets[1].bDiscriminator('pfDeepCSVJetTags:probb') : -100",
                        "jet_pt3                :=  ? jets.size()>2 ? jets[2].pt() : -100 ",
                        "jet_eta3               :=  ? jets.size()>2 ? jets[2].eta() : -100 ",
                        "jet_phi3               :=  ? jets.size()>2 ? jets[2].phi() : -100 ",
                        "jet_bdiscriminant3     :=  ? jets.size()>2 ? jets[2].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt4                :=  ? jets.size()>3 ? jets[3].pt() : -100 ",
                        "jet_eta4               :=  ? jets.size()>3 ? jets[3].eta() : -100 ",
                        "jet_phi4               :=  ? jets.size()>3 ? jets[3].phi() : -100 ",
                        "jet_bdiscriminant4     :=  ? jets.size()>3 ? jets[3].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt5                :=  ? jets.size()>4 ? jets[4].pt() : -100 ",
                        "jet_eta5               :=  ? jets.size()>4 ? jets[4].eta() : -100 ",
                        "jet_phi5               :=  ? jets.size()>4 ? jets[4].phi() : -100 ",
                        "jet_bdiscriminant5     :=  ? jets.size()>4 ? jets[4].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt6                :=  ? jets.size()>5 ? jets[5].pt() : -100 ",
                        "jet_eta6               :=  ? jets.size()>5 ? jets[5].eta() : -100 ",
                        "jet_phi6               :=  ? jets.size()>5 ? jets[5].phi() : -100 ",
                        "jet_bdiscriminant6     :=  ? jets.size()>5 ? jets[5].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt7                :=  ? jets.size()>6 ? jets[6].pt() : -100 ",
                        "jet_eta7               :=  ? jets.size()>6 ? jets[6].eta() : -100 ",
                        "jet_phi7               :=  ? jets.size()>6 ? jets[6].phi() : -100 ",
                        "jet_bdiscriminant7     :=  ? jets.size()>6 ? jets[6].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt8                :=  ? jets.size()>7 ? jets[7].pt() : -100 ",
                        "jet_eta8               :=  ? jets.size()>7 ? jets[7].eta() : -100 ",
                        "jet_phi8               :=  ? jets.size()>7 ? jets[7].phi() : -100 ",
                        "jet_bdiscriminant8     :=  ? jets.size()>7 ? jets[7].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
"jet_pt9                :=  ? jets.size()>8 ? jets[8].pt() : -100 ",
                        "jet_eta9               :=  ? jets.size()>8 ? jets[8].eta() : -100 ",
                        "jet_phi9               :=  ? jets.size()>8 ? jets[8].phi() : -100 ",
                        "jet_bdiscriminant9     :=  ? jets.size()>8 ? jets[8].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt10                :=  ? jets.size()>9 ? jets[9].pt() : -100 ",
                        "jet_eta10               :=  ? jets.size()>9 ? jets[9].eta() : -100 ",
                        "jet_phi10               :=  ? jets.size()>9 ? jets[9].phi() : -100 ",
                        "jet_bdiscriminant10     :=  ? jets.size()>9 ? jets[9].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt11                :=  ? jets.size()>10 ? jets[10].pt() : -100 ",
                        "jet_eta11               :=  ? jets.size()>10 ? jets[10].eta() : -100 ",
                        "jet_phi11               :=  ? jets.size()>10 ? jets[10].phi() : -100 ",
                        "jet_bdiscriminant11     :=  ? jets.size()>10 ? jets[10].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt12                :=  ? jets.size()>11 ? jets[11].pt() : -100 ",
                        "jet_eta12               :=  ? jets.size()>11 ? jets[11].eta() : -100 ",
                        "jet_phi12               :=  ? jets.size()>11 ? jets[11].phi() : -100 ",
                        "jet_bdiscriminant12     :=  ? jets.size()>11 ? jets[11].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt13                :=  ? jets.size()>12 ? jets[12].pt() : -100 ",
                        "jet_eta13               :=  ? jets.size()>12 ? jets[12].eta() : -100 ",
                        "jet_phi13               :=  ? jets.size()>12 ? jets[12].phi() : -100 ",
                        "jet_bdiscriminant13     :=  ? jets.size()>12 ? jets[12].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt14                :=  ? jets.size()>13 ? jets[13].pt() : -100 ",
                        "jet_eta14               :=  ? jets.size()>13 ? jets[13].eta() : -100 ",
                        "jet_phi14               :=  ? jets.size()>13 ? jets[13].phi() : -100 ",
                        "jet_bdiscriminant14     :=  ? jets.size()>13 ? jets[13].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",
                        "jet_pt15                :=  ? jets.size()>14 ? jets[14].pt() : -100 ",
                        "jet_eta15               :=  ? jets.size()>14 ? jets[14].eta() : -100 ",
                        "jet_phi15               :=  ? jets.size()>14 ? jets[14].phi() : -100 ",
                        "jet_bdiscriminant15     :=  ? jets.size()>14 ? jets[14].bDiscriminator('pfDeepCSVJetTags:probb') : -100 ",

			"jet_bbdiscriminant1  := ?(jets.size>0)? jets.at(0).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant2  := ?(jets.size>1)? jets.at(1).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant3  := ?(jets.size>2)? jets.at(2).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant4  := ?(jets.size>3)? jets.at(3).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant5  := ?(jets.size>4)? jets.at(4).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant6  := ?(jets.size>5)? jets.at(5).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant7  := ?(jets.size>6)? jets.at(6).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant8  := ?(jets.size>7)? jets.at(7).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant9  := ?(jets.size>8)? jets.at(8).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant10  := ?(jets.size>9)? jets.at(9).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant11  := ?(jets.size>10)? jets.at(10).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant12  := ?(jets.size>11)? jets.at(11).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant13  := ?(jets.size>12)? jets.at(12).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant14  := ?(jets.size>13)? jets.at(13).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",
                        "jet_bbdiscriminant15  := ?(jets.size>14)? jets.at(14).bDiscriminator('pfDeepCSVJetTags:probbb') : -1",      
			"jet_cdiscriminant1  := ?(jets.size>0)? jets.at(0).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant2  := ?(jets.size>1)? jets.at(1).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant3  := ?(jets.size>2)? jets.at(2).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant4  := ?(jets.size>3)? jets.at(3).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant5  := ?(jets.size>4)? jets.at(4).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant6  := ?(jets.size>5)? jets.at(5).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant7  := ?(jets.size>6)? jets.at(6).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant8  := ?(jets.size>7)? jets.at(7).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant9  := ?(jets.size>8)? jets.at(8).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant10  := ?(jets.size>9)? jets.at(9).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant11  := ?(jets.size>10)? jets.at(10).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant12  := ?(jets.size>11)? jets.at(11).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant13  := ?(jets.size>12)? jets.at(12).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant14  := ?(jets.size>13)? jets.at(13).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
                        "jet_cdiscriminant15  := ?(jets.size>14)? jets.at(14).bDiscriminator('pfDeepCSVJetTags:probc') : -1",
			"jet_udsgdiscriminant1  := ?(jets.size>0)? jets.at(0).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant2  := ?(jets.size>1)? jets.at(1).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant3  := ?(jets.size>2)? jets.at(2).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant4  := ?(jets.size>3)? jets.at(3).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant5  := ?(jets.size>4)? jets.at(4).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant6  := ?(jets.size>5)? jets.at(5).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant7  := ?(jets.size>6)? jets.at(6).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant8  := ?(jets.size>7)? jets.at(7).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant9  := ?(jets.size>8)? jets.at(8).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant10  := ?(jets.size>9)? jets.at(9).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant11  := ?(jets.size>10)? jets.at(10).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant12  := ?(jets.size>11)? jets.at(11).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant13  := ?(jets.size>12)? jets.at(12).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant14  := ?(jets.size>13)? jets.at(13).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",
                        "jet_udsgdiscriminant15  := ?(jets.size>14)? jets.at(14).bDiscriminator('pfDeepCSVJetTags:probudsg') : -1",    
			"jet_energy1  := ?(jets.size>0)? jets.at(0).energy : -1",
                        "jet_energy2  := ?(jets.size>1)? jets.at(1).energy : -1",
                        "jet_energy3  := ?(jets.size>2)? jets.at(2).energy : -1",
                        "jet_energy4  := ?(jets.size>3)? jets.at(3).energy : -1",
                        "jet_energy5  := ?(jets.size>4)? jets.at(4).energy : -1",
                        "jet_energy6  := ?(jets.size>5)? jets.at(5).energy : -1",
                        "jet_energy7  := ?(jets.size>6)? jets.at(6).energy : -1",
                        "jet_energy8  := ?(jets.size>7)? jets.at(7).energy : -1",
                        "jet_energy9  := ?(jets.size>8)? jets.at(8).energy : -1",
                        "jet_energy10  := ?(jets.size>9)? jets.at(9).energy : -1",
                        "jet_energy11  := ?(jets.size>10)? jets.at(10).energy : -1",
                        "jet_energy12  := ?(jets.size>11)? jets.at(11).energy : -1",
                        "jet_energy13  := ?(jets.size>12)? jets.at(12).energy : -1",
                        "jet_energy14  := ?(jets.size>13)? jets.at(13).energy : -1",
                        "jet_energy15  := ?(jets.size>14)? jets.at(14).energy : -1",


                   ],
                       ## histograms
                       histograms=[]
)



process.p = cms.Path(process.dataRequirements
                     #*process.flashggMetFilters
                     #*process.flashggDiPhotons # needed for 0th vertex from microAOD
                     *process.flashggUpdatedIdMVADiPhotons
                     *process.flashggDiPhotonSystematics
                     *process.flashggMetSystematics
                     *process.flashggMuonSystematics
                     *process.flashggElectronSystematics
                     *(process.flashggUnpackedJets*process.jetSystematicsSequence)
                     *(process.flashggTagSequence*process.systematicsTagSequences)
                     *process.flashggSystTagMerger
                     *process.flashggTagSequence
                     *process.flashggTagTester
                     *process.tthLeptonicTagDumper
)
