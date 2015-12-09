#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' # keep updated for JEC
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )


# input files
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-BetaV7-25ns/Spring15BetaV7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15-ReMiniAOD-BetaV7-25ns-Spring15BetaV7-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151021_151505/0000/myMicroAODOutputFile_10.root"
        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-FinalPrompt-BetaV7-25ns/Spring15BetaV7/DoubleEG/RunIISpring15-FinalPrompt-BetaV7-25ns-Spring15BetaV7-v0-Run2015D-PromptReco-v4/151124_234634/0000/myMicroAODOutputFile_1.root"
        )
                            )

#output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)



# customize job config
from flashgg.MetaData.JobConfig import customize
customize.parse()

customize.setDefault("maxEvents",1000)
customize.setDefault("targetLumi",1.e+3)
#customize.setDefault("processType","data")

#ANALYSIS

# setup to load scale and smearings
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   flashggDiPhotonScaleAndSmearings = cms.PSet(initialSeed = cms.untracked.uint32(664))
                                                   )

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel

if customize.processType == 'data':
    print 'data' 
    from flashgg.Systematics.flashggDiPhotonScale_cfi import flashggDiPhotonScale
    process.flashggDiPhotonScaleAndSmearings = flashggDiPhotonScale.clone()
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*") )
    #process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*") )
else:
    print 'mc'
    from flashgg.Systematics.flashggDiPhotonSmearings_cfi import flashggDiPhotonSmearings
    process.flashggDiPhotonScaleAndSmearings = flashggDiPhotonSmearings.clone()
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WP75_Gsf_v*") )


# preselction Zee
process.load("flashgg.MicroAOD.flashggPreselectedDiPhotons_cfi")
process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"
process.flashggPreselectedDiPhotons.src = cms.InputTag("flashggDiPhotonScaleAndSmearings")


# diphoton dumper
process.load("flashgg.Taggers.diphotonDumper_cfi") ##  import diphotonDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.diphotonDumper.src = "flashggPreselectedDiPhotons"

process.diphotonDumper.dumpTrees = True
process.diphotonDumper.dumpWorkspace = False
process.diphotonDumper.quietRooFit = True


# split tree, histogram and datasets by process
#process.diphotonDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"
## do not split by process
process.diphotonDumper.nameTemplate = "minitree_$SQRTS_$LABEL_$SUBCAT"


# interestng categories 
cfgTools.addCategories(process.diphotonDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [
                        ("All","leadingPhoton.pt>0",0)
                       ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=["mass               :=mass",
                                  "pho1_pt            :=leadingPhoton.pt",
                                  "pho2_pt            :=subLeadingPhoton.pt",
                                  "pho1_eta           :=leadingPhoton.superCluster.eta",
                                  "pho2_eta           :=subLeadingPhoton.superCluster.eta",
                                  "pho1_r9            :=leadingPhoton.r9",
                                  "pho2_r9            :=subLeadingPhoton.r9",
                                  "pho1_energy        :=leadingPhoton.energy",
                                  "pho2_energy        :=subLeadingPhoton.energy",
                                  "pho1_rawEnergy     :=leadingPhoton.superCluster.rawEnergy",
                                  "pho2_rawEnergy     :=subLeadingPhoton.superCluster.rawEnergy",
                                  "pho1_etrue         := ?leadingPhoton().hasMatchedGenPhoton()?leadingPhoton().matchedGenPhoton().energy():0",
                                  "pho2_etrue         := ?subLeadingPhoton().hasMatchedGenPhoton()?subLeadingPhoton().matchedGenPhoton().energy():0"
                                  ],
                       histograms=[]
                       )



if customize.processType == "data":
    process.p = cms.Path(process.hltHighLevel*
                         process.flashggDiPhotonScaleAndSmearings*
                         process.flashggPreselectedDiPhotons*
                         process.diphotonDumper)
else:
    process.p = cms.Path(process.hltHighLevel*
                         process.flashggDiPhotonScaleAndSmearings*
                         process.flashggPreselectedDiPhotons*
                         process.diphotonDumper)




customize(process)


