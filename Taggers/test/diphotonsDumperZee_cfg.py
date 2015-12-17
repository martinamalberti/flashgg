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



process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   flashggDiPhotonSystematics = cms.PSet(initialSeed = cms.untracked.uint32(664)),
                                                   flashggElectronSystematics = cms.PSet(initialSeed = cms.untracked.uint32(11)),
                                                   flashggMuonSystematics = cms.PSet(initialSeed = cms.untracked.uint32(13)),
                                                   flashggTagSystematics = cms.PSet(initialSeed = cms.untracked.uint32(999)),
                                                   flashggRandomizedPerPhotonDiPhotons = cms.PSet(initialSeed = cms.untracked.uint32(281765313))
                                                   )


# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()


# import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel

# load randomizer for smearings
process.load("flashgg.MicroAOD.flashggRandomizedPerPhotonDiPhotonProducer_cff")

# load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
process.flashggDiPhotonSystematics.src = cms.InputTag("flashggRandomizedPerPhotonDiPhotons") 

# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*") )
    newvpset = cms.VPSet()
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if pset.Label.value().count("Scale"):
            pset.NoCentralShift = cms.bool(False) # Turn on central shift for data (it is off for MC)                                                                              
            pset.NSigmas = cms.vint32() # Do not perform shift
            newvpset += [pset]
    process.flashggDiPhotonSystematics.SystMethods = newvpset
else:
    print 'mc'
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WP75_Gsf_v*") )
    newvpset = cms.VPSet()
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if pset.Label.value().count("Smear"):
            pset.NSigmas = cms.vint32() # Do not perform shift
            newvpset += [pset]
    process.flashggDiPhotonSystematics.SystMethods = newvpset
    


# preselection Zee
process.load("flashgg.MicroAOD.flashggPreselectedDiPhotons_cfi")
process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"
process.flashggPreselectedDiPhotons.src = cms.InputTag("flashggDiPhotonSystematics")


# diphoton dumper
process.load("flashgg.Taggers.diphotonDumper_cfi") ##  import diphotonDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.diphotonDumper.src = "flashggPreselectedDiPhotons"

process.diphotonDumper.dumpTrees = True
process.diphotonDumper.dumpWorkspace = False
process.diphotonDumper.quietRooFit = True


# split tree, histogram and datasets by process
#process.diphotonDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL_$SUBCAT"

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
                                  "vtxProbMVA         :=vtxProbMVA",
                                  "cosdphi            := cos(leadingPhoton.phi-subLeadingPhoton.phi)",
                                  
                                  "pho1_pt            :=leadingPhoton.pt",
                                  "pho1_eta           :=leadingPhoton.superCluster.eta",
                                  "pho1_r9            :=leadingPhoton.r9",
                                  "pho1_energy        :=leadingPhoton.energy",
                                  "pho1_rawEnergy     :=leadingPhoton.superCluster.rawEnergy",
                                  "pho1_eTrue         := ?leadingPhoton().hasMatchedGenPhoton()?leadingPhoton().matchedGenPhoton().energy():0",
                                  "pho1_sieie         :=leadingPhoton.full5x5_sigmaIetaIeta",
                                  "pho1_covieip       :=leadingPhoton.sieip",
                                  "pho1_etawidth      :=leadingPhoton.superCluster.etaWidth",
                                  "pho1_phiwidth      :=leadingPhoton.superCluster.phiWidth",
                                  "pho1_s4ratio       :=leadingPhoton.s4",
                                  "pho1_phoIso        :=leadingPhoton.pfPhoIso03",
                                  "pho1_ChgIso        :=leadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "pho1_esEffSigmaRR  :=leadingPhoton.esEffSigmaRR",
                                  "pho1_idmva         :=leadPhotonId",
                                  "pho1_sigmaEoE      :=leadingPhoton.sigEOverE",
                                  
                                  "pho2_pt            :=subLeadingPhoton.pt",
                                  "pho2_eta           :=subLeadingPhoton.superCluster.eta",
                                  "pho2_r9            :=subLeadingPhoton.r9",
                                  "pho2_energy        :=subLeadingPhoton.energy",
                                  "pho2_rawEnergy     :=subLeadingPhoton.superCluster.rawEnergy",
                                  "pho2_eTrue         := ?subLeadingPhoton().hasMatchedGenPhoton()?subLeadingPhoton().matchedGenPhoton().energy():0",
                                  "pho2_sieie         :=subLeadingPhoton.full5x5_sigmaIetaIeta",
                                  "pho2_covieip       :=subLeadingPhoton.sieip",
                                  "pho2_etawidth      :=subLeadingPhoton.superCluster.etaWidth",
                                  "pho2_phiwidth      :=subLeadingPhoton.superCluster.phiWidth",
                                  "pho2_s4ratio       :=subLeadingPhoton.s4",
                                  "pho2_phoIso        :=subLeadingPhoton.pfPhoIso03",
                                  "pho2_ChgIso        :=subLeadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "pho2_esEffSigmaRR  :=subLeadingPhoton.esEffSigmaRR",
                                  "pho2_idmva         :=subLeadPhotonId",
                                  "pho2_sigmaEoE      :=subLeadingPhoton.sigEOverE"
                                  ],
                       histograms=[]
                       )


process.p = cms.Path(process.hltHighLevel*
                     process.flashggRandomizedPerPhotonDiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggPreselectedDiPhotons*
                     process.diphotonDumper)


# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.5e+3)
# call the customization
customize(process)
