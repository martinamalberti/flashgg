#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )



# input files
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        #data
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/SingleElectron/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015D-16Dec2015-v1/160127_024003/0000/myMicroAODOutputFile_1.root"
        # mc
        "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/160210_050006/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160210_045908/0000/myMicroAODOutputFile_1.root"
        ))

#output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)



# import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# load module to recompute photon id on-the-fly
process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")

# import flashgg customization to check if we have data, signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

# load syst producer
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
process.flashggDiPhotonSystematics.src = "flashggUpdatedIdMVADiPhotons"

# load appropriate scale and smearing bins 
process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official  tool instead  ????????????????
#useEGMTools(process)

from flashgg.Systematics.SystematicsCustomize import *

#customize.processType = 'data'

# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*") )
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele22_eta2p1_WPLoose_Gsf_v*") )
    customizePhotonSystematicsForMC(process)
    vpset = process.flashggDiPhotonSystematics.SystMethods
    newvpset = cms.VPSet()
    for pset in vpset:
        pset.NSigmas = cms.vint32() # only central value, no syst. shifts
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            newvpset+= [pset]
    process.flashggDiPhotonSystematics.SystMethods = newvpset        

printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])


# preselection Zee
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")
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
                                  "cosdphi            :=cos(leadingPhoton.phi-subLeadingPhoton.phi)",
                                  "pt                 :=pt",

                                  "pho1_pt            :=leadingPhoton.pt",
                                  "pho1_eta           :=leadingPhoton.superCluster.eta",
                                  "pho1_energy        :=leadingPhoton.energy",
                                  "pho1_rawEnergy     :=leadingPhoton.superCluster.rawEnergy",
                                  "pho1_eTrue         := ?leadingPhoton().hasMatchedGenPhoton()?leadingPhoton().matchedGenPhoton().energy():0",
                                  "pho1_r9            :=leadingPhoton.r9",
                                  "pho1_full5x5_r9    :=leadingPhoton.full5x5_r9",
                                  "pho1_sieie         :=leadingPhoton.sigmaIetaIeta",
                                  "pho1_full5x5_sieie :=leadingPhoton.full5x5_sigmaIetaIeta",
                                  "pho1_covieip       :=leadingPhoton.sieip",
                                  "pho1_etawidth      :=leadingPhoton.superCluster.etaWidth",
                                  "pho1_phiwidth      :=leadingPhoton.superCluster.phiWidth",
                                  "pho1_s4            :=leadingPhoton.s4",
                                  "pho1_phoIso        :=leadingPhoton.pfPhoIso03",
                                  "pho1_ChgIso        :=leadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "pho1_ChgIsoWorstVtx:=leadingPhoton.pfChgIsoWrtWorstVtx03",
                                  "pho1_esEffSigmaRR  :=leadingPhoton.esEffSigmaRR",
                                  "pho1_idmva         :=leadPhotonId",
                                  "pho1_sigmaEoE      :=leadingPhoton.sigEOverE",
                                  
                                  "pho2_pt            :=subLeadingPhoton.pt",
                                  "pho2_eta           :=subLeadingPhoton.superCluster.eta",
                                  "pho2_energy        :=subLeadingPhoton.energy",
                                  "pho2_rawEnergy     :=subLeadingPhoton.superCluster.rawEnergy",
                                  "pho2_eTrue         := ?subLeadingPhoton().hasMatchedGenPhoton()?subLeadingPhoton().matchedGenPhoton().energy():0",
                                  "pho2_r9            :=subLeadingPhoton.r9",
                                  "pho2_full5x5_r9    :=subLeadingPhoton.full5x5_r9",
                                  "pho2_sieie         :=subLeadingPhoton.sigmaIetaIeta",
                                  "pho2_full5x5_sieie :=subLeadingPhoton.full5x5_sigmaIetaIeta",
                                  "pho2_covieip       :=subLeadingPhoton.sieip",
                                  "pho2_etawidth      :=subLeadingPhoton.superCluster.etaWidth",
                                  "pho2_phiwidth      :=subLeadingPhoton.superCluster.phiWidth",
                                  "pho2_s4            :=subLeadingPhoton.s4",
                                  "pho2_phoIso        :=subLeadingPhoton.pfPhoIso03",
                                  "pho2_ChgIso        :=subLeadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "pho2_ChgIsoWorstVtx:=subLeadingPhoton.pfChgIsoWrtWorstVtx03",
                                  "pho2_esEffSigmaRR  :=subLeadingPhoton.esEffSigmaRR",
                                  "pho2_idmva         :=subLeadPhotonId",
                                  "pho2_sigmaEoE      :=subLeadingPhoton.sigEOverE"
                                  ],
                       histograms=[]
                       )



# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
    process.dataRequirements += process.eeBadScFilter
    


process.p = cms.Path(process.hltHighLevel*
                     process.dataRequirements*
                     process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggPreselectedDiPhotons*
                     process.diphotonDumper)



#printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.7e+3)
# call the customization
customize(process)
