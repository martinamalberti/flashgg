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
        ))

#output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
                                   )



# import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# load module to recompute photon id on-the-fly  ### FIXME for single photons!!!!!!!
#process.load("flashgg/Taggers/flashggUpdatedIdMVADiPhotons_cfi")

# import flashgg customization to check if we have data, signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

# load syst producer - need to apply to single photons
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
#process.flashggDiPhotonSystematics.src = "flashggUpdatedIdMVADiPhotons"

from flashgg.Systematics.flashggDiPhotonSystematics_cfi import *
for isyst in [process.MCScaleHighR9EB, process.MCScaleLowR9EB, process.MCScaleHighR9EE, process.MCScaleLowR9EE, process.MCSmearHighR9EE, process.MCSmearLowR9EE, process.MCSmearHighR9EB, process.MCSmearLowR9EB, process.SigmaEOverESmearing]:
    isyst.MethodName = isyst.PhotonMethodName

process.flashggDiPhotonSystematics = cms.EDProducer('FlashggPhotonSystematicProducer',
                src = cms.InputTag("flashggRandomizedPhotons"),
                SystMethods2D = cms.VPSet(),
                SystMethods = cms.VPSet(
                    MCScaleHighR9EB,
                    MCScaleLowR9EB,
                    MCScaleHighR9EE,
                    MCScaleLowR9EE,
                    MCSmearHighR9EE,
                    MCSmearLowR9EE,
                    MCSmearHighR9EB,
                    MCSmearLowR9EB
                    )
)


# customize photon systematics for MC and data
from flashgg.Systematics.SystematicsCustomize import *

# use EGM tool
#useEGMTools(process)

#customize.processType = 'data' # for test

# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v*") ) #FIXME
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele22_eta2p1_WPLoose_Gsf_v*") )
    customizePhotonSystematicsForMC(process)
    #syst (1D) 
    vpset   = process.flashggDiPhotonSystematics.SystMethods
    newvpset = cms.VPSet()
    for pset in vpset:
        pset.NSigmas = cms.vint32() # no up/down syst shifts
        pset.ApplyCentralValue = cms.bool(False) # no central value
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
        newvpset+= [pset]
    process.flashggDiPhotonSystematics.SystMethods = newvpset        
    #syst (2D) : smearings with EGMTool
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

print 'ciao'

# preselection Zee
#process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")
#process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"
#process.flashggPreselectedDiPhotons.src = cms.InputTag("flashggDiPhotonSystematics")


# wenu producer + dumper
process.load("flashgg.Taggers.flashggWenu_cfi")
process.flashggWenu.PhotonTag=cms.InputTag('flashggDiPhotonSystematics')

process.load("flashgg.Taggers.wenuDumper_cfi") ##  import photonjetDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.wenuDumper.src = "flashggWenu"
process.wenuDumper.maxCandPerEvent = -1 # take them all
#process.wenuDumper.splitLumiWeight=cms.untracked.bool(True)
process.wenuDumper.dumpTrees = True
process.wenuDumper.dumpWorkspace = False
process.wenuDumper.quietRooFit = True
#process.wenuDumper.nameTemplate ="$PROCESS_$SQRTS_$LABEL"
process.wenuDumper.nameTemplate ="tree_$SQRTS_$LABEL"

## define categories and associated objects to dump
cfgTools.addCategory(process.wenuDumper,
                     "Reject",
                     "0",
                      -1 ## if nSubcat is -1 do not store anythings
                     )

# interestng categories 
cfgTools.addCategories(process.wenuDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [
                        ("All","1",0)
                       ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=["pho1_pt            :=photon.pt",
                                  "pho1_eta           :=photon.superCluster.eta",
                                  "pho1_energy        :=photon.energy",
                                  "pho1_rawEnergy     :=photon.superCluster.rawEnergy",
                                  "pho1_eTrue         := ?photon().hasMatchedGenPhoton()?photon().matchedGenPhoton().energy():0",
                                  "pho1_r9            :=photon.r9",
                                  "pho1_full5x5_r9    :=photon.full5x5_r9",
                                  #"pho1_etawidth      :=photon.superCluster.etaWidth",
                                  #"pho1_phiwidth      :=photon.superCluster.phiWidth",
                                  #"pho1_s4            :=photon.s4",
                                  #"pho1_phoIso        :=photon.pfPhoIso03",
                                  #"pho1_ChgIso        :=photon.pfChgIsoWrtChosenVtx03",
                                  #"pho1_ChgIsoWorstVtx:=photon.pfChgIsoWrtWorstVtx03",
                                  #"pho1_esEffSigmaRR  :=photon.esEffSigmaRR",
                                  #"pho1_idmva         :=photon.leadPhotonId",
                                  #"pho1_sigmaEoE      :=photon.sigEOverE",
                                  "met  := met.pt"
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
                     #process.flashggUpdatedIdMVADiPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggWenu*
                     process.wenuDumper
                     )



#printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents",100)
customize.setDefault("targetLumi",2.7e+3)
# call the customization
customize(process)
