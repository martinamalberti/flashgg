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
        "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/SingleElectron/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015C_25ns-16Dec2015-v1/160127_023910/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/SingleElectron/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015D-16Dec2015-v1/160127_024003/0000/myMicroAODOutputFile_1.root"
        # mc
        #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/160210_050006/0000/myMicroAODOutputFile_1.root"
        #"/store/group/phys_higgs/cmshgg/malberti/flashgg/RunIIFall15DR76-WJets/1_3_0/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-WJets-1_3_0-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/160420_104732/0000/myMicroAODOutputFile_1.root"
        ))

#output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
                                   )



# import trigger filter
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# load module to correct photon shower shapes on-the-fly  ###  for single photons!!!!!!!
process.load("flashgg/Taggers/flashggUpdatedShowerShapesPhotons_cfi")

# import flashgg customization to check if we have data, signal or background
from flashgg.MetaData.JobConfig import customize
customize.parse()

# load syst producer - need to apply to single photons
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")

for isyst in [process.MCScaleHighR9EB, process.MCScaleLowR9EB, process.MCScaleHighR9EE, process.MCScaleLowR9EE, process.MCSmearHighR9EE, process.MCSmearLowR9EE, process.MCSmearHighR9EB, process.MCSmearLowR9EB, process.SigmaEOverESmearing]:
    isyst.MethodName = isyst.PhotonMethodName

for isyst in [process.MCScaleHighR9EB_EGM, process.MCScaleLowR9EB_EGM, process.MCScaleHighR9EE_EGM, process.MCScaleLowR9EE_EGM, process.MCSmearHighR9EE_EGM, process.MCSmearLowR9EE_EGM, process.MCSmearHighR9EB_EGM, process.MCSmearLowR9EB_EGM, process.SigmaEOverESmearing_EGM]:
    isyst.MethodName = isyst.PhotonMethodName

process.flashggDiPhotonSystematics = cms.EDProducer('FlashggPhotonSystematicProducer',
                #src = cms.InputTag("flashggRandomizedPhotons"),
                src = cms.InputTag("flashggUpdatedShowerShapesPhotons"),
                SystMethods2D = cms.VPSet(),
                SystMethods = cms.VPSet(
                    process.MCScaleHighR9EB,
                    process.MCScaleLowR9EB,
                    process.MCScaleHighR9EE,
                    process.MCScaleLowR9EE,
                    process.MCSmearHighR9EE,
                    process.MCSmearLowR9EE,
                    process.MCSmearHighR9EB,
                    process.MCSmearLowR9EB,
                    process.SigmaEOverESmearing
                    )
)



# customize photon systematics for MC and data
from flashgg.Systematics.SystematicsCustomize import *

process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")
#useEGMTools(process) # non funziona...

####customize.processType = 'data' # for test

# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    process.hltHighLevel= hltHighLevel.clone(
        HLTPaths = cms.vstring(
            "HLT_Ele23_WPLoose_Gsf_v*",
            "HLT_Ele27_WPLoose_Gsf_v*",
            "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
            "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
            ),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False) # throw exception on unknown path names 
        )
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
    process.hltHighLevel= hltHighLevel.clone(
        HLTPaths = cms.vstring(
            "HLT_Ele23_WPLoose_Gsf_v*",
            "HLT_Ele27_WPLoose_Gsf_v*",
            "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
            "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
            ),
        andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
        throw = cms.bool(False) # throw exception on unknown path names 
        )
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
        pset.ApplyCentralValue = cms.bool(False) # no central value
        if ( pset.Label.value().count("MCSmear") or pset.Label.value().count("SigmaEOverESmearing")):
            pset.ApplyCentralValue = cms.bool(True)
            newvpset2D+= [pset]
    process.flashggDiPhotonSystematics.SystMethods2D = newvpset2D       

print 'syst 1D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods])
print 'syst 2D'
printSystematicVPSet([process.flashggDiPhotonSystematics.SystMethods2D])


# wenu producer + dumper
process.load("flashgg.Taggers.flashggWenu_cfi")
process.flashggWenu.PhotonTag=cms.InputTag('flashggDiPhotonSystematics')

process.load("flashgg.Taggers.wenuDumper_cfi") ##  import photonjetDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools

process.wenuDumper.src = "flashggWenu"
process.wenuDumper.maxCandPerEvent = -1 # take them all
process.wenuDumper.dumpTrees = True
process.wenuDumper.dumpWorkspace = False
process.wenuDumper.quietRooFit = True
process.wenuDumper.nameTemplate ="tree_$SQRTS_$LABEL"


process.wenuDumper.globalVariables.addTriggerBits = cms.PSet(
    tag = cms.InputTag("TriggerResults::HLT"),
    bits = cms.vstring(
    "HLT_Ele23_WPLoose_Gsf_v",
    "HLT_Ele27_WPLoose_Gsf_v",
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
    "HLT_Ele22_eta2p1_WPLoose_Gsf_v"
    )
)

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
                       ## if different variables wanted for different categories, can add  categorie one by one with cfgTools.addCategory
                       variables=["ele1_energy        :=photon.energy",  # photon regression
                                  "ele1_et            :=photon.pt", # photon regression
                                  "ele1_rawEnergy     :=electron.superCluster.rawEnergy",
                                  "ele1_esEnergy      :=electron.superCluster.preshowerEnergy",
                                  "ele1_scEta         :=electron.superCluster.eta",
                                  "ele1_scPhi         :=electron.superCluster.phi",
                                  "ele1_tkP           :=electron.trackMomentumAtVtx().R()",
                                  "ele1_tkPt          :=electron.trackMomentumAtVtx().rho()",
                                  "ele1_eta           :=electron.eta()",
                                  "ele1_phi           :=electron.phi()",
                                  #"ele1_eTrue         := ?electron.genLepton()?electron.genLepton().energy():0",
                                  "ele1_r9            :=photon.r9",
                                  "ele1_full5x5_r9    :=photon.full5x5_r9",
                                  "met                := met.pt",
                                  "metPhi             := met.phi",
                                  "mex                := met.px",
                                  "mey                := met.py"
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
                     process.flashggUpdatedShowerShapesPhotons*
                     process.flashggDiPhotonSystematics*
                     process.flashggWenu*
                     process.wenuDumper
                     )



#printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.7e+3)
# call the customization
customize(process)
