#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

from Configuration.AlCa.GlobalTag import GlobalTag
import os
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v13')
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_v11')
else:
    raise Exception,"The default setup for microAODstd.py does not support releases other than 76X and 80X"


# flags
doUpdatedIdMVADiPhotons = False
doDiphotonsSystematics = True


# input files
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        #data
        "/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring16_v1_p3/diphotons_80_v1/SingleElectron/EXOSpring16_v1_p3-diphotons_80_v1-v0-Run2016B-PromptReco-v1/160519_094908/0000/diphotonsMicroAOD_10.root" 
        #mc
        #"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring16_v1/diphotons_80_v1/DYToEE_NNPDF30_13TeV-powheg-pythia8/EXOSpring16_v1-diphotons_80_v1-v0-RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/160503_231305/0000/diphotonsMicroAOD_1.root"
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
if (doUpdatedIdMVADiPhotons):
    process.flashggDiPhotonSystematics.src = "flashggUpdatedIdMVADiPhotons"
else:
    process.flashggDiPhotonSystematics.src = "flashggDiPhotons"

from flashgg.Systematics.SystematicsCustomize import *

# load appropriate scale and smearing bins 
#process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official  tool instead  ????????????????
useEGMTools(process)

customize.processType = 'data' # for test

# if data, apply only energy scale corrections, if MC apply only energy smearings
if customize.processType == 'data':
    print 'data' 
    customizePhotonSystematicsForData(process)    # only central value, no syst. shifts 
else:
    print 'mc'
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


# preselection Zee
process.load("flashgg.Taggers.flashggPreselectedDiPhotons_cfi")
process.flashggPreselectedDiPhotons.variables[-1] = "-(passElectronVeto - 1)"
if (doDiphotonsSystematics): 
    process.flashggPreselectedDiPhotons.src = cms.InputTag("flashggDiPhotonSystematics")
else:
    process.flashggPreselectedDiPhotons.src = cms.InputTag("flashggDiPhotons")


## diphoton dumper
process.load("flashgg.Taggers.diphotonDumper_cfi") ##  import diphotonDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools
process.diphotonDumper.src = "flashggPreselectedDiPhotons"
process.diphotonDumper.dumpTrees = True
process.diphotonDumper.dumpWorkspace = False
process.diphotonDumper.quietRooFit = True
process.diphotonDumper.nameTemplate = "tree_$SQRTS_$LABEL_$SUBCAT"

process.diphotonDumper.globalVariables.addTriggerBits = cms.PSet(
    tag = cms.InputTag("TriggerResults::HLT"),
    bits = cms.vstring(
        #DoubleEG
        #"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*"
        #SingleElectron
        "HLT_Ele27_WPLoose_Gsf_v1",
        "HLT_Ele27_WPTight_Gsf_v1",
        #"HLT_Ele27_eta2p1_WPLoose_Gsf_v",
        #"HLT_Ele22_eta2p1_WPLoose_Gsf_v"
    )
)


## define categories and associated objects to dump
cfgTools.addCategory(process.diphotonDumper,
                     "Reject",
                     "0",
                      -1 ## if nSubcat is -1 do not store anythings
                     )

# interestng categories 
cfgTools.addCategories(process.diphotonDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [
                        ("All","1",0)
                       ],
                       ## variables to be dumped in trees/datasets. Same variables for all categories
                       ## if different variables wanted for different categories, can add categorie one by one with cfgTools.addCategory
                       variables=["mass               :=mass",
                                  "pt                 :=pt",
                                  #"vtxprobMVA         :=diPhotonMVA.vtxprob",
                                  #"cosdphi            :=diPhotonMVA.CosPhi",
                                  #"diphotonMVA        :=diPhotonMVA.result",

                                  "pho1_et            :=leadingPhoton.pt",
                                  "pho1_energy        :=leadingPhoton.energy",
                                  "pho1_rawEnergy     :=leadingPhoton.superCluster.rawEnergy",
                                  "pho1_esEnergy      :=leadingPhoton.superCluster.preshowerEnergy",
                                  "pho1_eTrue         := ?leadingPhoton().hasMatchedGenPhoton()?leadingPhoton().matchedGenPhoton().energy():0",
                                  "pho1_scEta         :=leadingPhoton.superCluster.eta",
                                  "pho1_scPhi         :=leadingPhoton.superCluster.phi",
                                  "pho1_eta           :=leadingPhoton.eta",
                                  "pho1_phi           :=leadingPhoton.phi",
                                  "pho1_r9            :=leadingPhoton.r9",
                                  "pho1_full5x5_r9    :=leadingPhoton.full5x5_r9",
                                  "pho1_idmva         :=leadPhotonId",
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
                                  "pho1_sigmaEoE      :=leadingPhoton.sigEOverE",

                                  "pho2_et            :=subLeadingPhoton.pt",
                                  "pho2_energy        :=subLeadingPhoton.energy",
                                  "pho2_rawEnergy     :=subLeadingPhoton.superCluster.rawEnergy",
                                  "pho2_esEnergy      :=subLeadingPhoton.superCluster.preshowerEnergy",
                                  "pho2_eTrue         := ?leadingPhoton().hasMatchedGenPhoton()?leadingPhoton().matchedGenPhoton().energy():0",
                                  "pho2_scEta         :=subLeadingPhoton.superCluster.eta",
                                  "pho2_scPhi         :=subLeadingPhoton.superCluster.phi",
                                  "pho2_eta           :=subLeadingPhoton.eta",
                                  "pho2_phi           :=subLeadingPhoton.phi",
                                  "pho2_r9            :=subLeadingPhoton.r9",
                                  "pho2_full5x5_r9    :=subLeadingPhoton.full5x5_r9",
                                  "pho2_idmva         :=subLeadPhotonId",
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
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
            #DoubleEG
            #"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*",
            #SingleEG
            "HLT_Ele27_WPLoose_Gsf_v*", # 7_6_X, 8_0_X
            "HLT_Ele27_WPLoose_Gsf_v*",
            #"HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
            #"HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
            ) )
#no trigger available in MC 80X
else:
    process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
            #DoubleEG
            #"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v*",
            #SingleEG
#            "HLT_Ele27_WPLoose_Gsf_v*" # 7_6_X
#            "HLT_Ele27_WPLoose_Gsf_v*",
#            "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
#            "HLT_Ele22_eta2p1_WPLoose_Gsf_v*"
            ) )


if (doDiphotonsSystematics and doUpdatedIdMVADiPhotons):
    process.p = cms.Path(process.hltHighLevel*
                         process.dataRequirements*
                         process.flashggUpdatedIdMVADiPhotons*
                         process.flashggDiPhotonSystematics*
                         process.diphotonDumper)
elif (doDiphotonsSystematics and not doUpdatedIdMVADiPhotons):
    process.p = cms.Path(process.hltHighLevel*
                         process.dataRequirements*
                         #process.flashggUpdatedIdMVADiPhotons*
                         process.flashggDiPhotonSystematics*
                         process.flashggPreselectedDiPhotons*
                         process.diphotonDumper)
else:
     process.p = cms.Path(process.hltHighLevel*
                          process.dataRequirements*
                          #process.flashggUpdatedIdMVADiPhotons* 
                          #process.flashggDiPhotonSystematics*
                          process.flashggPreselectedDiPhotons*
                          process.diphotonDumper)

#printSystematicInfo(process)

# set default options if needed
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",2.7e+3)
# call the customization
customize(process)
