import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.wenuDumpConfig_cff import wenuDumpConfig

wenuDumper = cms.EDAnalyzer('CutBasedWenuDumper',
                                **wenuDumpConfig.parameters_()
                                )


