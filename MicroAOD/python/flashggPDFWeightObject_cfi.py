import FWCore.ParameterSet.Config as cms

flashggPDFWeightObject = cms.EDProducer('FlashggPDFWeightProducer',
		LHEEventTag = cms.InputTag('externalLHEProducer'),
                GenTag      = cms.InputTag('generator'),
		tag = cms.untracked.string("initrwgt"),
                isStandardSample = cms.untracked.bool(True),
		pdfset = cms.untracked.string("NNPDF30_lo_as_0130_nf_4.LHgrid"), # for non centrally produced H->gg samples
		doAlphasWeights = cms.untracked.bool(True),
		doScaleWeights  = cms.untracked.bool(True),
		nPdfEigWeights = cms.uint32(60),
		mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv')
	)
