import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# fullsim
#-------------------------------------------------
process = cms.Process("LowMass")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories = ['TopGenEvent']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## configure process options
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

# load sequences
process.load("TopQuarkAnalysis.TopPairBSM.TopAnalysis_sequences")

# setup path
#process.p = cms.Path( process.TopAnalysisMuFilter ) # with muonic generator filter
#process.p = cms.Path( process.TopAnalysisNoMuFilter ) # with generator filter on all but muonci decays
process.p = cms.Path( process.TopAnalysis ) # no generator filter at all
#process.TopAnalyzer.debug = cms.bool(True)
#process.TopAnalyzer.jetSource = cms.InputTag('selectedLayer1JetsJPT')

# change defaults
#process.BooTopHLTFilter.HLTPaths = cms.vstring() # do not filter
process.BooTopHLTFilter.HLTPaths = [] # do not filter

# source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
	'/store/user/samvel/TTbarJets-madgraph/Spring10-v2/e79177c217b45596f158355f13b94b69/ljmet_1.root'
    )
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1) )


