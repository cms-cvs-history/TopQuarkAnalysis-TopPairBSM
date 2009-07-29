import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# fullsim
#-------------------------------------------------
process = cms.Process("LowMass")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

## configure process options
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")
## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_31X_V2::All')
## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

# load sequences
process.load("TopQuarkAnalysis.TopPairBSM.TopAnalysis_sequences")

# setup path
#process.p = cms.Path( process.TopAnalysisMuFilter ) # with muonic generator filter
#process.p = cms.Path( process.TopAnalysisNoMuFilter ) # with generator filter on all but muonci decays
process.p = cms.Path( process.TopAnalysis ) # no generator filter at all
process.TopAnalyzer.debug = cms.bool(True)
#process.TopAnalyzer.jetSource      = cms.InputTag('selectedLayer1JetsJPT')

# change defaults
process.BooTopHLTFilter.HLTPaths = cms.vstring() # do not filter

# source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'file:testPatTuple.root'
    'file:TTJets_madgraph_Summer09.root'
    )
                            )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1) )

    
