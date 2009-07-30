import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 production from
# fullsim
#-------------------------------------------------
process = cms.Process("BooTopPat")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0009/9885EAE5-876D-DE11-A0EE-001A9243D640.root',
       '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/DCAE7AE1-2A6D-DE11-9A56-001A9254460C.root',
       '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/D2FE0B42-2F6D-DE11-9FF3-001A9227D383.root',
       '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/C2F0B1DF-2C6D-DE11-8849-001A9243D62A.root',
       '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V2_preproduction_311-v1/0006/9CC08AE1-2A6D-DE11-96C9-001E8CCCE140.root'
    )
    , duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

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


#-------------------------------------------------
# patTuple configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopPairBSM.BooTopPatTuple")


#-------------------------------------------------
# pat tuple event content; first ALL objects
# are dropped in this process; then patTuple
# content is added
#-------------------------------------------------

## define pat tuple event content
from PhysicsTools.PatAlgos.patEventContent_cff import *

##----------------------------------
## Switch jet collection to SC5
##----------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process, 
                    cms.InputTag('sisCone5CaloJets'),             # jet collection; must be already in the event when patLayer0 sequence is executed
                    doJTA        = True,            # run jet-track association & JetCharge
                    doBTagging   = True,            # run b-tagging
                    jetCorrLabel = ('SC5','Calo'), # example jet correction name; set to None for no JEC
                    doType1MET   = False,             # recompute Type1 MET using these jets
		    genJetCollection=cms.InputTag("sisCone5GenJets")
                    )

# if you need to change JEC use the following
# FOR WINTER09 FASTSIM samples comment out the following line
# switchJECSet(process, newName='Winter09', oldName='Summer08Redigi')


#--------------------------------------
# PRE SELECTION
#

process.selectedLayer1Electrons.cut  = cms.string('pt > 15. & abs(eta) < 2.4')
process.selectedLayer1Muons.cut      = cms.string('pt > 15. & abs(eta) < 2.4')
process.selectedLayer1Jets.cut       = cms.string('pt > 20. & abs(eta) < 2.4 & nConstituents > 0')
#process.selectedLayer1METs.cut               = cms.string('et >= 0.') #no met selector at the moment GYJ 2009/07/28
#process.countLayer1Leptons.minNumber = 1
process.countLayer1Jets.minNumber      = 1
process.countLayer1Muons.minNumber     = 1

#-----------------------------------------
# Add JPT collection

process.load("TopQuarkAnalysis.TopPairBSM.ZPT_cff")
process.ZSPJetCorJetIcone5.src = cms.InputTag("sisCone5CaloJets")

addJetCollection(process,
                 cms.InputTag('JetPlusTrackZSPCorJetIcone5'),
                 'JPT',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=None,
                 doType1MET=False,
                 doL1Counters=False,
		 genJetCollection=cms.InputTag("sisCone5GenJets")
		 )

#------------------------------------------------
# Add tcMET collection

# load muon corrections for !MET module
#process.load("JetMETCorrections.Type1MET.MetMuonCorrections_cff")

# load track-corrected MET module
#process.load("RecoMET.METProducers.TCMET_cfi")

#def addAlso (label,value):
#        existing = getattr(process, label)
#        setattr( process, label + "tcMET", value)
#        process.patLayer0.replace( existing, existing * value )
#        process.patLayer1.replace( existing, existing * value )
#
#def addClone(label,**replaceStatements):
#        new      = getattr(process, label).clone(**replaceStatements)
#        addAlso(label, new)

##addClone('metJESCorIC5CaloJetMuons', metSource = cms.InputTag('tcMet')) # this does not work need to setup
## by hand
#process.metJESCorIC5CaloJetMuonstcMET = cms.EDFilter("PATBaseMETCleaner",
    ## Input MET from AOD
#    metSource = cms.InputTag('tcMet'), ## met corrected for jets and for muons
    #metSource = cms.InputTag('met'),                     ## NO MET CORRECTIONS

#    markItems = cms.bool(True),    ## write the status flags in the output items
#    bitsToIgnore = cms.vstring(),  ## You can specify some bit names, e.g. "Overflow/User1", "Core/Duplicate", "Isolation/All".
#    saveRejected = cms.string(''), ## set this to a non empty label to save the list of items which fail
#    saveAll = cms.string(''),      ## set this to a non empty label to save a list of all items both passing and failing
#)

#l0met = getattr(process, 'metJESCorIC5CaloJetMuons'+'tcMET')
#addAlso('metJESCorIC5CaloJetMuons', l0met)

#addClone('layer1METs', metSource = cms.InputTag('metJESCorIC5CaloJetMuons'+'tcMET'))
#l1MET = getattr(process, 'allLayer1METs'+'tcMET')
#addClone('layer1METs', src=cms.InputTag('layer1METs'+'tcMET'))


#--------------------------------------------------
# PATH

process.p = cms.Path(process.recoJPTJets+
#                     process.MetMuonCorrections*process.tcMet+
                     process.BooTopPatTuple)


#-------------------------------------------------
# process output; first the event selection is
# defined: only those events that have passed the
# full production path are selected and written
# to file; the event content has been defined
# above
#-------------------------------------------------

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    verbose = cms.untracked.bool(True),
    outputCommands = cms.untracked.vstring('drop *'),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                           
    fileName = cms.untracked.string('TTJets_madgraph_Summer09.root'),
    dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('USER'),
            filterName = cms.untracked.string('')
                )
)
process.out.outputCommands += patEventContent
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patEventContentNoLayer1Cleaning

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1.2.7.2.1 $'),
    annotation = cms.untracked.string('PAT tuple creation'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/TopQuarkAnalysis/TopPairBSM/test/Attic/BooTopPatTuple_cfg.py,v $')
)



#-------------------------------------------------
# output paths; in order not to write the
# persistent output to file comment the output
# path
#-------------------------------------------------

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.px=cms.Path(process.dump)

## output
process.outpath = cms.EndPath(process.out)
