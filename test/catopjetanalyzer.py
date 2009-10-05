
# import configurations
import FWCore.ParameterSet.Config as cms

print "About to process"

# define the process
process = cms.Process("TTBSM")

# This is an example PAT configuration showing the usage of PAT on full sim samples

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# note that you can use a bunch of core tools of PAT 
# to taylor your PAT configuration; for a few examples
# uncomment the following lines

from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, 'Muons')
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Electrons', 'Muons', 'Taus'])

print "Setting variables"

outputdir = './'
algorithm = 'ca'
idtag = '_330'

outputFileName = outputdir +  'ttbsm_' + algorithm + '_pat' + idtag + '.root'


print "Output file : " + outputFileName

# CATopJets
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.cambridge6GenJets_cff import cambridge6GenJets
process.load("TopQuarkAnalysis.TopPairBSM.caTopJets_cff")
process.load("TopQuarkAnalysis.TopPairBSM.CATopJetTagger_cfi")

process.cambridge8GenJets = cambridge6GenJets.clone( FJ_ktRParam = cms.double(0.8) )

# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *

print "About to switch jet collection"

switchJetCollection(process, 
        'antikt5CaloJets',     # Jet collection; must be already in the event when patLayer0 sequence is executed
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('AK5', 'Calo'),   # example jet correction name; set to None for no JEC
        doType1MET=False,
        genJetCollection = cms.InputTag("antikt5GenJets")
                    )

## ==== Example with CaloJets
addJetCollection(process, 
        'caTopCaloJets',         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'TopTagCalo',
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('KT6', 'Calo'),   # example jet correction name; set to None for no JEC
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        genJetCollection = cms.InputTag("cambridge8GenJets")
                 )

## ==== Example withPFJets
addJetCollection(process, 
        'caTopPFJets',         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'TopTagPF',
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel=('KT6', 'PF'),   # example jet correction name; set to None for no JEC
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        genJetCollection = cms.InputTag("cambridge8GenJets")
                 )

# Place appropriate jet cuts (NB: no cut on number of constituents)
process.selectedLayer1Jets.cut = cms.string('pt > 30. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagCalo.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagPF.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & muonID("TMLastStationLoose")')
process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5 & electronID("eidLoose")')
# reduce size of leptons
process.allLayer1Electrons.isoDeposits = cms.PSet()
process.allLayer1Electrons.embedGsfTrack = cms.bool(False)
process.allLayer1Electrons.embedSuperCluster = cms.bool(False)
process.allLayer1Electrons.embedPFCandidate = cms.bool(False)
process.allLayer1Muons.isoDeposits = cms.PSet()
process.allLayer1Muons.addTeVRefits = cms.bool(False)
process.allLayer1Muons.embedCombinedMuon = cms.bool(False)
process.allLayer1Muons.embedStandAloneMuon = cms.bool(False)
process.allLayer1Muons.embedPickyMuon = cms.bool(False)
process.allLayer1Muons.embedTpfmsMuon = cms.bool(False)


# Jets

# Turn off resolutions, they don't mean anything here
process.allLayer1JetsTopTagCalo.addResolutions = cms.bool(False)
# Add CATopTag info... piggy-backing on b-tag functionality
process.allLayer1JetsTopTagCalo.discriminatorSources = cms.VInputTag()
process.allLayer1JetsTopTagCalo.addBTagInfo = cms.bool(True)
process.allLayer1JetsTopTagCalo.addTagInfos = cms.bool(True)
process.allLayer1JetsTopTagCalo.tagInfoSources = cms.VInputTag( cms.InputTag('CATopCaloJetTagInfos') )
process.allLayer1JetsTopTagCalo.addDiscriminators = cms.bool(False)
# Add parton match to quarks and gluons
process.allLayer1JetsTopTagCalo.addGenPartonMatch = cms.bool(False)
process.allLayer1JetsTopTagCalo.embedGenPartonMatch = cms.bool(False)
process.allLayer1JetsTopTagCalo.embedGenJetMatch = cms.bool(False)
# Add jet MC flavour (custom built to capture tops)
process.allLayer1JetsTopTagCalo.getJetMCFlavour = cms.bool(True)

# Turn off resolutions, they don't mean anything here
process.allLayer1JetsTopTagPF.addResolutions = cms.bool(False)
# Add CATopTag info... piggy-backing on b-tag functionality
process.allLayer1JetsTopTagPF.discriminatorSources = cms.VInputTag()
process.allLayer1JetsTopTagPF.addBTagInfo = cms.bool(True)
process.allLayer1JetsTopTagPF.addTagInfos = cms.bool(True)
process.allLayer1JetsTopTagPF.tagInfoSources = cms.VInputTag( cms.InputTag('CATopPFJetTagInfos') )
process.allLayer1JetsTopTagPF.addDiscriminators = cms.bool(False)
# Add parton match to quarks and gluons
process.allLayer1JetsTopTagPF.addGenPartonMatch = cms.bool(False)
process.allLayer1JetsTopTagPF.embedGenPartonMatch = cms.bool(False)
process.allLayer1JetsTopTagPF.embedGenJetMatch = cms.bool(False)
# Add jet MC flavour (custom built to capture tops)
process.allLayer1JetsTopTagPF.getJetMCFlavour = cms.bool(True)

# jet flavor stuff
process.jetPartons.withTop = cms.bool(True)
process.jetPartonAssociation.coneSizeToAssociate = cms.double(0.8)
process.jetPartonAssociation.doPriority = cms.bool(True)
process.jetPartonAssociation.priorityList = cms.vint32(6)
process.jetFlavourAssociation.definition = cms.int32(4)
process.jetFlavourAssociation.physicsDefinition = cms.bool(False)
process.jetPartonMatch.mcPdgId = cms.vint32(1,2,3,4,5,6,21)
process.jetPartonMatch.maxDeltaR = cms.double(0.8)


#process.allLayer1Jets.JetPartonMapSource = cms.InputTag("CAJetFlavourIdentifier")

print "Done switching jet collection"

# only keep events that have at least one jet
process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("cleanLayer1Jets"),
                                  minNumber = cms.uint32( 1 )
                                  )

if algorithm == 'kt' :
    process.caTopCaloJets.algorithm = cms.int32(0)
    process.caTopGenJets.algorithm = cms.int32(0)
    process.caTopPFJets.algorithm = cms.int32(0)
elif algorithm == 'ca' :
    process.caTopCaloJets.algorithm = cms.int32(1)
    process.caTopGenJets.algorithm = cms.int32(1)
    process.caTopPFJets.algorithm = cms.int32(1)    
elif algorithm == 'antikt' :
    process.caTopCaloJets.algorithm = cms.int32(2)
    process.caTopGenJets.algorithm = cms.int32(2)
    process.caTopPFJets.algorithm = cms.int32(2)

# pythia output
process.printList = cms.EDAnalyzer( "ParticleListDrawer",
                                src = cms.InputTag( "genParticles" ),
                                maxEventsToPrint = cms.untracked.int32( 0 )
)

# In addition you usually want to change the following parameters:
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#   process.source.fileNames = [
#     '/store/relval/CMSSW_3_1_1/RelValCosmics/GEN-SIM-RECO/STARTUP31X_V1-v2/0002/7625DA7D-E36B-DE11-865A-000423D174FE.root'
#                               ]         ##  (e.g. 'file:AOD.root')
process.maxEvents.input = cms.untracked.int32(100)
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
process.out.fileName = outputFileName
process.options.wantSummary = True       ##  (to suppress the long output at the end of the job)    

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# extend event content to include PAT objects
process.out.outputCommands.extend(['drop *_genParticles_*_*',
                                   'drop *_generalTracks_*_*',
                                   'keep recoCaloJets_caTopCaloJets_*_*',
                                   'keep recoGenJets_ca8GenJets_*_*',
                                   'keep recoPFJets_caTopPFJets_*_*',
                                   'keep *_CATopCaloJetTagInfos_*_*',
                                   'keep *_CATopPFJetTagInfos_*_*',                                   
                                   "drop *_cleanLayer1Jets*_*_*",
                                   "keep *_selectedLayer1Jets*_*_*",
                                   'drop *_cleanLayer1Taus_*_*',
                                   'drop *_cleanLayer1Hemispheres_*_*',
                                   'drop *_cleanLayer1Photons_*_*'
                                   #'keep *_CAJetPartonMatcher_*_*',
                                   #'keep *_CAJetFlavourIdentifier_*_*'
                                   ]
                                  )

# drop the meta data for dropped data
process.out.dropMetaData = cms.untracked.string("DROPPED")


# define path 'p'
process.p = cms.Path(process.genJetParticles*
                     process.cambridge8GenJets*
                     process.caTopGenJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CATopCaloJetTagInfos*
                     process.CATopPFJetTagInfos*
                     process.patDefaultSequence*
                     process.jetFilter
                     )


# print process.dumpPython()
