# import configurations
import FWCore.ParameterSet.Config as cms

print "About to process"

# define the process
process = cms.Process("TTBSM")

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# note that you can use a bunch of core tools of PAT 
# to taylor your PAT configuration; for a few examples
# uncomment the following lines

from PhysicsTools.PatAlgos.tools.coreTools import *
#removeMCMatching(process, 'Muons')
#removeAllPATObjectsBut(process, ['Muons'])
#removeSpecificPATObjects(process, ['Electrons', 'Muons', 'Taus'])

# add the trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent


print "Setting variables"

outputdir = './'
algorithm = 'ca'
idtag = '_330'

outputFileName = outputdir +  'ttbsm_' + algorithm + '_pat' + idtag + '.root'


print "Output file : " + outputFileName

# CATopJets
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.load("TopQuarkAnalysis.TopPairBSM.caTopJets_cff")
process.load("TopQuarkAnalysis.TopPairBSM.CATopJetTagger_cfi")

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )

# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *

print "About to switch jet collection"

run33xOn31xMC( process,
               jetSrc = cms.InputTag("antikt5CaloJets"),
               jetIdTag = "antikt5"
               )

## ==== Example with CaloJets
addJetCollection(process, 
        cms.InputTag('caTopCaloJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'TopTagCalo',
        doJTA=False,            # Run Jet-Track association & JetCharge
        doBTagging=False,       # Run b-tagging
        jetCorrLabel=('KT6', 'Calo'),   # example jet correction name; set to None for no JEC
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        genJetCollection = cms.InputTag("ca8GenJets"),
        doJetID = False
                 )

## ==== Example withPFJets
addJetCollection(process, 
        cms.InputTag('caTopPFJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
        'TopTagPF',
        doJTA=False,            # Run Jet-Track association & JetCharge
        doBTagging=False,       # Run b-tagging
        jetCorrLabel=('KT6', 'PF'),   # example jet correction name; set to None for no JEC
        doType1MET=False,
        doL1Cleaning=False,
        doL1Counters=False,
        genJetCollection = cms.InputTag("ca8GenJets"),
        doJetID = False
                 )



# Place appropriate jet cuts (NB: no cut on number of constituents)
process.selectedLayer1Jets.cut = cms.string('pt > 30. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagCalo.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagPF.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5')
process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5')
# reduce size of leptons
process.allLayer1Electrons.isoDeposits = cms.PSet()
process.allLayer1Muons.isoDeposits = cms.PSet()

# Jets
from PhysicsTools.PatAlgos.tools.jetTools import *
setTagInfos(process,
            coll = "allLayer1Jets",
            tagInfos = cms.vstring("secondaryVertex")
            )
# Turn off resolutions, they don't mean anything here
process.allLayer1JetsTopTagCalo.addResolutions = cms.bool(False)
# Add CATopTag info... piggy-backing on b-tag functionality
process.allLayer1JetsTopTagCalo.discriminatorSources = cms.VInputTag()
process.allLayer1JetsTopTagCalo.addBTagInfo = cms.bool(True)
process.allLayer1JetsTopTagCalo.addTagInfos = cms.bool(True)
process.allLayer1JetsTopTagCalo.tagInfoSources = cms.VInputTag( cms.InputTag('CATopCaloJetTagInfos') )
process.allLayer1JetsTopTagCalo.addDiscriminators = cms.bool(False)
# Add parton match to quarks and gluons
process.allLayer1JetsTopTagCalo.addGenPartonMatch = cms.bool(True)
process.allLayer1JetsTopTagCalo.embedGenPartonMatch = cms.bool(True)
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
else:
    print "Error: algorithm '%s' unknown. Use one of kt, ca, antikt." % algorithm
    raise AttributeError()

# pythia output
process.printList = cms.EDAnalyzer( "ParticleListDrawer",
                                src = cms.InputTag( "genParticles" ),
                                maxEventsToPrint = cms.untracked.int32( 0 )
)


## produce ttGenEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

# prune gen particles

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *",
    "keep status = 3", #keeps all particles from the hard matrix element
    "+keep (abs(pdgId) = 11 | abs(pdgId) = 13) & status = 1" #keeps all stable muons and electrons and their (direct) mothers.
    )
)


# require >= 1 jets
process.countLayer1Jets.minNumber = cms.uint32(1)



# In addition you usually want to change the following parameters:
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
process.out.fileName = outputFileName
process.options.wantSummary = True       ##  (to suppress the long output at the end of the job)    

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# extend event content to include PAT objects
process.out.outputCommands.extend(['drop *_genParticles_*_*',
                                   'drop *_generalTracks_*_*',
                                   'keep *_prunedGenParticles_*_*',
                                   'keep *_decaySubset_*_*',
                                   'keep *_initSubset_*_*',
                                   'keep *_offlineBeamSpot_*_*',
                                   'keep recoCaloJets_caTopCaloJets_*_*',
                                   'keep recoGenJets_ca8GenJets_*_*',
                                   'keep recoPFJets_caTopPFJets_*_*',
                                   'keep *_CATopCaloJetTagInfos_*_*',
                                   'keep *_CATopPFJetTagInfos_*_*',                                   
                                   "drop *_cleanLayer1Jets*_*_*",
                                   "keep *_selectedLayer1Jets*_*_*",
                                   'drop *_cleanLayer1Taus_*_*',
                                   'drop *_cleanLayer1Hemispheres_*_*',
                                   'drop *_cleanLayer1Photons_*_*',
                                   'keep GenEventInfoProduct_generator_*_*'
                                   #'keep *_CAJetPartonMatcher_*_*',
                                   #'keep *_CAJetFlavourIdentifier_*_*'
                                   ]
                                  )
process.out.outputCommands += patTriggerEventContent

# drop the meta data for dropped data
process.out.dropMetaData = cms.untracked.string("DROPPED")


# define path 'p'
process.p = cms.Path(process.genJetParticles*
                     process.makeGenEvt *
                     process.prunedGenParticles*
                     process.ca8GenJets*
                     process.caTopGenJets*
                     process.caTopCaloJets*
                     process.caTopPFJets*
                     process.CATopCaloJetTagInfos*
                     process.CATopPFJetTagInfos*
                     process.patDefaultSequence*
                     process.countLayer1Jets
                     )

process.source.fileNames = [
    '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3-v1/0025/48AC6C31-AA88-DE11-B02C-0030487C6F54.root',
    '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3-v1/0025/9E80A46A-AA88-DE11-94BD-001E682F882A.root',
    '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3-v1/0025/6004FF4C-AA88-DE11-B0BF-001E68A9941C.root',
    '/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3-v1/0025/129A0B85-AA88-DE11-B08A-001E6837DFEA.root'
    ]
    
process.maxEvents.input = cms.untracked.int32(500)         ##  (e.g. -1 to run on all events)
