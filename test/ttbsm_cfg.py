# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################

# This will apply a dijet skim (2 jets with pt > 25)
skimDijets = True
# Set to true for running on data
useData = False
# Set to true if running on a ttbar sample
useTTHyp = False
# Set to true if running on a single top sample
useSTHyp = False


###############################
####### Global Setup ##########
###############################

if useData == False :
    # global tag for MC
    process.GlobalTag.globaltag = cms.string('START36_V7::All')
else :
    # global tag for 361 data
    process.GlobalTag.globaltag = cms.string('GR_R_35X_V8B::All')

# get the Spring10 jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Spring10")

# require physics declared
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )

# configure HLT
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

# Switch to the REDIGI trigger for MC
if useData == False :
    process.patTrigger.processName = "REDIGI"
    process.patTriggerEvent.processName = "REDIGI"


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )


###############################
########## MC Setup ###########
###############################

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

if useData == False :
    # run ak5 gen jets and b-tagging sequences
    run36xOn35xInput( process, "ak5GenJets")

    # produce ttGenEvent
    if useTTHyp :
        process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    elif useSTHyp :
        process.load("TopQuarkAnalysis.TopEventProducers.sequences.stGenEvent_cff")
        # prune gen particles

else :
    removeMCMatching( process, ['All'] )

###############################
########## PF Setup ###########
###############################

# Add PF Met
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )



###############################
###### Vanilla CA8 Jets #######
###############################


# Gen Jets
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )

from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *

from RecoJets.JetProducers.ca4CaloJets_cfi import ca4CaloJets
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8CaloJets = ca4CaloJets.clone( rParam = cms.double(0.8) )
process.ca8PFJets = ca4PFJets.clone( rParam = cms.double(0.8) )

# calo jet id
from RecoJets.JetProducers.ca4JetID_cfi import ca4JetID
process.ca8JetID = ca4JetID.clone( src = cms.InputTag('ca8CaloJets') )



# Add Calo jets
addJetCollection(process,cms.InputTag('ca8CaloJets'),
                 'CA8', 'Calo',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('KT6','Calo'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ca8GenJets"),
                 doJetID      = True,
                 jetIdLabel = "ca8"                 
                 )



# Add PF jets
addJetCollection(process,cms.InputTag('ca8PFJets'),
                 'CA8', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('KT6','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ca8GenJets"),
                 doJetID      = False                
                 )

###############################
###### Jet Pruning Setup ######
###############################


# Pruned PF Jets
process.caPrunedPFJets = cms.EDProducer(
    "SubJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    SubJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )
process.caPrunedPFJets.nSubjets = cms.int32(2) 
process.caPrunedPFJets.inputEMin = cms.double(1.0)
process.caPrunedPFJets.jetCollInstanceName = cms.string("subjets")


addJetCollection(process, 
                 cms.InputTag('caPrunedPFJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
                 'CA8Pruned', 'PF',
                 doJTA=True,            # Run Jet-Track association & JetCharge
                 doBTagging=True,       # Run b-tagging
                 jetCorrLabel=('KT6', 'PF'),   # example jet correction name; set to None for no JEC
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = False
                 )

# Pruned Calo Jets
process.caPrunedCaloJets = cms.EDProducer(
    "SubJetProducer",
    CaloJetParameters,
    AnomalousCellParameters,
    SubJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )
process.caPrunedCaloJets.nSubjets = cms.int32(2) 
process.caPrunedCaloJets.inputEMin = cms.double(1.0)
process.caPrunedCaloJets.jetCollInstanceName = cms.string("subjets")



addJetCollection(process, 
                 cms.InputTag('caPrunedCaloJets'),
                 'CA8Pruned', 'Calo',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=('KT6', 'Calo'),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = True,
                 jetIdLabel = "ca8"
                 )



###############################
#### CATopTag Setup ###########
###############################

# CATopJet PF Jets
# with adjacency 
process.caTopTagPFJets = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    CATopJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )

process.CATopPFJetTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("caTopTagPFJets"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )


addJetCollection(process, 
                 cms.InputTag('caTopTagPFJets'),
                 'CATopTag', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=('KT6', 'PF'),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = False
                 )





# CATopJet PF Jets
# without adjacency
process.caTopTagPFJetsNoAdj = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    CATopJetParameters.clone( useAdjacency = cms.int32(0) ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )

process.CATopPFJetTagInfosNoAdj = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("caTopTagPFJetsNoAdj"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )


addJetCollection(process, 
                 cms.InputTag('caTopTagPFJetsNoAdj'),
                 'CATopTagNoAdj', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=('KT6', 'PF'),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = False
                 )

# CATopJet Calo Jets
# with adjacency 
process.caTopTagCaloJets = cms.EDProducer(
    "CATopJetProducer",
    CaloJetParameters,
    AnomalousCellParameters,
    CATopJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )

process.CATopCaloJetTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("caTopTagCaloJets"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )


addJetCollection(process, 
                 cms.InputTag('caTopTagCaloJets'),
                 'CATopTag', 'Calo',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=('KT6', 'Calo'),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = True,
                 jetIdLabel = "ca8"
                 )





# CATopJet Calo Jets
# without adjacency
process.caTopTagCaloJetsNoAdj = cms.EDProducer(
    "CATopJetProducer",
    CaloJetParameters,
    AnomalousCellParameters,
    CATopJetParameters.clone( useAdjacency = cms.int32(0) ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )

process.CATopCaloJetTagInfosNoAdj = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("caTopTagCaloJetsNoAdj"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )


addJetCollection(process, 
                 cms.InputTag('caTopTagCaloJetsNoAdj'),
                 'CATopTagNoAdj', 'Calo',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=('KT6', 'Calo'),
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJets"),
                 doJetID = True,
                 jetIdLabel = "ca8"
                 )



###############################
### TagInfo and Matching Setup#
###############################

# Do some configuration of the jet substructure things
for jetcoll in (process.patJets,
                process.patJetsAK5PF,
                process.patJetsCA8Calo,
                process.patJetsCA8PF,
                process.patJetsCA8PrunedCalo,
                process.patJetsCA8PrunedPF,
                process.patJetsCATopTagCalo,
                process.patJetsCATopTagNoAdjCalo,
                process.patJetsCATopTagPF,
                process.patJetsCATopTagNoAdjPF
                ) :
    jetcoll.embedGenJetMatch = cms.bool(False)
    #getJetMCFlavour uses jetFlavourAssociation*, see below
    jetcoll.getJetMCFlavour = cms.bool(True)
    #those two use jetPartonMatch*, see below
    jetcoll.addGenPartonMatch = cms.bool(True)
    # Add CATopTag info... piggy-backing on b-tag functionality
    jetcoll.addBTagInfo = cms.bool(True)
    jetcoll.addTagInfos = cms.bool(True)

#adapt jet parton stuff for topjet collections. There is only one "jetPartons" module ...
process.patJetPartons.withTop = cms.bool(True)
# ... but many "jetPartonAssociation", "jetFlavourAssociation" and "jetPartonMatch":
for jetcollname in ('CATopTagCalo', 'CATopTagPF', 'CATopTagNoAdjCalo', 'CATopTagNoAdjPF' ):
	# definition 4 = match heaviest flavor 
	getattr(process, 'patJetFlavourAssociation'+jetcollname).definition = cms.int32(4)
	getattr(process, 'patJetFlavourAssociation'+jetcollname).physicsDefinition = cms.bool(False)
	# jetFlavourAssociation depends on jetPartonAssociation, so re-configure that as well:
	getattr(process, 'patJetPartonAssociation'+jetcollname).coneSizeToAssociate = cms.double(0.8)
	getattr(process, 'patJetPartonAssociation'+jetcollname).doPriority = cms.bool(True)
	getattr(process, 'patJetPartonAssociation'+jetcollname).priorityList = cms.vint32(6)
	getattr(process, 'patJetPartonMatch'+jetcollname).mcPdgId = cms.vint32(1,2,3,4,5,6,21)
	getattr(process, 'patJetPartonMatch'+jetcollname).maxDeltaR = cms.double(0.8)




###############################
#### Selections Setup #########
###############################


# AK5 Jets
#   Calo
process.selectedPatJets.cut = cms.string("pt > 25")
process.patJets.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfos")
    )
process.patJets.addTagInfos = cms.bool(True)
process.patJets.embedCaloTowers = cms.bool(False)
process.patJets.embedPFCandidates = cms.bool(False)
#   PF
process.selectedPatJetsAK5PF.cut = cms.string("pt > 25")
process.patJetsAK5PF.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfos")
    )
process.patJetsAK5PF.addTagInfos = cms.bool(True)
process.patJetsAK5PF.embedCaloTowers = cms.bool(False)
process.patJetsAK5PF.embedPFCandidates = cms.bool(False)


# CA8 Jets
#   Calo
process.selectedPatJetsCA8Calo.cut = cms.string("pt > 25")
process.patJetsCA8Calo.tagInfoSources = cms.VInputTag()
process.patJetsCA8Calo.embedCaloTowers = cms.bool(False)
process.patJetsCA8Calo.embedPFCandidates = cms.bool(False)
#   PF
process.selectedPatJetsCA8PF.cut = cms.string("pt > 25")
process.patJetsCA8PF.tagInfoSources = cms.VInputTag()
process.patJetsCA8PF.embedCaloTowers = cms.bool(False)
process.patJetsCA8PF.embedPFCandidates = cms.bool(False)


# CA8 Pruned jets
#  Calo
process.selectedPatJetsCA8PrunedCalo.cut = cms.string("pt > 25")
process.patJetsCA8PrunedCalo.tagInfoSources = cms.VInputTag()
#  PF
process.selectedPatJetsCA8PrunedPF.cut = cms.string("pt > 25")
process.patJetsCA8PrunedPF.tagInfoSources = cms.VInputTag()

# CA8 TopJets
#   Calo
process.selectedPatJetsCATopTagCalo.cut = cms.string("pt > 25")
process.patJetsCATopTagCalo.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopCaloJetTagInfos')
    )
#   Calo, no adjacency
process.selectedPatJetsCATopTagNoAdjCalo.cut = cms.string("pt > 25")
process.patJetsCATopTagNoAdjCalo.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopCaloJetTagInfosNoAdj')
    )
#   PF
process.selectedPatJetsCATopTagPF.cut = cms.string("pt > 25")
process.patJetsCATopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopPFJetTagInfos')
    )
#   PF, no adjacency
process.selectedPatJetsCATopTagNoAdjPF.cut = cms.string("pt > 25")
process.patJetsCATopTagNoAdjPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopPFJetTagInfosNoAdj')
    )


# electrons
process.selectedPatElectrons.cut = cms.string('pt > 3.0')
process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.embedTrack = cms.bool(True)
# muons
process.selectedPatMuons.cut = cms.string("pt > 3.0")
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.embedTrack = cms.bool(True)
# taus
process.selectedPatTaus.cut = cms.string("pt > 5 & abs(eta) < 3")
# photons
process.patPhotons.isoDeposits = cms.PSet()
#taus
process.patTaus.isoDeposits = cms.PSet()



# FILTERS:
# For signal samples, no selection
# For high-stats multijet samples, require >= 2 jets

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("selectedPatJets"),
                                  minNumber = cms.uint32(2),
                                  )

process.pfJetFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("selectedPatJetsAK5PF"),
                                    minNumber = cms.uint32(2),
                                    )



# let it run

process.patseq = cms.Sequence(
    process.hltLevel1GTSeed*
    process.scrapingVeto*
    process.hltPhysicsDeclared*
    process.primaryVertexFilter*
    process.genJetParticles*
    process.ca8GenJets*
    process.ca8CaloJets*
    process.ca8JetID*
    process.ca8PFJets*    
    process.caPrunedPFJets*
    process.caPrunedCaloJets*
    process.caTopTagPFJets*
    process.CATopPFJetTagInfos*
    process.caTopTagPFJetsNoAdj*
    process.CATopPFJetTagInfosNoAdj*
    process.caTopTagCaloJets*
    process.CATopCaloJetTagInfos*
    process.caTopTagCaloJetsNoAdj*
    process.CATopCaloJetTagInfosNoAdj*
    process.patDefaultSequence
)

if useData == False :
    process.patseq.remove( process.hltLevel1GTSeed )
    process.patseq.remove( process.hltPhysicsDeclared )

if useTTHyp or useSTHyp :
    process.patseq.append( process.makeGenEvt )

process.p0 = cms.Path(
    process.patseq
    )

process.p1 = cms.Path(
    process.patseq*process.jetFilter
    )

process.p2 = cms.Path(
    process.patseq*process.pfJetFilter
    )

# "or" the preceeding paths
if skimDijets :
    process.out.SelectEvents.SelectEvents = cms.vstring('p1', 'p2')
else :
    process.out.SelectEvents.SelectEvents = cms.vstring('p0')
    
# rename output file
process.out.fileName = cms.untracked.string('ttbsm_361.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = 1000
process.options.wantSummary = True
process.out.dropMetaData = cms.untracked.string("DROPPED")



# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0073/500DAE6F-2B4B-DF11-ABC4-E41F13181564.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/F4CD6579-7346-DF11-837D-E41F13181CF8.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/EAAE740C-7D46-DF11-82C3-00215E222442.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/EA12C324-7B46-DF11-90DF-00215E21DBFA.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/E6DB0D77-7346-DF11-85C1-00215E2211AC.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/E603D5C8-7846-DF11-9903-00215E221FB0.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/DE858FB9-7646-DF11-998C-00215E21DAC8.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/DA54DE19-7F46-DF11-B262-00215E21D9F6.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/D8F67342-7646-DF11-9EBD-00215E22175E.root',
    '/store/mc/Spring10/QCDDiJet_Pt80to120/GEN-SIM-RECO/START3X_V26_S09-v1/0013/D0A9CB03-7646-DF11-82F3-00215E21D86A.root'
        ] );
process.source.fileNames = readFiles


from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
#process.out.outputCommands = patEventContentNoCleaning
#process.out.outputCommands += patExtraAodEventContent
#process.out.outputCommands += patTriggerEventContent

process.out.outputCommands = [
    'drop *_cleanPat*_*_*',
    'keep *_selectedPat*_*_*',
    'keep *_patMETs*_*_*',
    'keep recoPFCandidates_particleFlow_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep recoTracks_generalTracks_*_*',
    'drop patPFParticles_*_*_*',
    'keep recoPFJets_caPrunedPFJets_*_*',
    'keep recoPFJets_caTopTagPFJets_*_*',
    'keep recoPFJets_caTopTagPFJetsNoAdj_*_*',
    'keep recoCaloJets_caPrunedCaloJets_*_*',
    'keep recoCaloJets_caTopTagCaloJets_*_*',
    'keep recoCaloJets_caTopTagCaloJetsNoAdj_*_*',
    'keep patTriggerObjects_patTrigger*_*_*',
    'keep patTriggerFilters_patTrigger*_*_*',
    'keep patTriggerPaths_patTrigger*_*_*',
    'keep patTriggerEvent_patTriggerEvent*_*_*',
    'keep *_cleanPatPhotonsTriggerMatch*_*_*',
    'keep *_cleanPatElectronsTriggerMatch*_*_*',
    'keep *_cleanPatMuonsTriggerMatch*_*_*',
    'keep *_cleanPatTausTriggerMatch*_*_*',
    'keep *_cleanPatJetsTriggerMatch*_*_*',
    'keep *_patMETsTriggerMatch*_*_*'
    ]

if useData :
    process.out.outputCommands += ['keep *_towerMaker_*_*',
                                   'drop *_MEtoEDMConverter_*_*']
else :
    process.out.outputCommands += ['keep *_genParticles_*_*',    
                                   'keep recoGenJets_ca8GenJets_*_*',
                                   'keep recoGenJets_ak5GenJets_*_*',
                                   'keep *_decaySubset_*_*',
                                   'keep *_initSubset_*_*',
                                   'keep *_genEvt*_*_*'
                                   ]
