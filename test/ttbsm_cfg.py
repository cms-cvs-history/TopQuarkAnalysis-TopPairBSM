# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################

# This will apply a dijet skim (2 jets with pt > 25)
skimDijets = True
# Set to true for running on data
useData = True
# Set to true if running on a ttbar sample
useTTHyp = False
# Set to true if running on a single top sample
useSTHyp = False


###############################
####### Global Setup ##########
###############################

if useData == False :
    # global tag for MC
    process.GlobalTag.globaltag = cms.string('START38_V8::All')
else :
    # global tag for 361 data
    process.GlobalTag.globaltag = cms.string('GR_R_38X_V7::All')

# get the Spring10 jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Spring10")


# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )


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
    # run b-tagging sequences
    # run36xOn35xInput( process )    
    removeMCMatching( process, ['All'] )

###############################
########## PF Setup ###########
###############################

process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")


from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
if useData == False :
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix, explicitExclude=True)
else :
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix, explicitExclude=True)

# turn to false when running on data
if useData :
    getattr(process, "patElectrons"+postfix).embedGenMatch = False
    getattr(process, "patMuons"+postfix).embedGenMatch = False




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
process.ca8PFJets = ca4PFJets.clone( rParam = cms.double(0.8),
                                     src = cms.InputTag('pfNoElectron'+postfix) )


###############################
###### Jet Pruning Setup ######
###############################


# Pruned PF Jets
process.caPrunedPFJets = cms.EDProducer(
    "SubJetProducer",
    PFJetParameters.clone( src = cms.InputTag('pfNoElectron'+postfix) ),
    AnomalousCellParameters,
    SubJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = SubJetParameters.jetSize
    )
process.caPrunedPFJets.nSubjets = cms.int32(2) 
process.caPrunedPFJets.inputEMin = cms.double(1.0)
process.caPrunedPFJets.jetCollInstanceName = cms.string("subjets")



###############################
#### CATopTag Setup ###########
###############################

# CATopJet PF Jets
# with adjacency 
process.caTopTagPFJets = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag('pfNoElectron'+postfix) ),
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






# CATopJet PF Jets
# without adjacency
process.caTopTagPFJetsNoAdj = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag('pfNoElectron'+postfix) ),
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


process.PF2PAT += ( process.caPrunedPFJets*
                    process.caTopTagPFJets*
                    process.CATopPFJetTagInfos*
                    process.caTopTagPFJetsNoAdj*
                    process.CATopPFJetTagInfosNoAdj)

# add the modules to the sequence

for module in (
    process.CATopPFJetTagInfosNoAdj,
    process.caTopTagPFJetsNoAdj,
    process.CATopPFJetTagInfos,
    process.caTopTagPFJets,
    process.caPrunedPFJets,
    process.ca8PFJets  
    ) :
    getattr(process,"patPF2PATSequence"+postfix).replace( getattr(process,"pfNoElectron"+postfix), getattr(process,"pfNoElectron"+postfix)*module )

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



###############################
### TagInfo and Matching Setup#
###############################

# Do some configuration of the jet substructure things
for jetcoll in (process.patJets,
                process.patJetsPFlow,
                process.patJetsCA8PF,
                process.patJetsCA8PrunedPF,
                process.patJetsCATopTagPF,
                process.patJetsCATopTagNoAdjPF
                ) :
    if useData == False :
        jetcoll.embedGenJetMatch = cms.bool(True)
        jetcoll.getJetMCFlavour = cms.bool(True)
        jetcoll.addGenPartonMatch = cms.bool(True)
    # Add CATopTag info... piggy-backing on b-tag functionality
    jetcoll.addBTagInfo = cms.bool(True)
    jetcoll.addTagInfos = cms.bool(True)
    jetcoll.embedCaloTowers = cms.bool(False)
    jetcoll.embedPFCandidates = cms.bool(False)
    jetcoll.embedGenJetMatch = cms.bool(False)

#adapt jet parton stuff for topjet collections. There is only one "jetPartons" module ...
process.patJetPartons.withTop = cms.bool(True)
# ... but many "jetPartonAssociation", "jetFlavourAssociation" and "jetPartonMatch":
for jetcollname in ('CATopTagPF', 'CATopTagNoAdjPF' ):
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
#   PF
process.selectedPatJetsPFlow.cut = cms.string("pt > 10")
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfos")
    )
process.patJetsPFlow.addTagInfos = cms.bool(True)
process.patJetsPFlow.embedCaloTowers = cms.bool(False)
process.patJetsPFlow.embedPFCandidates = cms.bool(False)
process.selectedPatJetsPFlowClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsPFlow'),
                                   cut = cms.string('pt > 10 & abs(eta) < 5')
                                   )

# CA8 Jets
#   PF
process.selectedPatJetsCA8PF.cut = cms.string("pt > 10")
process.patJetsCA8PF.tagInfoSources = cms.VInputTag()
process.selectedPatJetsCA8PFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCA8PF'),
                                   cut = cms.string('pt > 10 & abs(eta) < 5')
                                   )

# CA8 Pruned jets
#  PF
process.selectedPatJetsCA8PrunedPF.cut = cms.string("pt > 10")
process.patJetsCA8PrunedPF.tagInfoSources = cms.VInputTag()
process.selectedPatJetsCA8PrunedPFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCA8PrunedPF'),
                                   cut = cms.string('pt > 10 & abs(eta) < 5')
                                   )

# CA8 TopJets
#   PF
process.selectedPatJetsCATopTagPF.cut = cms.string("pt > 10")
process.patJetsCATopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopPFJetTagInfos')
    )
process.selectedPatJetsCATopTagPFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCATopTagPF'),
                                   cut = cms.string('pt > 10 & abs(eta) < 5')
                                   )
#   PF, no adjacency
process.selectedPatJetsCATopTagNoAdjPF.cut = cms.string("pt > 10")
process.patJetsCATopTagNoAdjPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopPFJetTagInfosNoAdj')
    )
process.selectedPatJetsCATopTagNoAdjPFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCATopTagNoAdjPF'),
                                   cut = cms.string('pt > 10 & abs(eta) < 5')
                                   )

# electrons
process.selectedPatElectrons.cut = cms.string('pt > 3.0')
process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.embedTrack = cms.bool(True)
process.patElectrons.usePV = cms.bool(False)
# muons
process.selectedPatMuons.cut = cms.string("pt > 3.0")
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.embedTrack = cms.bool(True)
process.patMuons.usePV = cms.bool(False)
# taus
process.selectedPatTaus.cut = cms.string("pt > 5 & abs(eta) < 3")
# photons
process.patPhotons.isoDeposits = cms.PSet()
#taus
process.patTaus.isoDeposits = cms.PSet()



# electrons
process.selectedPatElectronsPFlow.cut = cms.string('pt > 3.0')
process.patElectronsPFlow.isoDeposits = cms.PSet()
process.patElectronsPFlow.embedTrack = cms.bool(True)
process.patElectronsPFlow.usePV = cms.bool(False)
# muons
process.selectedPatMuonsPFlow.cut = cms.string("pt > 3.0")
process.patMuonsPFlow.isoDeposits = cms.PSet()
process.patMuonsPFlow.embedTrack = cms.bool(True)
process.patMuonsPFlow.usePV = cms.bool(False)
# taus
process.selectedPatTausPFlow.cut = cms.string("pt > 5 & abs(eta) < 3")
# photons
process.patPhotonsPFlow.isoDeposits = cms.PSet()
#taus
process.patTausPFlow.isoDeposits = cms.PSet()


# FILTERS:
# For signal samples, no selection
# For high-stats multijet samples, require >= 2 jets

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("selectedPatJets"),
                                  minNumber = cms.uint32(2),
                                  )

process.pfJetFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("selectedPatJetsPFlow"),
                                    minNumber = cms.uint32(2),
                                    )



# remove trigger matching for PF2PAT as that is currently broken
process.patPF2PATSequencePFlow.remove(process.patTriggerSequencePFlow)
    
# let it run

process.patseq = cms.Sequence(
    process.scrapingVeto*
    process.primaryVertexFilter*
    process.genJetParticles*
    process.ca8GenJets*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.patDefaultSequence*
    process.selectedPatJetsPFlowClones*
    process.selectedPatJetsCA8PFClones*
    process.selectedPatJetsCA8PrunedPFClones*
    process.selectedPatJetsCATopTagPFClones*
    process.selectedPatJetsCATopTagNoAdjPFClones
)

if useData == True :
    process.patseq.remove( process.genJetParticles )
    process.patseq.remove( process.ca8GenJets )

if (useTTHyp or useSTHyp) and useData == False :
    process.patseq += process.makeGenEvt

if (skimDijets == False or useData == False) :
    process.patseq.remove( process.JetMETTau_1e28 )

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
process.out.fileName = cms.untracked.string('ttbsm_381.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# process all the events
process.maxEvents.input = 100
process.options.wantSummary = True
process.out.dropMetaData = cms.untracked.string("DROPPED")



# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FEFBA245-1D6A-DF11-A909-002618943901.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FCA432B0-C769-DF11-9402-002618943963.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FC75616A-C769-DF11-884C-0026189437FE.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FC309771-D969-DF11-B2BE-00304867C1B0.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FC14372B-196A-DF11-98FE-003048678E24.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FAB5C1DB-1B6A-DF11-99FA-003048679046.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FA90C3AB-EC69-DF11-A8D7-003048679296.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/FA81B0C2-1F6A-DF11-920A-002354EF3BDB.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F8E89DD0-0F6A-DF11-8125-001A92971B62.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F8565F96-106A-DF11-9A6C-002618943916.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F809CFF5-D469-DF11-98D1-00261894395B.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F6CCAD64-B769-DF11-B457-0018F3D09680.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F6492DA0-9D69-DF11-8CD6-001A9281173A.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F61D9497-CA69-DF11-9E91-002618943927.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F60273F3-FA69-DF11-A9EB-003048678E94.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F48E25CF-D869-DF11-8521-0026189438F2.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F47619B6-296A-DF11-88FA-003048D3C010.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F415E3F4-9B69-DF11-9156-001BFCDBD1BC.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F2DE2005-9B69-DF11-92C5-003048678F92.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F211DFED-CA69-DF11-A6CB-00304867C034.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F08985EC-276A-DF11-B408-0030486792B8.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F071FBF9-DF69-DF11-8AAB-00261894383A.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F03DCDD5-9D69-DF11-B6B8-001A92811738.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/F016CC3D-9F69-DF11-8F54-0030486792BA.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/EE0752E5-DC69-DF11-AE82-001A92971B36.root',
'/store/data/Commissioning10/MinimumBias/RECO/May27thReReco_v1/0017/EC7FFD93-9969-DF11-A8A7-003048678E92.root'
        ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")


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
