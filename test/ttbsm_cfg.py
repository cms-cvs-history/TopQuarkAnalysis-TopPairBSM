# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('skimDijets',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Apply a skim of >= 2jets")

options.register ('useData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run this on real data")

options.register ('useTTHyp',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "On MC, run the ttbar event hypothesis")

options.register ('useSTHyp',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "On MC, run the single top event hypothesis")

options.register ('runOn35x',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Run on samples produced with <= 35x")
                  
options.parseArguments()

mytrigs = ['HLT_Jet30U']

print options

import sys

# This will apply a dijet skim (2 jets with pt > 25)
skimDijets = options.skimDijets
# Set to true for running on data
useData = options.useData
# Set to true if running on a ttbar sample
useTTHyp = options.useTTHyp
# Set to true if running on a single top sample
useSTHyp = options.useSTHyp
# Set to true to run on < 36x samples
runOn35x = options.runOn35x

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
if useData == False and runOn35x:
    process.patTrigger.processName = "REDIGI"
    process.patTriggerEvent.processName = "REDIGI"

if useData == True and skimDijets == True :
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    process.hltSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = mytrigs)

    



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

if runOn35x == True:
    if useData == False :
        # run ak5 gen jets and b-tagging sequences
        run36xOn35xInput( process, "ak5GenJets")
    else :
        run36xOn35xInput( process )    

if useData == False:
    # produce ttGenEvent
    if useTTHyp :
        process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    elif useSTHyp :
        process.load("TopQuarkAnalysis.TopEventProducers.sequences.stGenEvent_cff")
        # prune gen particles


###############################
########## PF Setup ###########
###############################

from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
if useData == False :
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
else :
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix)

# turn to false when running on data
if useData :
    getattr(process, "patPhotons"+postfix).embedGenMatch = False
    getattr(process, "patElectrons"+postfix).embedGenMatch = False
    getattr(process, "patMuons"+postfix).embedGenMatch = False

    removeMCMatching( process, ['All'] )

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


process.PF2PAT += ( process.ca8PFJets*
                    process.caPrunedPFJets*
                    process.caTopTagPFJets*
                    process.CATopPFJetTagInfos)

# add the modules to the sequence

for module in (
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


###############################
### TagInfo and Matching Setup#
###############################

# Do some configuration of the jet substructure things
for jetcoll in (process.patJets,
                process.patJetsPFlow,
                process.patJetsCA8PF,
                process.patJetsCA8PrunedPF,
                process.patJetsCATopTagPF
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
#for jetcollname in ('CATopTagPF') :
jetcollname = 'CATopTagPF'
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
                                   cut = cms.string('pt > 10')
                                   )

# CA8 Jets
#   PF
process.selectedPatJetsCA8PF.cut = cms.string("pt > 30")
process.patJetsCA8PF.tagInfoSources = cms.VInputTag()
process.selectedPatJetsCA8PFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCA8PF'),
                                   cut = cms.string('pt > 30')
                                   )

# CA8 Pruned jets
#  PF
process.selectedPatJetsCA8PrunedPF.cut = cms.string("pt > 30")
process.patJetsCA8PrunedPF.tagInfoSources = cms.VInputTag()
process.selectedPatJetsCA8PrunedPFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCA8PrunedPF'),
                                   cut = cms.string('pt > 30')
                                   )

# CA8 TopJets
#   PF
process.selectedPatJetsCATopTagPF.cut = cms.string("pt > 30")
process.patJetsCATopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopPFJetTagInfos')
    )
process.selectedPatJetsCATopTagPFClones = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag('selectedPatJetsCATopTagPF'),
                                   cut = cms.string('pt > 30')
                                   )

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




# remove trigger matching for PF2PAT as that is currently broken
process.patPF2PATSequencePFlow.remove(process.patTriggerSequencePFlow)
    
# let it run

process.patseq = cms.Sequence(
    process.hltSelection*
    process.scrapingVeto*
    process.primaryVertexFilter*
    process.genJetParticles*
    process.ca8GenJets*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.patDefaultSequence*
    process.selectedPatJetsPFlowClones*
    process.selectedPatJetsCA8PFClones*
    process.selectedPatJetsCA8PrunedPFClones*
    process.selectedPatJetsCATopTagPFClones
)

if useData == True :
    process.patseq.remove( process.genJetParticles )
    process.patseq.remove( process.ca8GenJets )
else :
    process.patseq.remove( process.hltSelection )

if (useTTHyp or useSTHyp) and useData == False :
    process.patseq += process.makeGenEvt


process.p0 = cms.Path(
    process.patseq
    )

process.out.SelectEvents.SelectEvents = cms.vstring('p0')
    
# rename output file
process.out.fileName = cms.untracked.string('ttbsm_381.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.MessageLogger.suppressWarning = cms.untracked.vstring("patTrigger");


# process all the events
process.maxEvents.input = 100
process.options.wantSummary = True
process.out.dropMetaData = cms.untracked.string("DROPPED")



# Add the files 
readFiles38x = cms.untracked.vstring()
readFiles36x = cms.untracked.vstring()
secFiles = cms.untracked.vstring()



readFiles38x.extend( [
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0160/B249DC90-2798-DF11-96D3-001A92811728.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FACE6E6F-B297-DF11-8C03-002618FDA21D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FABF5E30-E897-DF11-8249-002618943836.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/FA906A3C-AA97-DF11-B010-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/F650EB46-D697-DF11-A374-00261894396E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/ECE61B70-B497-DF11-953A-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/EA5216EA-AD97-DF11-B7D2-002618943860.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E86D8349-A597-DF11-92FE-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E868BC6C-E497-DF11-B07C-002618943967.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E6D0EE42-A897-DF11-AB1E-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E2F27541-AA97-DF11-AFD3-00304867BFB0.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E2EBD86E-B297-DF11-9C98-0026189438D4.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/E081EC45-E397-DF11-86D9-00261894396E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/DAB2CB4A-B197-DF11-AAB6-003048678B14.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/D869396E-A797-DF11-802D-0026189438B5.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/D682A16A-A997-DF11-9489-0026189438B5.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/D4A0101E-E397-DF11-AED9-002618943967.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/D2522E3F-E597-DF11-82A1-00261894393A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/D0822148-BB97-DF11-854E-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/CCE8FCB2-A497-DF11-9B16-00248C55CC40.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/CC81C150-E097-DF11-B71A-00261894393A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/CC30DDBF-D397-DF11-94CB-00261894396E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C89C8342-AC97-DF11-AF82-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C8011413-AF97-DF11-BDAD-002618FDA21D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C6FDD670-BC97-DF11-A97C-002618943964.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C27D3737-E797-DF11-BB5A-00261894393A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C2390E14-AF97-DF11-9260-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/C0312FB1-A597-DF11-896C-0026189438D4.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/BEEE5E52-E297-DF11-B51D-002618943926.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/B040AB56-A897-DF11-9DAD-0030486792DE.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/AED0FC67-AD97-DF11-B124-002618FDA21D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/AC7C5944-AC97-DF11-A45C-003048678A78.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/AA737E0B-DA97-DF11-9A6B-0026189438A0.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/A2412344-B497-DF11-930E-002618FDA21D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/9E0A1470-B997-DF11-B673-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/9C71293C-DC97-DF11-AE96-002618FDA287.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/9A69592A-E197-DF11-8A70-002618FDA287.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/98DBC699-A897-DF11-AF4A-003048678B08.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/98254051-E297-DF11-BB55-002618FDA208.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/967C255F-DB97-DF11-8FDF-002618FDA262.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/940B070B-DA97-DF11-8617-002618943926.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/9281D65A-DD97-DF11-A291-002618FDA265.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/88DA882B-DE97-DF11-BC66-00248C55CC97.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/888AA340-AE97-DF11-A4F6-003048678A78.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/86341F45-B897-DF11-B123-003048678B14.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/846FB748-B697-DF11-9424-002618FDA216.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/74E90F3D-AA97-DF11-B477-0026189438D4.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/74CBB955-DD97-DF11-8A29-00261894393A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/74422470-B997-DF11-B5B6-002618943964.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/7435EB3D-AA97-DF11-9542-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/72DDC354-DD97-DF11-901F-002618FDA211.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/72686F6C-E497-DF11-9426-002618943836.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/6A237350-E097-DF11-A811-00261894396E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/688CC9A8-D497-DF11-A82A-002618943926.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/66BB531E-D597-DF11-9B8A-0026189438BA.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/62AC2C6B-A997-DF11-98A4-002618FDA210.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/6299FC2A-E197-DF11-8773-002618FDA265.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/6091AA48-B197-DF11-9BCF-00261894396A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/5CC39069-B397-DF11-893F-003048678A78.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/58D9079B-A897-DF11-8B6B-003048678A78.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/563D9DA8-AD97-DF11-8CC7-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/54127B46-B497-DF11-8758-00261894393C.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/52445FC5-AD97-DF11-B021-002618FDA216.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/520B2B2F-DE97-DF11-80B3-002618FDA287.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4CB31A7F-DE97-DF11-B046-002618FDA262.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4C589749-B697-DF11-8059-00261894396A.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4AB60F1F-D597-DF11-9AC0-00248C55CC97.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4A3AC46E-A797-DF11-BB72-003048678B5E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4809E146-B397-DF11-94B3-003048678B14.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4217C067-AB97-DF11-B4B0-002618FDA21D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/4035A012-AF97-DF11-91C1-003048678B5E.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/3A8F3B49-AC97-DF11-9E4B-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/3A75E46D-DA97-DF11-B845-002618943836.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/367EDB45-B197-DF11-84A9-002618FDA216.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/32571795-A697-DF11-8320-00248C55CC40.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/2C935C68-AD97-DF11-9519-0026189438D4.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/264BE66D-B297-DF11-82B3-002618943860.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/261BF88D-E297-DF11-9F26-002618FDA287.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/247AEBA5-9D97-DF11-B5F4-00261894385D.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/22744443-AF97-DF11-8450-002618943964.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1E62C747-B697-DF11-A581-002618943949.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1CE7213F-D797-DF11-A0C9-00248C55CC97.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1A368544-B897-DF11-A6DB-00261894393C.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1A11F06D-B797-DF11-857F-002618943964.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1898FD6F-AF97-DF11-AF4D-002354EF3BE1.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/162AD99E-E297-DF11-87E2-002618943985.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/161AF945-B397-DF11-9FAE-002618FDA216.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/14DB0EA6-A697-DF11-B754-0026189438D4.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/14507168-AD97-DF11-99B4-0030486792DE.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/1273FA44-B397-DF11-BFB6-002618943946.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/10B71087-AF97-DF11-83CC-003048678A78.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/10A5BF6E-B797-DF11-84A5-002618943915.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/0C3A4055-E297-DF11-8321-00248C55CC97.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/06C5A183-F097-DF11-885C-002618FDA211.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/00E19D5E-DB97-DF11-A4D3-002618943926.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/009B303F-E597-DF11-9FC6-002618943811.root',
'/store/data/Run2010A/JetMETTau/RECO/Jul23ReReco_PreProd_v1/0157/0003F63E-E597-DF11-94EB-002618FDA208.root'

        ] );




readFiles36x.extend( [
'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTau-Jun14thSkim_v1/0149/5000F2BB-7D84-DF11-9BB8-00261894386C.root',
'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTau-Jun14thSkim_v1/0129/54A7C1C6-4280-DF11-8388-0018F3D096AE.root',
'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTau-Jun14thSkim_v1/0107/34957F4B-437E-DF11-94B8-00261894387C.root'

        ] );


process.source.fileNames = readFiles36x

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
    'keep recoCaloJets_caPrunedCaloJets_*_*',
    'keep recoCaloJets_caTopTagCaloJets_*_*',
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
