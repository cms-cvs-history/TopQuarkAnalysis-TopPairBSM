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
#set 'algorithm' to 'kt', 'antiki' or 'ca':
algorithm = 'ca'
# this is an ID tag to specify versions, etc
idtag = '_332'
# this is the name of the output file.
outputFileName = outputdir +  'ttbsm_' + algorithm + '_pat' + idtag + '.root'
print "Output file : " + outputFileName
#set 'runon' to '31x' if you intent to run on data which was reconstructed with CMSSW 31X, and to '33x' if you want
# to run it on 33x (fastsim). -- Jochen
runon = '332rereco'
#set 'type' to 'ttbar' if you want to run the ttbar gen event
type = 'qcd'

# CATopJets
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.load("TopQuarkAnalysis.TopPairBSM.caTopJets_cff")
process.load("TopQuarkAnalysis.TopPairBSM.CATopJetTagger_cfi")

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )

# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *

print "About to switch jet collection"

if runon=='31x':
    run33xOn31xMC( process,
               jetSrc = cms.InputTag("antikt5CaloJets"),
               jetIdTag = "antikt5")

if runon=='332rereco':
    run33xOnReRecoMC( process, "ak5GenJets" )

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

#jet ID does not work reliably in PAT for now (PAT as of 2009-10-20), so switch it off -- Jochen:
process.allLayer1Jets.addJetID = cms.bool(False)
# Place appropriate jet cuts (NB: no cut on number of constituents)
process.selectedLayer1Jets.cut = cms.string('pt > 20. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagCalo.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1JetsTopTagPF.cut = cms.string('pt > 250. & abs(eta) < 5.0')
process.selectedLayer1Muons.cut = cms.string('pt > 20. & abs(eta) < 2.5')
process.selectedLayer1Electrons.cut = cms.string('pt > 20. & abs(eta) < 2.5')
# reduce size of leptons
process.allLayer1Electrons.isoDeposits = cms.PSet()
process.allLayer1Muons.isoDeposits = cms.PSet()
#embed the inner track of the muon: --Jochen
process.allLayer1Muons.embedTrack = cms.bool(True)

# Jets
from PhysicsTools.PatAlgos.tools.jetTools import *
setTagInfos(process,
            coll = "allLayer1Jets",
            tagInfos = cms.vstring("secondaryVertex")
            )

for jetcoll in (process.allLayer1JetsTopTagCalo, process.allLayer1JetsTopTagPF):
	jetcoll.embedGenJetMatch = cms.bool(False)
	#getJetMCFlavour uses jetFlavourAssociation*, see below
	jetcoll.getJetMCFlavour = cms.bool(True)
	#those two use jetPartonMatch*, see below
	jetcoll.addGenPartonMatch = cms.bool(True)
	jetcoll.embedGenPartonMatch = cms.bool(True)
	# Add CATopTag info... piggy-backing on b-tag functionality
	jetcoll.discriminatorSources = cms.VInputTag()
	jetcoll.addBTagInfo = cms.bool(True)
	jetcoll.addTagInfos = cms.bool(True)
	jetcoll.addDiscriminators = cms.bool(False)

process.allLayer1JetsTopTagCalo.tagInfoSources = cms.VInputTag( cms.InputTag('CATopCaloJetTagInfos') )
process.allLayer1JetsTopTagPF.tagInfoSources = cms.VInputTag( cms.InputTag('CATopPFJetTagInfos') )

#adapt jet parton stuff for topjet collections. There is only one "jetPartons" module ...
process.jetPartons.withTop = cms.bool(True)
# ... but many "jetPartonAssociation", "jetFlavourAssociation" and "jetPartonMatch":
for jetcollname in ('TopTagCalo', 'TopTagPF'):
	# definition 4 = match heaviest flavor 
	getattr(process, 'jetFlavourAssociation'+jetcollname).definition = cms.int32(4)
	getattr(process, 'jetFlavourAssociation'+jetcollname).physicsDefinition = cms.bool(False)
	# jetFlavourAssociation depends on jetPartonAssociation, so re-configure that as well:
	getattr(process, 'jetPartonAssociation'+jetcollname).coneSizeToAssociate = cms.double(0.8)
	getattr(process, 'jetPartonAssociation'+jetcollname).doPriority = cms.bool(True)
	getattr(process, 'jetPartonAssociation'+jetcollname).priorityList = cms.vint32(6)
	getattr(process, 'jetPartonMatch'+jetcollname).mcPdgId = cms.vint32(1,2,3,4,5,6,21)
	getattr(process, 'jetPartonMatch'+jetcollname).maxDeltaR = cms.double(0.8)

#Note: the original jetPartonAssociation, jetFlavourAssociation, jetPartonMatch which is used for antiktcalo
# is not changed.

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
				   #(cmp. TopQuarkAnalysis/TopEventProducers/python/tqafEventContent_cff.py):
                                   'keep *_decaySubset_*_*',
                                   'keep *_initSubset_*_*',
				   'keep *_genEvt_*_*',
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
				   'drop *_cleanLayer1Electrons_*_*',
				   'drop *_cleanLayer1Muons_*_*',
				   'keep *_selectedLayer1Electrons_*_*',
				   'keep *_selectedLayer1Muons_*_*',
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
		     #addJetCollection only adds the "allLayer1Jets*" to the sequence, so 
                     # add the selection manually: -- Jochen
                     process.selectedLayer1JetsTopTagCalo *
                     process.selectedLayer1JetsTopTagPF *
                     process.countLayer1Jets
                     )

if type!='ttbar' :
    process.p.remove( process.makeGenEvt )

process.source.fileNames = [
    #TODO: fill in your dataset filenames here
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/F8C20DCD-5FCA-DE11-AFF0-001E0BE922E2.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/C626BA4F-F4CA-DE11-89B4-00237DA1FC56.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/B6A540D3-30CA-DE11-9CDE-001E0B5FC57A.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/7EFDA264-63CA-DE11-9DD1-001F296B9576.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/6C730F68-C6CA-DE11-8F8F-001CC4AADC6E.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/507CAF87-4CCA-DE11-A4DB-001E0B5FD4A6.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0001/22C12A62-80CA-DE11-8228-00237DA1CD92.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/FE24C6A9-1CCA-DE11-8842-00237DA1FC56.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/EAB6EEDB-11CA-DE11-9FF8-001F29C4D35E.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/CC984B69-18CA-DE11-A3C0-001F296A7698.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/9AE0D101-1BCA-DE11-89DF-00215AAC88D2.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/8C71AFC4-29CA-DE11-9233-001CC4AADC6E.root',
       '/store/mc/Summer09/QCDDiJet_Pt2600to3000/GEN-SIM-RECO/MC_31X_V9_7TeV-v1/0000/36203CFC-16CA-DE11-B2D4-00237DA15C96.root'
    ]

#On MC, there are often non-unique run and event ids. Safeguard against skipping in that case: --Jochen
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents.input = cms.untracked.int32(500)         ##  (e.g. -1 to run on all events)

#override settings in CMSSW/ PhysicsTools/ PatAlgos/ python/ recoLayer0/ photonIsolation_cff.py, where
#it is (wrongly) assumed that the reconstruction process label is called "RECO": for fastsim, it is usually called "HLT". Therefore,
# omit the process label. Maybe this breaks some things for 31x (TODO: someone could check that), so do it only for run on 33x: --Jochen
if runon=='33x':
    process.gamIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
    process.gamIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")
