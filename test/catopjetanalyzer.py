# import configurations
import FWCore.ParameterSet.Config as cms

print "About to process"

# define the process
process = cms.Process("CATopJets")

print "Creating message logger"
# input message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.include( "SimGeneral/HepPDTESSource/data/pythiapdt.cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

print "Setting variables"

dataset = 'qcd_470'
algorithm = 'ca'
output_dst = True
nevents = 200
idtag = 'default'
outputdir = '/uscms_data/d1/rappocc/'
inputtype = 'CaloJet'

# this defines the input files

if dataset == 'zprime' :
    from TopQuarkAnalysis.TopPairBSM.RecoInput_ZPrime2000_cfi import *
elif dataset == 'qcd_smallstats' :
    from TopQuarkAnalysis.TopPairBSM.RecoInput_QCD_500_1000_cfi import *
elif dataset == 'qcd_470' :
    from TopQuarkAnalysis.TopPairBSM.RecoInput_QCD_470_cfi import *
elif dataset == 'qcd_600' :
    from TopQuarkAnalysis.TopPairBSM.RecoInput_QCD_600_cfi import *
else :
    from TopQuarkAnalysis.TopPairBSM.RecoInput_ttbar_cfi import *

print "Dataset = " + dataset

# get generator sequences
#process.load("Configuration.StandardSequences.Generator_cff")
process.load("RecoJets.Configuration.GenJetParticles_cff")

# Load geometry
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.GeometryPilot1_cff")

# Kt6 jets
process.load("RecoJets.JetProducers.kt6CaloJets_cff")

# CATopJets
#process.load("TopQuarkAnalysis.TopPairBSM.caTopJets_cff")
process.load("TopQuarkAnalysis.TopPairBSM.CATopJetTagger_cfi")

from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
process.load("RecoJets.Configuration.GenJetParticles_cff")

print "Switching input collections"

parameters = cms.PSet()

if inputtype == 'GenJet' :
    parameters = GenJetParameters
else :
    parameters = CaloJetParameters

parameters.jetPtMin = cms.double(500.)
parameters.correctInputToSignalVertex = cms.bool(False)
parameters.inputEtMin = cms.double(5.0)


process.caTopJetsProducer = cms.EDProducer("CATopJetProducer",
                                           CATopJetParameters,
                                           parameters
                                           )

print "About to input pat sequences"

# input pat sequences
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

# switch jet collection to our juets
from PhysicsTools.PatAlgos.tools.jetTools import *

print "About to switch jet collection"

## ==== Example with CaloJets
switchJetCollection(process, 
        'caTopJetsProducer',   # Jet collection; must be already in the event when patLayer0 sequence is executed
        layers=[0,1],          # If you're not runnint patLayer1, set 'layers=[0]' 
        runCleaner="BasicJet", # =None if not to clean
        doJTA=True,            # Run Jet-Track association & JetCharge
        doBTagging=True,       # Run b-tagging
        jetCorrLabel='FKt6',   # example jet correction name; set to None for no JEC
        doType1MET=False)      # recompute Type1 MET using these jets

# Place appropriate jet cuts (NB: no cut on number of constituents)
process.selectedLayer1Jets.cut = cms.string('et > 500. & abs(eta) < 5.0')
# Turn off resolutions, they don't mean anything here
process.allLayer1Jets.addResolutions = cms.bool(False)
# Add CATopTag info... piggy-backing on b-tag functionality
process.layer0TagInfos.associations.append( cms.InputTag("CATopJetTagger") )
process.allLayer1Jets.addBTagInfo = cms.bool(True)
process.allLayer1Jets.addTagInfoRefs = cms.bool(True)
process.allLayer1Jets.tagInfoNames = cms.vstring('CATopJetTagger')
process.allLayer1Jets.addDiscriminators = cms.bool(False)
# Add parton match to quarks and gluons
process.allLayer1Jets.addGenPartonMatch = cms.bool(True)
process.allLayer1Jets.embedGenPartonMatch = cms.bool(True)
process.jetPartons.withTop = cms.bool(True)
process.jetPartonAssociation.coneSizeToAssociate = cms.double(0.8)
process.jetPartonAssociation.doPriority = cms.bool(True)
process.jetPartonAssociation.priorityList = cms.vint32(6)
process.jetFlavourAssociation.definition = cms.int32(4)
process.jetFlavourAssociation.physicsDefinition = cms.bool(False)

process.jetPartonMatch.mcPdgId = cms.vint32(1,2,3,4,5,6,21)
process.jetPartonMatch.maxDeltaR = cms.double(0.8)
#process.allLayer1Jets.genPartonMatch = cms.InputTag("CAJetPartonMatcher")
# Add jet MC flavour (custom built to capture tops)
process.allLayer1Jets.getJetMCFlavour = cms.bool(True)


#process.allLayer1Jets.JetPartonMapSource = cms.InputTag("CAJetFlavourIdentifier")




print "Done switching jet collection"

# input pat analyzer sequence
process.load("TopQuarkAnalysis.TopPairBSM.CATopJetKit_cfi")

# load the pat layer 1 event content
process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")


#process.L2L3JetCorrectorFKt6 = cms.ESSource("JetCorrectionServiceChain",
#    correctors = cms.vstring('L2RelativeJetCorrectorFKt6', 
#        'L3AbsoluteJetCorrectorFKt6'),
#    label = cms.string('L2L3JetCorrectorFKt6')
#)
#process.es_prefer_L2L3JetCorrectorFKt6 = cms.ESPrefer("JetCorrectionServiceChain","L2L3JetCorrectorFKt6")
#process.L2JetCorJetFKt6.src = cms.InputTag("kt6CaloJets")




if algorithm == 'kt' :
    process.caTopJetsProducer.algorithm = cms.int32(0)
elif algorithm == 'ca' :
    process.caTopJetsProducer.algorithm = cms.int32(1)
elif algorithm == 'antikt' :
    process.caTopJetsProducer.algorithm = cms.int32(2)

# pythia output
process.printList = cms.EDAnalyzer( "ParticleListDrawer",
                                src = cms.InputTag( "genParticles" ),
                                maxEventsToPrint = cms.untracked.int32( 0 )
)



# define the source, from reco input
process.source = RecoInput()

# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevents)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histotesting_' + dataset + '_' + algorithm + '_' + inputtype + '_' + idtag + '.root')
)


# setup event content
process.patEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *')
    )

# extend event content to include PAT objects
process.patEventContent.outputCommands.extend(process.patLayer1EventContent.outputCommands)
process.patEventContent.outputCommands.extend(['keep *_genParticles_*_*',
                                               'keep *_genParticlesForJets_*_*',
                                               'keep *_towerMaker_*_*',
                                               'keep *_caTopJetsProducer_*_*',
                                               'keep *_CATopJetTagger_*_*',
                                               'drop *_selectedLayer1Taus_*_*',
                                               'drop *_selectedLayer1Hemispheres_*_*',
                                               'drop *_selectedLayer1Photons_*_*',
                                               'keep *_CAJetPartonMatcher_*_*',
                                               'keep *_CAJetFlavourIdentifier_*_*'])


# define event selection to be that which satisfies 'p'
process.patEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )
    )

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
                               process.patEventSelection,
                               process.patEventContent,
                               verbose = cms.untracked.bool(False),
                               fileName = cms.untracked.string(outputdir + dataset + '_' + algorithm + '_' + inputtype + '_' + idtag + '_testing.root')
                               )

# define path 'p'
process.p = cms.Path(process.genParticlesForJets*
                     process.kt6CaloJets*
                     process.printList*
                     process.caTopJetsProducer*
                     process.CATopJetTagger*
                     process.patLayer0*
                     process.patLayer1
#                     process.CATopJetKit
                     )
# define output path
if output_dst == True :
    process.outpath = cms.EndPath(process.out)

# Set the threshold for output logging to 'info'
process.MessageLogger.cerr.threshold = 'INFO'
