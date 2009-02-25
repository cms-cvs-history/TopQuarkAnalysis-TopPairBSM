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
    #
    'file:/uscmst1b_scratch/lpc1/cmsroc/yumiceva/TQAF/challange2009/ChallengeFile.root',
    'file:/uscmst1b_scratch/lpc1/cmsroc/yumiceva/TQAF/challange2009/ChallengeFile2.root'
#    '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/00437D8F-04A9-DD11-8D99-0015C5E9B2AB.root',
#    '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/00BD3644-90A8-DD11-8D9E-001CC47BCFDC.root',
#    '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/02DCB67E-9FA9-DD11-B550-00E081402E8B.root',
#    '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/0419137A-9BA8-DD11-B5AC-0015C5E9C17C.root',
#    '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/046A072B-31A9-DD11-A931-00E08140EAB7.root'
    #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/16AAC418-218A-DD11-AC33-001F2908F0E4.root',
    #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/1E19C1C2-EF89-DD11-A6AB-001E0B1C74DA.root',
    #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/2AE099C8-1F8A-DD11-B30F-00144F2031D4.root',
    #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/2CA4D6BE-EF89-DD11-AEFA-001F290860E4.root',
    #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/32EE38A2-188A-DD11-827E-00144F283544.root'
    )
    , duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")
## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")


#-------------------------------------------------
# patTuple configuration
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopPairBSM.BooTopPatTuple")

## necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

## switch from clusters to rec hits in ECAL electron isolation
#from PhysicsTools.PatAlgos.recoLayer0.electronIsolation_cff import useElectronRecHitIsolation
#useElectronRecHitIsolation(process)
## switch from clusters to rec hits in ECAL photon isolation
#from PhysicsTools.PatAlgos.recoLayer0.photonIsolation_cff import usePhotonRecHitIsolation
#usePhotonRecHitIsolation(process)


#-------------------------------------------------
# process paths;
#-------------------------------------------------

## process path
process.p = cms.Path(process.BooTopPatTuple)

#-------------------------------------------------
# pat tuple event content; first ALL objects
# are dropped in this process; then patTuple
# content is added
#-------------------------------------------------

## define pat tuple event content
from TopQuarkAnalysis.TopObjectProducers.patTuple_EventContent_cff import *
makePatTupleEventContent(process)

## change jet collection
from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process, 
                    'sisCone5CaloJets',             # jet collection; must be already in the event when patLayer0 sequence is executed
                    layers       = [0,1],           # if you're not running patLayer1, set 'layers=[0]' 
                    runCleaner   = "CaloJet",       # =None if not to clean
                    doJTA        = True,            # run jet-track association & JetCharge
                    doBTagging   = True,            # run b-tagging
                    jetCorrLabel = ('SC5', 'Calo'), # example jet correction name; set to None for no JEC
                    doType1MET   = True             # recompute Type1 MET using these jets
                    )

#addJetCollection(process, "", "JPT", layers....  )


# selection
#
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi  import selectedLayer1Electrons
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi      import selectedLayer1Muons
from PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi       import selectedLayer1Taus
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi       import selectedLayer1Jets
from PhysicsTools.PatAlgos.selectionLayer1.metSelector_cfi       import selectedLayer1METs
from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import countLayer1Leptons
from PhysicsTools.PatAlgos.selectionLayer1.jetMinFilter_cfi      import minLayer1Jets
from PhysicsTools.PatAlgos.selectionLayer1.muonMinFilter_cfi     import minLayer1Muons
#

#selectedLayer1Electrons.cut  = cms.string('pt > 15. & abs(eta) < 2.5')
selectedLayer1Muons.cut      = cms.string('pt > 15. & abs(eta) < 2.4')
selectedLayer1Jets.cut       = cms.string('pt > 20. & abs(eta) < 2.4 & nConstituents > 0')
selectedLayer1METs.cut       = cms.string('et >= 0.')
#countLayer1Leptons.minNumber = 1
minLayer1Jets.minNumber      = 1
minLayer1Muons               = 1

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
    process.patTupleEventContent,
    verbose = cms.untracked.bool(True),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                           
    fileName = cms.untracked.string('/uscmst1b_scratch/lpc1/cmsroc/yumiceva/TQAF/challange2009/PatChallengeTestall.root'),
    dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('USER'),
            filterName = cms.untracked.string('')
                )
)

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1.2.1 $'),
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
