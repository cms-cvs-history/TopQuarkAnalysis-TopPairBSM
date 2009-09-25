import FWCore.ParameterSet.Config as cms
from PhysicsTools.StarterKit.kinAxis_cfi  import *


CATopJetKit = cms.EDProducer("CATopJetKit",
    ntuplize = cms.string('all'),
    outputTextName = cms.string('CATopJetKit_output.txt'),
    enable = cms.string(''),
    disable = cms.string(''),
    doMuon         = cms.bool(True),
    doElectron     = cms.bool(True),
    doTau          = cms.bool(True),
    doJet          = cms.bool(True),
    doMET          = cms.bool(True),
    doPhoton       = cms.bool(True),
    doTrack        = cms.bool(False),
    doGenParticles = cms.bool(True),
    electronSrc = cms.InputTag("selectedLayer1Electrons"),
    tauSrc = cms.InputTag("selectedLayer1Taus"),
    muonSrc = cms.InputTag("selectedLayer1Muons"),
    jetSrc = cms.InputTag("selectedLayer1Jets"),
    photonSrc = cms.InputTag("selectedLayer1Photons"),
    METSrc = cms.InputTag("layer1METs"),
    trackSrc = cms.InputTag("allTrackCandidates"),
    genParticleSrc = cms.InputTag("genParticles"), 
    muonAxis     = kinAxis(0, 200, 0, 200),
    electronAxis = kinAxis(0, 200, 0, 200),
    tauAxis      = kinAxis(0, 200, 0, 200),
    jetAxis      = kinAxis(0, 200, 0, 200),
    METAxis      = kinAxis(0, 200, 0, 200),
    photonAxis   = kinAxis(0, 200, 0, 200),
    trackAxis    = kinAxis(0, 200, 0, 200)
)


