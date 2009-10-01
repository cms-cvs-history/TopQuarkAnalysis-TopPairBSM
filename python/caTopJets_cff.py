import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *


caTopCaloJets = cms.EDProducer("CATopJetProducer",
                               CATopJetParameters,
                               CaloJetParameters
                               )


caTopGenJets = cms.EDProducer("CATopJetProducer",
                              CATopJetParameters,
                              GenJetParameters
                              )



caTopPFJets = cms.EDProducer("CATopJetProducer",
                             CATopJetParameters,
                             PFJetParameters
                             )

