import FWCore.ParameterSet.Config as cms

TopAnalyzer = cms.EDAnalyzer("BooLowMAnalyzer",
                             debug = cms.bool( False ),
                             IsMCTop = cms.bool( True ),
                             leptonFlavor = cms.int32( 13 ),
                             genEventSource = cms.InputTag('genEvt'),
                             muonSource    = cms.InputTag('selectedLayer1Muons'),
                             electronSource = cms.InputTag('selectedLayer1Electrons'),
                             metSource      = cms.InputTag('selectedLayer1METs'),
                             jetSource      = cms.InputTag('selectedLayer1Jets'),
                             EvtSolution    = cms.InputTag('solutions::TEST'),
                             rootFilename = cms.string('TopAnalysis.root'),
                             jetCuts = cms.PSet(
                                       MinJetEt        = cms.double( 30. ),
                                       MinJetEta       = cms.double( 2.4)
                                       ),
                             muonCuts = cms.PSet(
                                       MinPt  = cms.double( 20. ),
                                       MinEta = cms.double( 2.1 ),
                                       RelIso = cms.double( 0.95 ),
                                       MinCaloEnergy = cms.double( 0. )
                                       ),
                             muonIsolation = cms.PSet(
                                       RelIso = cms.double( 0.95 ),
                                       MaxVetoEm = cms.double( 4.0 ),
                                       MaxVetoHad = cms.double( 6.0 )
                                       ),
                             electronCuts = cms.PSet(
                                       MinPt  = cms.double( 20. ),
                                       MinEta = cms.double( 2.4 ),
                                       RelIso = cms.double( 0.95 )
                                       ),
                             METCuts = cms. PSet(
                                       MinMET = cms.double( 0. )
                                       ),
                             writeAscii = cms.bool( False),
                             asciiFilename = cms.string('TopAnalysis.txt'),
                             processOnlyEvent = cms.int32( -1 ),
                             makeJetLegoPlots = cms.bool( False ),
                             )

