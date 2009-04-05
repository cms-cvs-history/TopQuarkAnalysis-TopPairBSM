import FWCore.ParameterSet.Config as cms
#
# pat tuple sequence for fullsim
#

from PhysicsTools.PatAlgos.patLayer0_cff                           import *
from PhysicsTools.PatAlgos.patLayer1_cff                           import *
from TopQuarkAnalysis.TopObjectProducers.patTuple_Defaults_cff     import *
#from TopQuarkAnalysis.TopObjectProducers.patTuple_isoDeposits_cff  import *
#from TopQuarkAnalysis.TopObjectProducers.patTuple_caloTaus_cff    import *
from PhysicsTools.HepMCCandAlgos.genEventProcID_cfi                import *
from PhysicsTools.HepMCCandAlgos.genEventRunInfo_cfi               import *

## std sequence for patTuple production
BooTopPatTuple = cms.Sequence(genEventProcID +             ## needs HepMCProduct in the event content
                        patLayer0_patTuple *         ## to be used from PhysicsTools/PatAlgos 
                        patLayer1
                        )

## reduced sequence for madgraph samples which do not contain the HepMCProduct
BooTopPatTuple_reduced = cms.Sequence(
    patLayer0_patTuple *         ## to be used from PhysicsTools/PatAlgos
    patLayer1# *                 ## V04-14-03 onwards
    #                       patCaloTaus                  ## skipped as long as common taus are caloTaus
    )
