# Configuration file fragment used for GenHToAATo4TauFilter module (GeneratorInterface/GenFilters/src/GenHToAATo4TauFilter.cc) initalisation
#genHToAATo4TauFilter_cfi GeneratorInterface/GenFilters/python/genHToTo4TauFilter_cfi.py

import FWCore.ParameterSet.Config as cms




process.genHToAATo4TauFilter = cms.EDFilter("GenHToAATo4TauFilter",
   src       = cms.InputTag("genParticles"), #GenParticles collection as input
   nHiggs    = cms.double(1),    #Number of higgs candidates
   tauPtCut  = cms.double(10.0), #at least a GenTau with this minimum pT
   tauEtaCut = cms.double(2.4),    #GenTau eta max value
   taudRCut  = cms.double(0.4)   #GenTauTau dR max value
)
