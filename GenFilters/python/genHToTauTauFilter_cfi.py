#Configuration file fragment used for GenHToTauTauFilter module (GeneratorInterface/GenFilters/src/GenHToTauTauFilter.cc) initalisation
#genHToTauTauFilter_cfi GeneratorInterface/GenFilters/python/genHToTauTauFilter_cfi.py

import FWCore.ParameterSet.Config as cms

genHToTauTauFilter = cms.EDFilter("GenHToTauTauFilter",
   src       = cms.InputTag("genParticles"), #GenParticles collection as input
   tauPtCut  = cms.double(10.0), #at least a GenTau with this minimum pT on each pair
   tauEtaCut = cms.double(2.4),  #GenTau eta max value
   nHiggs    = cms.double(2),    #number of A (higgs)
   taudRCut  = cms.double(0.4)   #GenTauTau dR cut
)
