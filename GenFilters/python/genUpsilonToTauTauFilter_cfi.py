#Configuration file fragment used for GenUpsilonToTauTauFilter module (GeneratorInterface/GenFilters/src/GenUpsilonToTauTauFilter.cc) initalisation
#genUpsilonToTauTauFilter_cfi GeneratorInterface/GenFilters/python/genUpsilonToTauTauFilter_cfi.py

import FWCore.ParameterSet.Config as cms

genUpsilonToTauTauFilter = cms.EDFilter("GenUpsilonToTauTauFilter",
   src       = cms.InputTag("genParticles"), #GenParticles collection as input
   tauPtCut  = cms.double(10.0), #at least a GenTau with this minimum pT 
   tauEtaCut = cms.double(2.4),  #GenTau eta
   nUpsilons = cms.double(1),    #GenTau eta
   taudRCut  = cms.double(0.4)   #GenTauTau cut
)
