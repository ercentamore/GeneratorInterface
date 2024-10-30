

//System include files
#include <memory>
#include <vector>

//User include files
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

double assign_dr_value(double mass, double deltsR_cut, double a_pt) {
    if (mass >= 1.2 && mass < 1.6 && a_pt <= 50) {
        return deltsR_cut+0.01;
    }
    else if (mass >= 1.6 && mass < 2.0 && a_pt <= 90) {
        return deltsR_cut+0.03;
    }
    else if (mass >= 2.0 && mass < 2.4 && a_pt <= 130) {
        return deltsR_cut+0.04;
    }
    else if (mass >= 2.4 && mass < 2.8 && a_pt <= 150) {
        return deltsR_cut+0.10;
    }
    else if (mass >= 2.8 && mass < 3.2 && a_pt <= 180) {
        return deltsR_cut+0.14;
    }
    else if (mass >= 3.2 && mass <= 3.6 && a_pt <= 200) {
        return deltsR_cut+0.16;
    }
    else {
        std::cerr << "Default value selected " << std::endl;
        return deltsR_cut; 
    }
}

//Class declaration
class GenAToTauTauFilter_unphysical : public edm::global::EDFilter<> {
public:
  explicit GenAToTauTauFilter_unphysical(const edm::ParameterSet&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  //Member data
  const edm::EDGetTokenT<reco::GenParticleCollection> token_;
  const double tauPtCut_, tauEtaCut_, nHiggs_, taudRCut_ ;
};

//Constructor
GenAToTauTauFilter_unphysical::GenAToTauTauFilter_unphysical(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
      tauPtCut_(params.getParameter<double>("tauPtCut")),
      tauEtaCut_(params.getParameter<double>("tauEtaCut")),
      nHiggs_(params.getParameter<double>("nHiggs")),
      taudRCut_(params.getParameter<double>("taudRCut")) {}

bool GenAToTauTauFilter_unphysical::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
  using namespace std;
  using namespace edm;
  using namespace reco;

  //Read GenParticles Collection from Event
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  //Loop over all taus in Event
  unsigned HToTauTauCandidate = 0;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Initialize filter GenAToTauTauFilter_unphysical  :   " << endl;
    if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->status()) != 2 || abs(iGen->daughter(1)->status()) != 2 || iGen->daughter(0)->numberOfMothers() < 1 || iGen->daughter(1)->numberOfMothers() < 1) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
    if ( iGen->daughter(0)->pt() < tauPtCut_ && iGen->daughter(1)->pt() < tauPtCut_ ) continue;
    if ( abs(iGen->daughter(0)->eta()) > tauEtaCut_ && abs(iGen->daughter(1)->eta()) > tauEtaCut_ ) continue;
    float deltaR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi());
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  A mass  :   " << iGen->mass()  << endl;
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  A mass  :   " << iGen->mass()<<" <<< Tau1 mass  "<<iGen->daughter(0)->mass()<<" <<< Tau2 mass:  "<<iGen->daughter(1)->mass() << endl;
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> A charge  <<<<<"<<iGen->charge()<<" <<< Tau1 charge  "<<iGen->daughter(0)->charge()<<" <<< Tau2 charge:  "<<iGen->daughter(1)->charge()<< std::endl;
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> A pdgId  <<<<<"<<iGen->pdgId()<<" <<< Tau1 pdgId  "<<iGen->daughter(0)->pdgId()<<" <<< Tau2 pdgId:  "<<iGen->daughter(1)->pdgId()<< std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  deltaR  :   " << deltaR  << endl;
    double final_taudRCut = assign_dr_value(iGen->mass(), taudRCut_,iGen->pt() );
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  final_taudRCut  :   " << final_taudRCut  << endl;
    if ( deltaR > final_taudRCut ) continue; //merged
    ++HToTauTauCandidate;

    std::cout<<"========================================================event passed filter======================================================="<<endl;
    // std::cout << "deltsR  :   " << deltaR  << endl;
  }
  return (HToTauTauCandidate >= nHiggs_);  //Return boolean whether event passes cut values


}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenAToTauTauFilter_unphysical);
