/*
Original Author:  Davide Di Croce
         Created:  Fev 2021
*/

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

//Class declaration
class GenUpsilonToTauTauFilter : public edm::global::EDFilter<> {
public:
  explicit GenUpsilonToTauTauFilter(const edm::ParameterSet&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  //Member data
  const edm::EDGetTokenT<reco::GenParticleCollection> token_;
  const double tauPtCut_, tauEtaCut_, nUpsilons_, taudRCut_ ;
};

//Constructor
GenUpsilonToTauTauFilter::GenUpsilonToTauTauFilter(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
      tauPtCut_(params.getParameter<double>("tauPtCut")),
      tauEtaCut_(params.getParameter<double>("tauEtaCut")),
      nUpsilons_(params.getParameter<double>("nUpsilons")),
      taudRCut_(params.getParameter<double>("taudRCut")) {}

bool GenUpsilonToTauTauFilter::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
  using namespace std;
  using namespace edm;
  using namespace reco;

  //Read GenParticles Collection from Event
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  //Loop over all taus in Event
  unsigned UpsilonToTauTauCandidate = 0;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 553 || iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
    if ( iGen->daughter(0)->pt() < tauPtCut_ && iGen->daughter(1)->pt() < tauPtCut_ ) continue;
    if ( abs(iGen->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(1)->eta()) > tauEtaCut_ ) continue;
    float deltaeta = fabs(iGen->daughter(0)->eta()-iGen->daughter(1)->eta());
    float deltaphi = fabs(iGen->daughter(0)->phi()-iGen->daughter(1)->phi());
    float deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    if ( deltaR > taudRCut_ ) continue;
    ++UpsilonToTauTauCandidate;
    std::cout<<"========================================================event passed filter======================================================="<<endl;
    std::cout << "deltsR  :   " << deltaR << endl;
  }

  return (UpsilonToTauTauCandidate >= nUpsilons_);  //Return boolean whether event passes cut values
}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenUpsilonToTauTauFilter);
