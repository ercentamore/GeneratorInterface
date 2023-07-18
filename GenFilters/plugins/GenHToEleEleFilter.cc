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
class GenHToEleEleFilter : public edm::global::EDFilter<> {
public:
  explicit GenHToEleEleFilter(const edm::ParameterSet&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  //Member data
  const edm::EDGetTokenT<reco::GenParticleCollection> token_;
  const double elePtCut_, eleEtaCut_, nHiggs_, eledRCut_ ;
};

//Constructor
GenHToEleEleFilter::GenHToEleEleFilter(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
      elePtCut_(params.getParameter<double>("elePtCut")),
      eleEtaCut_(params.getParameter<double>("eleEtaCut")),
      nHiggs_(params.getParameter<double>("nHiggs")),
      eledRCut_(params.getParameter<double>("eledRCut")) {}

bool GenHToEleEleFilter::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
  using namespace std;
  using namespace edm;
  using namespace reco;

  //Read GenParticles Collection from Event
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  //Loop over all eles in Event
  unsigned HToEleEleCandidate = 0;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 11 || abs(iGen->daughter(1)->pdgId()) != 11 ) continue;
    if ( iGen->daughter(0)->pt() < elePtCut_ && iGen->daughter(1)->pt() < elePtCut_ ) continue;
    if ( iGen->daughter(0)->eta() > eleEtaCut_ || iGen->daughter(1)->eta() > eleEtaCut_ ) continue;
    float deltaeta = fabs(iGen->daughter(0)->eta()-iGen->daughter(1)->eta());
    float deltaphi = fabs(iGen->daughter(0)->phi()-iGen->daughter(1)->phi());
    float deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    if ( deltaR > eledRCut_ ) continue;
    ++HToEleEleCandidate;
  }
  return (HToEleEleCandidate >= nHiggs_);  //Return boolean whether event passes cut values
}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenHToEleEleFilter);
