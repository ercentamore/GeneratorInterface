

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
class GenHToAATo4TauFilter : public edm::global::EDFilter<> {
public:
  explicit GenHToAATo4TauFilter(const edm::ParameterSet&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  //Member data
  const edm::EDGetTokenT<reco::GenParticleCollection> token_;
  const double tauPtCut_, tauEtaCut_, nHiggs_, taudRCut_ ;
};

//Constructor
GenHToAATo4TauFilter::GenHToAATo4TauFilter(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
      tauPtCut_(params.getParameter<double>("tauPtCut")),
      tauEtaCut_(params.getParameter<double>("tauEtaCut")),
      nHiggs_(params.getParameter<double>("nHiggs")),
      taudRCut_(params.getParameter<double>("taudRCut")) {}

bool GenHToAATo4TauFilter::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
  using namespace std;
  using namespace edm;
  using namespace reco;

  //Read GenParticles Collection from Event
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  //Loop over all taus in Event
  unsigned HToTauTauCandidate = 0;
  /*for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
    if ( iGen->daughter(0)->pt() < tauPtCut_ && iGen->daughter(1)->pt() < tauPtCut_ ) continue;
    if ( iGen->daughter(0)->eta() > tauEtaCut_ || iGen->daughter(1)->eta() > tauEtaCut_ ) continue;
    float deltaeta = fabs(iGen->daughter(0)->eta()-iGen->daughter(1)->eta());
    float deltaphi = fabs(iGen->daughter(0)->phi()-iGen->daughter(1)->phi());
    float deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    if ( deltaR > taudRCut_ ) continue;
    ++HToTauTauCandidate;
  }
  return (HToTauTauCandidate >= nHiggs_);  //Return boolean whether event passes cut values*/
//--------------------------------------------------------------------------------------------------------------------------------------
  // std::cout << " Applying GenHToAATo4TauFilter " << endl;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
    if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
    // if ( (iGen->daughter(0)->daughter(0)->pt() < tauPtCut_ && iGen->daughter(0)->daughter(1)->pt() < tauPtCut_) || (iGen->daughter(1)->daughter(0)->pt() < tauPtCut_ && iGen->daughter(1)->daughter(1)->pt() < tauPtCut_ )) continue;
    if ( (iGen->daughter(0)->daughter(0)->pt() < tauPtCut_ || iGen->daughter(0)->daughter(1)->pt() < tauPtCut_) || (iGen->daughter(1)->daughter(0)->pt() < tauPtCut_ || iGen->daughter(1)->daughter(1)->pt() < tauPtCut_ )) continue;
    // if ( (abs(iGen->daughter(0)->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(0)->daughter(1)->eta()) > tauEtaCut_) && (abs(iGen->daughter(1)->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(1)->daughter(1)->eta()) > tauEtaCut_)) continue;
    if ( (abs(iGen->daughter(0)->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(0)->daughter(1)->eta()) > tauEtaCut_) || (abs(iGen->daughter(1)->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(1)->daughter(1)->eta()) > tauEtaCut_)) continue;
    float deltaeta1 = fabs(iGen->daughter(0)->daughter(0)->eta()-iGen->daughter(0)->daughter(1)->eta());
    float deltaphi1 = fabs(iGen->daughter(0)->daughter(0)->phi()-iGen->daughter(0)->daughter(1)->phi());
    float deltaeta2 = fabs(iGen->daughter(1)->daughter(0)->eta()-iGen->daughter(1)->daughter(1)->eta());
    float deltaphi2 = fabs(iGen->daughter(1)->daughter(0)->phi()-iGen->daughter(1)->daughter(1)->phi());

    float deltaR1 = sqrt(deltaeta1*deltaeta1 + deltaphi1*deltaphi1);
    float deltaR2 = sqrt(deltaeta2*deltaeta2 + deltaphi2*deltaphi2);
    if ( deltaR1 > taudRCut_ || deltaR2 > taudRCut_) continue; // meged taus
    // if !( deltaR1 > taudRCut_ && deltaR2 > taudRCut_) continue; // unmerged taus
    // std::cout << " deltaR1:------- " << deltaR1 << "deltaR2:------- "<<deltaR2 << endl;
    // std::cout << " A1 mass:------- " << iGen->daughter(0).mass() << "A2 mass:------- "<<iGen->daughter(1).mass() << endl;
    ++HToTauTauCandidate;
    // std::cout << " Passed Through GenHToAATo4TauFilter " << endl;
    std::cout<<"========================================================event passed filter======================================================="<<endl;
    std::cout << " eta1  : " << iGen->daughter(0)->daughter(0)->eta()<< " eta2  : "<< iGen->daughter(0)->daughter(1)->eta() <<" eta3  : "<< iGen->daughter(1)->daughter(0)->eta() <<" eta4  : "<< iGen->daughter(0)->daughter(1)->eta() << std::endl;
    std::cout << " status1  : " << iGen->daughter(0)->daughter(0)->status()<< " status2  : "<< iGen->daughter(0)->daughter(1)->status() <<" status3  : "<< iGen->daughter(1)->daughter(0)->status() <<" status4  : "<< iGen->daughter(0)->daughter(1)->status() << std::endl;
    std::cout << "deltsR1  :   " << deltaR1 <<"            "<< "deltar2  :   " << deltaR2 << endl;
  }
  std::cout << "HToTauTauCandidate  :   " << HToTauTauCandidate << endl;
  return (HToTauTauCandidate >= nHiggs_);  //Return boolean whether event passes cut values
//---------------------------------------------------------------------------------------------------------------------------------------

}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenHToAATo4TauFilter);
