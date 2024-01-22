

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

#include "TH2D.h"

TH2D *A_mass_pt;
float a_mass;
float a_pt;
//Class declaration
class GenAToTauTauFilter : public edm::global::EDFilter<> {
public:
  explicit GenAToTauTauFilter(const edm::ParameterSet&);

private:
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  //Member data
  const edm::EDGetTokenT<reco::GenParticleCollection> token_;
  const double tauPtCut_, tauEtaCut_, nHiggs_, taudRCut_, a_mass_min_, a_mass_max_, mass_bins_,a_pt_min_, a_pt_max_, pt_bins_, n_Events_ ;
};

//Constructor
GenAToTauTauFilter::GenAToTauTauFilter(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
      tauPtCut_(params.getParameter<double>("tauPtCut")),
      tauEtaCut_(params.getParameter<double>("tauEtaCut")),
      nHiggs_(params.getParameter<double>("nHiggs")),
      taudRCut_(params.getParameter<double>("taudRCut")),

      a_mass_min_(params.getParameter<double>("a_mass_min")),
      a_mass_max_(params.getParameter<double>("a_mass_max")),
      mass_bins_(params.getParameter<double>("mass_bins")),
      a_pt_min_(params.getParameter<double>("a_pt_min")),
      a_pt_max_(params.getParameter<double>("a_pt_max")),
      pt_bins_(params.getParameter<double>("pt_bins")),
      n_Events_(params.getParameter<double>("n_Events"))
      {
        A_mass_pt  = new TH2D("h_a_mass_pt", "Mass vs pT of A;mass;pt" , mass_bins_,  a_mass_min_, a_mass_max_ , pt_bins_,  a_pt_min_, a_pt_max_ );
      }

bool GenAToTauTauFilter::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
  using namespace std;
  using namespace edm;
  using namespace reco;
  a_mass = -1;
  a_pt = -1;
  //Read GenParticles Collection from Event
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(token_, genParticles);

  //Loop over all taus in Event
  unsigned AToTauTauCandidate = 0;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
    if ( iGen->daughter(0)->pt() < tauPtCut_ && iGen->daughter(1)->pt() < tauPtCut_ ) continue;
    if ( abs(iGen->daughter(0)->eta()) > tauEtaCut_ || abs(iGen->daughter(1)->eta()) > tauEtaCut_ ) continue;
    float deltaeta = fabs(iGen->daughter(0)->eta()-iGen->daughter(1)->eta());
    float deltaphi = fabs(iGen->daughter(0)->phi()-iGen->daughter(1)->phi());
    float deltaR = sqrt(deltaeta*deltaeta+deltaphi*deltaphi);
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  A mass  :   " << iGen->mass()<<" <<< Tau1 mass  "<<iGen->daughter(0)->mass()<<" <<< Tau2 mass:  "<<iGen->daughter(1)->mass() << endl;
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> A charge  <<<<<"<<iGen->charge()<<" <<< Tau1 charge  "<<iGen->daughter(0)->charge()<<" <<< Tau2 charge:  "<<iGen->daughter(1)->charge()<< std::endl;
    // std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> A pdgId  <<<<<"<<iGen->pdgId()<<" <<< Tau1 pdgId  "<<iGen->daughter(0)->pdgId()<<" <<< Tau2 pdgId:  "<<iGen->daughter(1)->pdgId()<< std::endl;
    if ( deltaR > taudRCut_ ) continue;

    a_mass = iGen->mass();
    a_pt = iGen->pt();
    A_mass_pt->Fill(a_mass, a_pt);
    int binX = A_mass_pt->GetXaxis()->FindBin(a_mass);
    int binY = A_mass_pt->GetYaxis()->FindBin(a_pt);
    double binContent = A_mass_pt->GetBinContent(binX, binY);
    // std::cout << "binX :" << binX <<"   binY :  "<< binY<<std::endl;
    // std::cout << "Number of events in the bin: " << binContent <<"   mass:  "<< a_mass<<"  pt :  "<<a_pt<< std::endl;
    if (binContent > n_Events_) continue;
    ++AToTauTauCandidate;

    std::cout<<"========================================================event passed filter======================================================="<<endl;
    std::cout << "deltsR  :   " << deltaR  << endl;

  }
  return (AToTauTauCandidate >= nHiggs_);  //Return boolean whether event passes cut values

}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenAToTauTauFilter);
