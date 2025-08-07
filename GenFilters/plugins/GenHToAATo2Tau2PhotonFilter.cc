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

//Class declaration
class GenHToAATo2Tau2PhotonFilter : public edm::global::EDFilter<> {
public:
    explicit GenHToAATo2Tau2PhotonFilter(const edm::ParameterSet&);

private:
    bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

    //Member data
    const edm::EDGetTokenT<reco::GenParticleCollection> token_;
    const double tauPtCut_, tauEtaCut_, phoPtCut_, phoEtaCut_, diPhoPtCut_, diTauPtCut_, nHiggs_, phoDrCut_ ;
};

//Constructor
GenHToAATo2Tau2PhotonFilter::GenHToAATo2Tau2PhotonFilter(const edm::ParameterSet& params)
    : token_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("src"))),
        tauPtCut_(params.getParameter<double>("tauPtCut")),
        tauEtaCut_(params.getParameter<double>("tauEtaCut")),
        phoPtCut_(params.getParameter<double>("phoPtCut")),
        phoEtaCut_(params.getParameter<double>("phoEtaCut")),
        diPhoPtCut_(params.getParameter<double>("diPhoPtCut")),
        diTauPtCut_(params.getParameter<double>("diTauPtCut")),
        nHiggs_(params.getParameter<int>("nHiggs")),
        phoDrCut_(params.getParameter<double>("phoDrCut")) {}

bool GenHToAATo2Tau2PhotonFilter::filter(edm::StreamID, edm::Event& evt, const edm::EventSetup& params) const {
    using namespace std;
    using namespace edm;
    using namespace reco;

    //Read GenParticles Collection from Event
    edm::Handle<reco::GenParticleCollection> genParticles;
    evt.getByToken(token_, genParticles);

    //Loop over all taus in Event
    unsigned HTo2Tau2PhotonCandidate = 0;
    // std::cout << " Applying GenHToAATo2Tau2PhotonFilter " << endl;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
        if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2) continue;
        if ( !((abs(iGen->daughter(0)->pdgId()) == 36 && abs(iGen->daughter(1)->pdgId()) == 35)|| (abs(iGen->daughter(0)->pdgId()) == 35 && abs(iGen->daughter(1)->pdgId()) == 36)) )  continue;
        if ( abs(iGen->daughter(0)->pdgId()) == 36 ){
            if (abs(iGen->daughter(0)->pt() < diPhoPtCut_)) continue;
            const reco::Candidate* pho1 = iGen->daughter(0)->daughter(0);
            const reco::Candidate* pho2 = iGen->daughter(0)->daughter(1);
            if ((abs(iGen->daughter(0)->daughter(0)->pdgId()) != 22) || (abs(iGen->daughter(0)->daughter(1)->pdgId()) != 22)) continue;
            if ((iGen->daughter(0)->daughter(0)->pt() < phoPtCut_) || (iGen->daughter(0)->daughter(1)->pt() < phoPtCut_)) continue;
            if ((abs(iGen->daughter(0)->daughter(0)->eta()) > phoEtaCut_) || (abs(iGen->daughter(0)->daughter(1)->eta()) > phoEtaCut_)) continue;
            double dR = reco::deltaR( pho1->eta(), pho1->phi(), pho2->eta(), pho2->phi() );
            if ( dR < phoDrCut_ ) continue;
          
        }
        if ( abs(iGen->daughter(0)->pdgId()) == 35 ){
            if (abs(iGen->daughter(0)->pt() < diTauPtCut_)) continue;
            if ((abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15) || (abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15)) continue;
            if ((iGen->daughter(0)->daughter(0)->pt() < tauPtCut_) || (iGen->daughter(0)->daughter(1)->pt() < tauPtCut_)) continue;
            if ((abs(iGen->daughter(0)->daughter(0)->eta()) > tauEtaCut_) || (abs(iGen->daughter(0)->daughter(1)->eta()) > tauEtaCut_)) continue;
        }
        if ( abs(iGen->daughter(1)->pdgId()) == 36 ){
            if (abs(iGen->daughter(0)->pt() < diPhoPtCut_)) continue;
            const reco::Candidate* pho1 = iGen->daughter(1)->daughter(0);
            const reco::Candidate* pho2 = iGen->daughter(1)->daughter(1);
            if ((abs(iGen->daughter(1)->daughter(0)->pdgId()) != 22) || (abs(iGen->daughter(1)->daughter(1)->pdgId()) != 22)) continue;
            if ((iGen->daughter(1)->daughter(0)->pt() < phoPtCut_) || (iGen->daughter(1)->daughter(1)->pt() < phoPtCut_)) continue;
            if ((abs(iGen->daughter(1)->daughter(0)->eta()) > phoEtaCut_) || (abs(iGen->daughter(1)->daughter(1)->eta()) > phoEtaCut_)) continue;
            double dR = reco::deltaR( pho1->eta(), pho1->phi(), pho2->eta(), pho2->phi() );
            if ( dR < phoDrCut_ ) continue;
        }
        if ( abs(iGen->daughter(1)->pdgId()) == 35 ){
            if (abs(iGen->daughter(0)->pt() < diTauPtCut_)) continue;
            if ((abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15) || (abs(iGen->daughter(1)->daughter(1)->pdgId() != 15))) continue;
            if ((iGen->daughter(1)->daughter(0)->pt() < tauPtCut_) || (iGen->daughter(1)->daughter(1)->pt() < tauPtCut_)) continue;
            if ((abs(iGen->daughter(1)->daughter(0)->eta()) > tauEtaCut_) || (abs(iGen->daughter(1)->daughter(1)->eta()) > tauEtaCut_)) continue;
        }
        ++HTo2Tau2PhotonCandidate;
        std::cout<<"========================================================event passed filter======================================================="<<endl;
    }
    return (HTo2Tau2PhotonCandidate >= nHiggs_);  //Return boolean whether event passes cut values
}

// Define module as a plug-in
DEFINE_FWK_MODULE(GenHToAATo2Tau2PhotonFilter);
