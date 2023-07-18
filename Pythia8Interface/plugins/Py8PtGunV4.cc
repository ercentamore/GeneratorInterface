#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"

#include "TLorentzVector.h"

namespace gen {

class Py8PtGunV4 : public Py8GunBase {

   public:

      Py8PtGunV4( edm::ParameterSet const& );
      ~Py8PtGunV4() {}

      bool generatePartonsAndHadronize() override;
      const char* classname() const override;

   private:

      // PtGun particle(s) characteristics
      double  fMinEta;
      double  fMaxEta;
      double  fMinPt ;
      double  fMaxPt ;
      double  fMinMass ;
      double  fMaxMass ;
      double  fMaxdR ;
      bool    fAddAntiParticle;

};

// implementation
//
Py8PtGunV4::Py8PtGunV4( edm::ParameterSet const& ps )
   : Py8GunBase(ps)
{

   // ParameterSet defpset ;
   edm::ParameterSet pgun_params =
   ps.getParameter<edm::ParameterSet>("PGunParameters"); // , defpset ) ;
   fMinEta     = pgun_params.getParameter<double>("MinEta"); // ,-2.2);
   fMaxEta     = pgun_params.getParameter<double>("MaxEta"); // , 2.2);
   fMinPt      = pgun_params.getParameter<double>("MinPt"); // ,  0.);
   fMaxPt      = pgun_params.getParameter<double>("MaxPt"); // ,  0.);
   fMinMass    = pgun_params.getParameter<double>("MinMass"); // ,  0.);
   fMaxMass    = pgun_params.getParameter<double>("MaxMass"); // ,  0.);
   fMaxdR      = pgun_params.getParameter<double>("MaxdR"); // ,  0.);
   fAddAntiParticle = pgun_params.getParameter<bool>("AddAntiParticle"); //, false) ;

}

bool Py8PtGunV4::generatePartonsAndHadronize()
{

   fMasterGen->event.reset();

   for ( size_t i=0; i<fPartIDs.size(); i++ )
   {

      int particleID = fPartIDs[i]; // this is PDG - need to convert to Py8 ???
   
      double pt = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt; // flat

      double mass = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass; // flat mass

      // Calculate angles
      double dR = 99;
      double phi = 99;
      double eta = 99;
      double the = 99;
      while ( dR > fMaxdR )
      {
         phi = (fMaxPhi-fMinPhi) * randomEngine().flat() + fMinPhi;
         eta = (fMaxEta-fMinEta) * randomEngine().flat() + fMinEta;
         the = 2.*atan(exp(-eta));
         TLorentzVector AntiCandidate;
         AntiCandidate.SetPxPyPzE(-pt*cos(phi), -pt*sin(phi), -pt*cos(the)/sin(the), sqrt( pow(pt / sin(the),2) + mass*mass ));
         dR = sqrt( pow((eta-AntiCandidate.Eta()),2) + pow((phi-AntiCandidate.Phi()),2) );
      }

      // Calculate momenta
      double pp = pt / sin(the); // sqrt( ee*ee - mass*mass );
      double ee = sqrt( pp*pp + mass*mass );

      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pp * cos(the);

      if ( !((fMasterGen->particleData).isParticle( particleID )) )
      {
         particleID = std::fabs(particleID) ;
      }

      if( 1<= fabs(particleID) && fabs(particleID) <= 6) // quarks
        (fMasterGen->event).append( particleID, 23, 101, 0, px, py, pz, ee, mass );
      else if (fabs(particleID) == 21)                   // gluons
        (fMasterGen->event).append( 21, 23, 101, 102, px, py, pz, ee, mass );
      else                                               // other
        (fMasterGen->event).append( particleID, 1, 0, 0, px, py, pz, ee, mass );

      // Here also need to add anti-particle (if any)
      // otherwise just add a 2nd particle of the same type
      // (for example, gamma)
      if ( fAddAntiParticle )
      {
        if( 1 <= fabs(particleID) && fabs(particleID) <= 6){ // quarks
          (fMasterGen->event).append( -particleID, 23, 0, 101, -px, -py, -pz, ee, mass );
        }
        else if (fabs(particleID) == 21){                   // gluons
          (fMasterGen->event).append( 21, 23, 102, 101, -px, -py, -pz, ee, mass );
        }
        else if ( (fMasterGen->particleData).isParticle( -particleID ) ){
          (fMasterGen->event).append( -particleID, 1, 0, 0, -px, -py, -pz, ee, mass );
        }
        else {
          (fMasterGen->event).append( particleID, 1, 0, 0, -px, -py, -pz, ee, mass );
        }

      } // antiparticle

   } // fPartIDs

   if ( !fMasterGen->next() ) return false;

   event().reset(new HepMC::GenEvent);
   return toHepMC.fill_next_event( fMasterGen->event, event().get() );

} // generatePartonsAndHadronize()

const char* Py8PtGunV4::classname() const
{
   return "Py8PtGunV4";
}

typedef edm::GeneratorFilter<gen::Py8PtGunV4, gen::ExternalDecayDriver> Pythia8PtGunV4;

} // end namespace

using gen::Pythia8PtGunV4;
DEFINE_FWK_MODULE(Pythia8PtGunV4);
