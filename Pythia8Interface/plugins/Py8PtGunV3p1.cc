#include "GeneratorInterface/Pythia8Interface/interface/Py8PtGun_added.h"
namespace gen {

class Py8PtGunV3p1 : public Py8GunBase {

   public:

      Py8PtGunV3p1( edm::ParameterSet const& );
      ~Py8PtGunV3p1() {}

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
      bool    fAddAntiParticle;

};

// implementation
//
Py8PtGunV3p1::Py8PtGunV3p1( edm::ParameterSet const& ps )
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
   fAddAntiParticle = pgun_params.getParameter<bool>("AddAntiParticle"); //, false) ;

}

std::vector <int> pT_bins_   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};
std::vector <double> m_bins_ = {4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.4, 8.8, 9.2, 9.6, 10.0, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.8, 13.2, 13.6, 14.0, 14.4, 14.8, 15.2, 15.6, 16.0, 16.4, 16.8, 17.2, 17.6, 18, 18.4, 18.8, 19.2, 19.6, 20 };
std::vector <int> occ_ = {
//35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  4.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  4.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  4.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  5.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  5.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  6.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  6.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  6.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  7.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  7.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  8.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  8.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  8.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  9.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  9.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  10.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  10.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  10.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  11.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  11.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  12.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  12.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  12.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  13.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  13.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  14.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  14.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  14.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  15.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  15.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  16.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  16.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  16.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  17.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  17.6
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  18.0
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  18.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  18.8
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  19.2
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1, // ->  19.4
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,    1,  1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,  1,   1,   1 // ->  20.0


};
std::vector <double> invpdf_ = get_inverse_pdf(occ_);


bool Py8PtGunV3p1::generatePartonsAndHadronize()
{

   //std::cout << " ! ! UPDATED PYTHIA VERSION ! ! " << std::endl;
   fMasterGen->event.reset();

   for ( size_t i=0; i<fPartIDs.size(); i++ )
   {

      int particleID = fPartIDs[i]; // this is PDG - need to convert to Py8 ???

      double rand_sampler = rand() / double(RAND_MAX);
      double pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
      double mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
      double weight       = lookup_invpdf(mass, m_bins_, pt, pT_bins_, invpdf_);
      while ( rand_sampler > weight ) {
         rand_sampler = rand() / double(RAND_MAX);
         pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
         mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
         weight       = lookup_invpdf(mass, m_bins_, pt, pT_bins_, invpdf_);
      }
      // Calculate angles
      double phi = (fMaxPhi-fMinPhi) * randomEngine().flat() + fMinPhi;
      double eta = (fMaxEta-fMinEta) * randomEngine().flat() + fMinEta;
      double the = 2.*atan(exp(-eta));

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

const char* Py8PtGunV3p1::classname() const
{
   return "Py8PtGunV3p1";
}

typedef edm::GeneratorFilter<gen::Py8PtGunV3p1, gen::ExternalDecayDriver> Pythia8PtGunV3p1;

} // end namespace

using gen::Pythia8PtGunV3p1;
DEFINE_FWK_MODULE(Pythia8PtGunV3p1);
