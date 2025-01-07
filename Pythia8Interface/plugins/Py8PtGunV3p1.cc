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

std::vector <int> pT_bins_   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300};
std::vector <double> m_bins_ = {1.6, 2.0, 2.4, 2.8, 3.2, 3.6};
std::vector <int> occ_ = {
//35,   40,   45,   50,   55,   60,   65,   70,   75,   80,  85,  90,   95,  100,    105,  110,  115,  120,  125,  130,  135,  140, 145,  150,   155,  160,  165,  170,  175,  180,  185,  190,  195,  200  205   210   215    220   225   230   235   240  245    250   255  260    265  270    275   280   285   290    295  300
6910, 6620, 7058, 6600, 6434, 6772, 6690, 6510, 6580, 6812, 6534, 6840, 6470, 6640, 6828, 6652, 6496, 6948, 6324, 6808, 6550, 6824, 6538, 6586, 6398, 6770, 6522, 6682, 6560, 6492, 6328, 6688, 6656, 6892, 6524, 6540, 6720, 6890, 6784, 6762, 6532, 6778, 6688, 6968, 6554, 6930, 6500, 6598, 6720, 6696, 6554, 6418, 6682, 6850,
1768, 12438, 10000, 8758, 6930, 7000, 6592, 6982, 6566, 6690, 6712, 6798, 6828, 6682, 6676, 7048, 6820, 6608, 6652, 6676, 6556, 6564, 6678, 6574, 6546, 6790, 6838, 6498, 6314, 6852, 6636, 6262, 7000, 7018, 6686, 6774, 6448, 6598, 6584, 6428, 6498, 6644, 6608, 6498, 6614, 6674, 6808, 6544, 6996, 6722, 6440, 6674, 6458, 6634,
324, 510, 26770, 43228, 33202, 13494, 7742, 7038, 6754, 6508, 6386, 6764, 6720, 6622, 6798, 6472, 6558, 6520, 6678, 6708, 6836, 6874, 6540, 6682, 6718, 6220, 6530, 6552, 6700, 6374, 6896, 6524, 6792, 6638, 6666, 6604, 6438, 6694, 6710, 6304, 6594, 6528, 6726, 6864, 6944, 6616, 6592, 6782, 6786, 6618, 6792, 6684, 6682, 6514,
1824, 5508, 6458, 24172, 8118, 6496, 6702, 6674, 6688, 6508, 6580, 6616, 6750, 7022, 6672, 6622, 6620, 6492, 6866, 6914, 6460, 6836, 6692, 6756, 6736, 6802, 6548, 6800, 6522, 6856, 6728, 6222, 6460, 6436, 6554, 6708, 6450, 6556, 6848, 6546, 6934, 6718, 6482, 6300, 6812, 6688, 6572, 6796, 6922, 6298, 6344, 6686, 6880, 6896,
3632, 6044, 6794, 6666, 12420, 7242, 7108, 6668, 6562, 6258, 6690, 6764, 6390, 6610, 6424, 6404, 6668, 6854, 6626, 6520, 7004, 6782, 6816, 6750, 6446, 6832, 6716, 6744, 6910, 6982, 6092, 6710, 6406, 6394, 6620, 6696, 6670, 6798, 6494, 6264, 6680, 6656, 6778, 6408, 6658, 6686, 6576, 6572, 6872, 6652, 6780, 6418, 6656, 6706,
2286, 5534, 6172, 6602, 6792, 19190, 8582, 7058, 6688, 6312, 6542, 6568, 6598, 6612, 6746, 6604, 6902, 6578, 6610, 6296, 6580, 6544, 6556, 7190, 6636, 6602, 6966, 6488, 6530, 6448, 6808, 6488, 6774, 6390, 6182, 6190, 6244, 6306, 6454, 7306, 7036, 7206, 6928, 6938, 6674, 7028, 6828, 6640, 6648, 6820, 6408, 6444, 6558, 6514

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
      while ( rand_sampler > weight )
      {
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
