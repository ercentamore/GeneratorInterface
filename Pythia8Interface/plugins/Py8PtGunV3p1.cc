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
3390, 4096, 4470, 4822, 4588, 4790, 4936, 4980, 4818, 4744, 4894, 4764, 4990, 4780, 4764, 4844, 5016, 4704, 4920, 4720, 4970, 4808, 4892, 4888, 5004, 4830, 4990, 4904, 4834, 4938, 4956, 4862, 4832, 4804, 4884, 4950, 4858, 4766, 4762, 4782, 4904, 4826, 4744, 4896, 4798, 4884, 4814, 4800, 4852, 4886, 4926, 4754, 4738, 4752,
1278, 2380, 3504, 4718, 4770, 4878, 4772, 4864, 4858, 4738, 4790, 4794, 4700, 4798, 4794, 4706, 4906, 4794, 4726, 4830, 4820, 4880, 4872, 4856, 4818, 4728, 4910, 4962, 4824, 4986, 5060, 4594, 4682, 4828, 4704, 4934, 4910, 4798, 4842, 4906, 4838, 4914, 5010, 4800, 4778, 4876, 4648, 4864, 4880, 4898, 4954, 4854, 4868, 4842,
210, 372, 810, 2342, 4110, 4592, 4830, 4820, 4876, 4832, 4722, 4896, 4794, 4928, 4834, 4816, 4812, 4822, 4772, 4846, 4836, 4756, 4790, 5026, 4896, 4854, 4862, 4958, 4784, 5034, 4654, 4832, 4830, 4924, 4922, 4810, 4846, 5026, 4912, 4914, 4782, 4634, 4888, 4858, 4824, 4670, 4734, 4726, 4804, 4800, 4782, 4820, 4898, 4846,
1324, 3908, 4920, 4880, 4888, 4806, 4984, 4878, 4866, 4786, 4738, 4780, 4836, 4832, 5028, 4748, 4686, 4866, 4620, 4764, 4790, 4840, 4708, 4864, 4684, 4856, 4704, 4902, 4974, 4920, 4950, 4914, 4758, 4808, 4860, 4732, 4788, 4780, 4838, 4950, 4844, 4904, 4828, 4716, 4682, 4928, 4980, 4758, 4732, 4736, 4814, 4888, 4778, 4874,
2612, 4416, 4530, 4910, 4844, 4976, 4968, 4806, 4842, 4896, 4822, 4816, 4898, 4822, 5028, 5054, 4754, 4796, 4696, 4906, 4790, 4864, 4738, 4766, 4758, 4702, 5122, 4850, 4848, 4826, 4746, 4834, 4872, 4702, 4806, 5002, 4682, 4922, 4886, 4908, 4892, 4888, 4808, 4744, 4738, 4754, 4974, 4896, 4900, 4862, 4910, 4926, 4814, 4792,
1660, 3806, 4540, 4746, 5064, 4790, 4900, 4772, 4818, 4792, 4802, 4796, 4856, 4936, 4952, 4912, 4810, 4948, 4620, 4912, 4864, 4620, 4876, 4872, 4864, 4746, 4818, 4754, 4978, 4784, 4778, 4866, 4888, 4992, 4298, 4342, 4380, 4616, 4820, 4686, 4750, 4764, 4706, 4770, 4882, 4924, 4922, 4870, 4782, 4926, 4778, 4800, 4980, 4864

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
