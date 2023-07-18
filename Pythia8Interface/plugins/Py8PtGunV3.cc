#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"

#include <numeric>

namespace gen {

class Py8PtGunV3 : public Py8GunBase {

   public:

      Py8PtGunV3( edm::ParameterSet const& );
      ~Py8PtGunV3() {}

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
Py8PtGunV3::Py8PtGunV3( edm::ParameterSet const& ps )
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

int sum(std::vector <int> dist) {
    return std::accumulate(dist.begin(), dist.end(), 0);
}

double max_element(std::vector <double> dist) {
    double max = 0;
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        double el = dist[i];
        if (max < el){max = el;}
    }
    return max;
}

std::vector <double> get_inverse_pdf(std::vector <int> dist) {
    std::vector <double> invpdf(dist.size());
    double sum_hist = sum(dist);
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        if (dist[i] != 0 ) {
            invpdf[i] = sum_hist / dist[i];
            //std::cout << "Bin " << i << " -> " << invpdf[i] << std::endl;
        }
        else {invpdf[i] = 1;}
    }
    double max_invpdf = max_element(invpdf);
    for (int i = 0; i < s; i++) {
        invpdf[i] = invpdf[i] / max_invpdf;
    }
    return invpdf;
}

double lookup_mass_invpdf(double mgen, std::vector <double> m_bins, std::vector <double> m_invpdf) {
    int im = 0;
    int s1 = m_bins.size();
    int s2 = m_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        im = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (mgen <= m_bins[ib]) { break; }
    }
    return m_invpdf[im];
}

double lookup_pt_invpdf(double pTgen, std::vector <int> pT_bins, std::vector <double> pT_invpdf) {
    int ipt = 0;
    int s1 = pT_bins.size();
    int s2 = pT_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        ipt = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (pTgen <= pT_bins[ib]) { break; }
    }
    return pT_invpdf[ipt];
}

double lookup_invpdf(double Mgen, std::vector <double> M_bins, double pTgen, std::vector <int> pT_bins, std::vector <double> invpdf) {
    unsigned int ibin = 0;
    unsigned int m1  = M_bins.size();
    unsigned int pt1 = pT_bins.size();
    unsigned int inv = invpdf.size();
    bool found_mass = false;
    bool found_end  = false;
    for (unsigned int ibx = 0; ibx < m1; ibx++) {
        if (found_mass || found_end) { break; }
        for (unsigned int iby = 0; iby < pt1; iby++) {
            ibin = (ibx*pt1)+ iby;
            if ( ((ibx*pt1) + iby + 1) >  (inv - 1) ) {
                found_end = true;
                break;
            }
            if ( (Mgen  <= M_bins[ibx]) && (pTgen <= pT_bins[iby]) ) {
                found_mass = true;
                break;
            }
        }
    }
    return invpdf[ibin];
}

double get_rand_el(std::vector <int> dist) {
    int randomIndex = rand() % dist.size();
      return dist[randomIndex];
}

std::vector <int> pT_bins   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180};
std::vector <double> m_bins = {4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.4, 8.8, 9.2, 9.6, 10.0, 10.4, 10.8, 11.2, 11.6, 12.0, 12.4, 12.8, 13.2, 13.6, 14.0, 14.4, 14.8, 15.2, 15.6, 16.0};
std::vector <int> occ = {
//  35,     40,     45,     50,     55,     60,     65,     70,     75,     80,     85,     90,     95,    100,    105,    110,   115,     120,    125,    130,    135,    140,    145,    150,    155,    160,    165,    170,    175,    180
481450, 618290, 702250, 762020, 806340, 823970, 831230, 840120, 850110, 861730, 855320, 859300, 861390, 867200, 864970, 879600, 868390, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, 875880, // -> 4.0
546460, 698700, 801370, 853890, 898670, 921090, 935890, 950630, 970280, 966420, 974310, 977630, 982480, 980070, 988750, 984460, 993560, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, 993610, // ->  4.4
543940, 687860, 782700, 848770, 885130, 912710, 930350, 940750, 943770, 959790, 967690, 962740, 973890, 963690, 971000, 981920, 970810, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, 975820, // ->  4.8
531550, 677070, 768110, 825250, 866860, 902180, 912200, 928310, 937360, 945350, 947720, 955360, 954110, 956270, 958990, 958720, 959160, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, 960000, // ->  5.2
514470, 662660, 755980, 818090, 853570, 882080, 897260, 914130, 926590, 935440, 939010, 937210, 943750, 942670, 948760, 961920, 951950, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, 955180, // ->  5.6
496650, 652700, 744450, 809750, 838160, 872260, 879320, 902950, 915380, 925460, 918770, 935470, 939320, 942280, 937260, 931860, 937310, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, 943040, // ->  6.0
478400, 639740, 733160, 800200, 836550, 860560, 877940, 896280, 907870, 913130, 916780, 922190, 923250, 938350, 919400, 925930, 928060, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, 930860, // ->  6.4
419760, 613150, 724320, 785680, 826620, 845570, 871750, 883330, 895980, 900980, 896330, 920020, 910170, 916150, 919630, 923500, 926300, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, 925630, // ->  6.8
264370, 576570, 707680, 773920, 812100, 838970, 857380, 876860, 883270, 889260, 904720, 903560, 908730, 911500, 916450, 914050, 910560, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, 910580, // ->  7.2
106660, 453370, 674800, 761410, 799370, 838010, 856310, 871340, 875430, 884150, 887870, 894790, 903000, 899220, 911240, 907720, 907370, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, 907390, // ->  7.6
 22850, 290860, 579840, 737040, 791480, 825730, 843270, 853780, 880730, 873370, 881720, 891790, 893440, 899060, 899200, 901230, 902230, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, 909090, // ->  8.0
  6700, 135750, 447860, 658580, 778940, 818180, 837630, 855120, 866980, 874490, 884620, 879120, 893710, 892010, 892960, 891500, 895950, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, 892820, // ->  8.4
  4280,  36080, 318270, 546180, 718010, 801580, 835400, 847420, 854710, 870300, 868480, 881640, 881390, 882970, 882300, 892100, 898010, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, 884270, // ->  8.8
  2730,   5580, 173820, 441630, 624580, 741300, 806960, 838050, 851980, 859150, 861930, 877460, 876330, 871250, 884690, 883570, 884840, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, 886240, // ->  9.2
  1600,   3710,  59380, 336300, 528540, 660230, 763020, 819650, 848860, 857150, 862440, 863150, 868800, 879090, 880570, 874570, 884470, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, 876020, // ->  9.6
   970,   2110,   6560, 217860, 441860, 582800, 692050, 766570, 823170, 848610, 858290, 860010, 877650, 872340, 869680, 875430, 872600, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, 875990, // -> 10.0
   660,   1540,   2810,  95200, 359770, 513980, 623320, 713370, 772100, 822980, 837990, 850330, 852150, 862630, 863000, 871730, 864510, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, 870800, // -> 10.4
   430,    960,   1920,  16570, 262340, 445940, 562970, 650970, 717990, 779510, 820290, 848220, 855310, 856930, 860800, 867870, 861340, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, 861050, // -> 10.8
   220,    620,   1120,   2020, 136670, 376520, 504210, 596230, 669280, 726800, 775670, 810580, 838760, 852920, 851890, 864280, 860380, 865950, 860380, 860380, 860380, 860380, 860380, 860380, 860380, 860380, 860380, 860380, 860380, 860380, // -> 11.2
   130,    420,    790,   1340,  39220, 295130, 444510, 547470, 624450, 678310, 729610, 775790, 805590, 838280, 846620, 852810, 854960, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, 859890, // -> 11.6
    87,    220,    490,    950,   1940, 189470, 389040, 497840, 574930, 639810, 686980, 733440, 762850, 807440, 829890, 840990, 848920, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 12.0
//Extra bins
    47,    135,    314,    527,   1053,  59280, 277885, 435254, 574930, 639810, 662445, 643331, 675667, 748612, 772983, 840990, 848920, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 12.4
    27,    100,    207,    426,    673,   7580, 235091, 386181, 461586, 488997, 611412, 643331, 675667, 682863, 772983, 798940, 807686, 835579, 835579, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 12.8
    17,     46,    141,    247,    488,    857, 101706, 322173, 409021, 488997, 564305, 577322, 618998, 682863, 684066, 768905, 776155, 776155, 776155, 776155, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 13.2
     8,     40,     95,    151,    316,    524,  28900, 265988, 403272, 483513, 499532, 567892, 612459, 682863, 628345, 764099, 771304, 771304, 771304, 771304, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 13.6
     5,     21,     53,    109,    227,    325,    669, 140106, 331816, 413134, 460276, 517599, 552521, 682863, 615304, 621131, 686412, 730071, 730071, 730071, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 14.0
     2,     10,     41,     78,    144,    245,    361,  60595, 300293, 413134, 460276, 517599, 552521, 682863, 615304, 621131, 686412, 730071, 730071, 730071, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 14.4
     1,      7,     22,     45,    100,    209,    311,   5254, 209044, 413134, 460276, 517599, 552521, 682863, 615304, 621131, 686412, 730071, 730071, 730071, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, 849220, // -> 14.8
     1,      2,      8,     23,     40,     60,     72,    350,  61385, 182811, 222083, 252329, 233440, 317531, 346280, 346280, 346280, 346280, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, // -> 15.2
     1,      2,      8,     23,     40,     60,     72,    350,  61385, 182811, 222083, 252329, 233440, 317531, 346280, 346280, 346280, 346280, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, // -> 15.6
     1,      2,      8,     23,     40,     60,     72,    350,  61385, 182811, 222083, 252329, 233440, 317531, 346280, 346280, 346280, 346280, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, // -> 16.0
     1,      2,      8,     23,     40,     60,     72,    350,  61385, 182811, 222083, 252329, 233440, 317531, 346280, 346280, 346280, 346280, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549, 354549 // -> inf
};
std::vector <double> invpdf = get_inverse_pdf(occ);


bool Py8PtGunV3::generatePartonsAndHadronize()
{

   //std::cout << " ! ! UPDATED PYTHIA VERSION ! ! " << std::endl;
   fMasterGen->event.reset();

   for ( size_t i=0; i<fPartIDs.size(); i++ )
   {

      int particleID = fPartIDs[i]; // this is PDG - need to convert to Py8 ???

      double rand_sampler = rand() / double(RAND_MAX);
      double pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
      double mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
      double weight       = lookup_invpdf(mass, m_bins, pt, pT_bins, invpdf);
      while ( rand_sampler > weight ) {
         rand_sampler = rand() / double(RAND_MAX);
         pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
         mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
         weight       = lookup_invpdf(mass, m_bins, pt, pT_bins, invpdf);
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

const char* Py8PtGunV3::classname() const
{
   return "Py8PtGunV3";
}

typedef edm::GeneratorFilter<gen::Py8PtGunV3, gen::ExternalDecayDriver> Pythia8PtGunV3;

} // end namespace

using gen::Pythia8PtGunV3;
DEFINE_FWK_MODULE(Pythia8PtGunV3);
