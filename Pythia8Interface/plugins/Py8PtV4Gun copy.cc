#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"

#include <numeric>
#include <algorithm>
#include <cmath>
#include <memory>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace gen {
    // Declare helper functions
    // 1) Sum of integers
    int sum(const std::vector<int>& dist) {
        return std::accumulate(dist.begin(), dist.end(), 0);
    }

    // 2) Max over doubles
    double max_element(const std::vector<double>& dist) {
        double mx = 0.0;
        for (double v : dist) {
            if (v > mx) mx = v;
        }
        return mx;
    }

    // 3) Build the inverse-PDF and normalize to [0,1]
    std::vector<double> get_inverse_pdf(const std::vector<int>& dist) {
        std::vector<double> invpdf(dist.size());
        double total = sum(dist);
        for (size_t i = 0; i < dist.size(); ++i) {
            if (dist[i] > 0) invpdf[i] = total / double(dist[i]);
            else            invpdf[i] = 1.0;
        }
        double mx = max_element(invpdf);
        for (double& v : invpdf) v /= mx;
        return invpdf;
    }

    // 4) Find bin index for a value x in a sorted list of bin-edges
    static size_t find_bin(double x, const std::vector<double>& bins) {
        auto it = std::upper_bound(bins.begin(), bins.end(), x);
        if (it == bins.begin()) return 0;
        size_t idx = std::distance(bins.begin(), it) - 1;
        return std::min(idx, bins.size() - 1);
    }

    // 5) Lookup the 3-D inverse-PDF weight
    double lookup_invpdf2d(
        double m1, double m2,
        const std::vector<double>& mBins,
        const std::vector<double>& invpdf2d
    ) {
        const size_t nMass = mBins.size();

        size_t ib1 = find_bin(m1,    mBins);
        size_t ib2 = find_bin(m2,    mBins);

        // flatten (ib1,ip1,ib2,ip2) → single index
        size_t idx = (ib1 * nMass + ib2);
        return invpdf2d[idx];
    }

    // 6) Load or initialize a flat 2-D histogram of ints
    std::vector<int> load_or_init_csv2d(
        const std::string& filename,
        size_t nMass1, size_t nMass2
    ) {
        // total number of cells
        size_t Ntot = nMass1 * nMass2;
        // default occupancy (you used 100 before)
        std::vector<int> occ(Ntot, 100);
    
        std::ifstream ifs(filename);
        if (!ifs) {
            // file not found → return flat default
            return occ;
        }
    
        std::string line;
        size_t idx = 0;
        // The file should have Ntot comma-separated ints,
        // we just read them in row-major order.
        while (idx < Ntot && std::getline(ifs, line)) {
            std::stringstream ss(line);
            std::string cell;
            while (idx < Ntot && std::getline(ss, cell, ',')) {
                occ[idx++] = std::stoi(cell);
            }
        }
        return occ;
    }

    // 7) Equivalent of numpy.arange(start, stop, step)
    std::vector<double> arange(double start, double stop, double step) {
        size_t n = static_cast<size_t>(std::ceil((stop - start) / step));
        std::vector<double> v;
        v.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            double val = start + i * step;
            if (val >= stop) break;
            v.push_back(val);
        }
        return v;
    }

    static double sqr(double x) { return x * x; }

    class Py8PtGunV4 : public Py8GunBase {
        public:

            Py8PtGunV4( edm::ParameterSet const& );
            ~Py8PtGunV4() {}

            bool generatePartonsAndHadronize() override;
            const char* classname() const override;

        private:

            // PtGun particle(s) characteristics
            double fMinEta;
            double fMaxEta;
            double fMinMass;
            double fMaxMass;
            double fMassRes;
            double fParentMass;
            bool fAddAntiParticle;
            std::vector<int> fDaughterIDs;
            std::vector<double> m0_bins_;
            std::vector<double> pT_bins_;
            std::vector<double> invpdf2d_;
    };

    // implementation
    //
    Py8PtGunV4::Py8PtGunV4( edm::ParameterSet const& ps )
        : Py8GunBase(ps)
    {
        // ParameterSet defpset ;
        auto p  = ps.getParameter<edm::ParameterSet>("PGunParameters");
        fMinEta          = p.getParameter<double>("MinEta");
        fMaxEta          = p.getParameter<double>("MaxEta");
        fMinMass         = p.getParameter<double>("MinMass");
        fMaxMass         = p.getParameter<double>("MaxMass");
        fMassRes         = p.getParameter<double>("MassRes");
        fParentMass      = p.getParameter<double>("ParentMass");
        fAddAntiParticle = p.getParameter<bool>("AddAntiParticle");
        fDaughterIDs     = p.getParameter<std::vector<int>>("DaughterIDs");
        m0_bins_ = arange(fMinMass, fMaxMass + fMassRes, fMassRes);
        size_t nMass = m0_bins_.size();
        auto occ2d  = load_or_init_csv2d("occ2d.csv", nMass, nMass);
        invpdf2d_   = get_inverse_pdf(occ2d);

    }

    bool Py8PtGunV4::generatePartonsAndHadronize()
    {
        if (!fMasterGen->next()) return false;      // <-- generates the Higgs

        // Find the Higgs entry
        int iH = -1;
        for (int i = 0; i < fMasterGen->event.size(); ++i) {
            if (fMasterGen->event[i].id() == 25 && fMasterGen->event[i].status() == 22) {
                iH = i; break;
            };
        };
        if (iH < 0) return false;                   // should never happen
      
        const Pythia8::Particle& H = fMasterGen->event[iH];
    
        // sample (m1, m2) from your inverse-PDF grid
        double m1, m2, u, w;
        do {
          u     = randomEngine().flat();
          m1 = (fMaxMass - fMinMass)*randomEngine().flat() + fMinMass;
          m2 = (fMaxMass - fMinMass)*randomEngine().flat() + fMinMass;
    
          w = lookup_invpdf2d(
            m1, m2,
            m0_bins_,
            invpdf2d_
          );
        } while(u > w);
    
        // after the loop:
        const double p2cm = (M * M - sqr(m1 + m2)) * (M * M - sqr(m1 - m2));
        if (p2cm <= 0.0)
          return false;  // numerical safety
    
        const double p = 0.5 * std::sqrt(p2cm) / M;
        const double costh = 2.0 * rng_(engine_) - 1.0;
        const double sinth = std::sqrt(1.0 - costh * costh);
        const double phi = 2.0 * M_PI * rng_(engine_);
    
        Vec4 p1LAB(p * sinth * std::cos(phi), p * sinth * std::sin(phi), p * costh,
                   std::sqrt(p * p + m1 * m1));
        Vec4 p2LAB = -p1LAB;
        p2LAB.e(std::sqrt(p * p + m2 * m2));
    
        // Boost to lab frame
        const Particle &par = event[iRes];
        const double bx = par.px() / par.e();
        const double by = par.py() / par.e();
        const double bz = par.pz() / par.e();
        const double b2 = bx * bx + by * by + bz * bz;
        const double gamma = 1.0 / std::sqrt(1.0 - b2);
        const double gp = (b2 > 0.0) ? (gamma - 1.0) / b2 : 0.0;
    
        auto boost = [&](Vec4 &v) {
          const double bp = bx * v.px() + by * v.py() + bz * v.pz();
          v.px(v.px() + gp * bp * bx + gamma * bx * v.e());
          v.py(v.py() + gp * bp * by + gamma * by * v.e());
          v.pz(v.pz() + gp * bp * bz + gamma * bz * v.e());
          v.e(gamma * (v.e() + bp));
        };
    
        boost(p1LAB);
        boost(p2LAB);
    
        auto pH = Pythia8::Vec4(pxH, pyH, pzH, EH);
        // append parent and daughters
        fMasterGen->event.append(
            /* id        */ fPartIDs[0],
            /* status    */ 22,
            /* m1,m2     */ 0, 0,
            /* d1,d2     */ 2, 3,
            /* col1,col2 */ 0, 0,
            /* p         */ pH,
            /* m0        */ mH
        );

        auto pA = Pythia8::Vec4(px1, py1, pz1, E1);
        fMasterGen->event.append(
          fDaughterIDs[0], 23, 1, 0, 0, 0, 0, 0,
          pA, mass1
        );
        auto pB = Pythia8::Vec4(px2, py2, pz2, E2);
        fMasterGen->event.append(
          fDaughterIDs[1], 23, 1, 0, 0, 0, 0, 0,
          pB, mass2
        );
    
        // let Pythia finish the job
        if(!fMasterGen->next()) return false;
        event().reset(new HepMC::GenEvent);
        return toHepMC.fill_next_event(fMasterGen->event, event().get());
    }

    const char* Py8PtGunV4::classname() const {
        return "Py8PtGunV4";
    }


    typedef edm::GeneratorFilter<gen::Py8PtGunV4, gen::ExternalDecayDriver> Pythia8PtGunV4;

} // end namespace

using gen::Pythia8PtGunV4;
DEFINE_FWK_MODULE(Pythia8PtGunV4);
