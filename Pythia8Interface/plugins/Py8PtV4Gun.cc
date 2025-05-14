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
    double lookup_invpdf3d(
        double m1, double m2,
        double pt1,
        const std::vector<double>& mBins,
        const std::vector<double>& ptBins,
        const std::vector<double>& invpdf3d
    ) {
        const size_t nMass = mBins.size();
        const size_t nPt   = ptBins.size();

        size_t ib1 = find_bin(m1,    mBins);
        size_t ip1 = find_bin(pt1,   ptBins);
        size_t ib2 = find_bin(m2,    mBins);

        // flatten (ib1,ip1,ib2,ip2) → single index
        size_t idx = ((ib1 * nPt + ip1) * nMass + ib2) * nPt;
        return invpdf3d[idx];
    }

    // 6) Load or initialize a flat 3-D histogram of ints
    std::vector<int> load_or_init_csv3d(
        const std::string& filename,
        size_t nMass1, size_t nMass2,
        size_t nPt1
    ) {
        // total number of cells
        size_t Ntot = nMass1 * nMass2 * nPt1;
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
            double fMinPt;
            double fMaxPt;
            double fPtRes;
            double fMinMass;
            double fMaxMass;
            double fMassRes;
            double fParentMass;
            bool fAddAntiParticle;
            std::vector<int> fDaughterIDs;
            std::vector<double> m0_bins_;
            std::vector<double> pT_bins_;
            std::vector<double> invpdf3d_;
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
        fMinPt           = p.getParameter<double>("MinPt");
        fMaxPt           = p.getParameter<double>("MaxPt");
        fPtRes           = p.getParameter<double>("PtRes");
        fMinMass         = p.getParameter<double>("MinMass");
        fMaxMass         = p.getParameter<double>("MaxMass");
        fMassRes         = p.getParameter<double>("MassRes");
        fParentMass      = p.getParameter<double>("ParentMass");
        fParentMass      = p.getParameter<double>("ParentMass");
        fAddAntiParticle = p.getParameter<bool>("AddAntiParticle");
        fDaughterIDs     = p.getParameter<std::vector<int>>("DaughterIDs");
        m0_bins_ = arange(fMinMass, fMaxMass + fMassRes, fMassRes);
        pT_bins_ = arange(fMinPt,   fMaxPt   + fPtRes,  fPtRes);
        size_t nMass = m0_bins_.size(), nPt = pT_bins_.size();
        auto occ3d  = load_or_init_csv3d("occ3d.csv", nMass, nPt, nMass);
        invpdf3d_   = get_inverse_pdf(occ3d);
    }

    bool Py8PtGunV4::generatePartonsAndHadronize()
    {
        // clear the old event
        fMasterGen->event.reset();     
    
        // sample Higgs production kinematics via the two PDF x's
        // this gives you rapidity y_H and therefore a boost along z,
        // plus a transverse kick pT_H vs. phi_H.
        const double yMin  = -2.4, yMax = +2.4;
        double yH   = yMin  + (yMax - yMin)*randomEngine().flat();
        double phiH = 2*M_PI * randomEngine().flat();
    
        // choose pT_H from Pythia’s built-in spectrum, or fall back to flat:
        double ptH  = (fMaxPt - fMinPt)*randomEngine().flat() + fMinPt;
        double mH   = fParentMass;                      // e.g. 125 GeV
        double mT   = std::sqrt(mH*mH + ptH*ptH);
        double EH   = mT * std::cosh(yH);
        double pzH  = mT * std::sinh(yH);
        double pxH  = ptH * std::cos(phiH);
        double pyH  = ptH * std::sin(phiH);
    
        // compute beta = p/E and gamma for the boost
        double bx = pxH/EH, by = pyH/EH, bz = pzH/EH;
        double b2 = bx*bx + by*by + bz*bz;
        double gamma = 1.0/std::sqrt(1.0 - b2);
        // helper to boost a four-vector (E,px,py,pz) -> lab
        auto boostToLab = [&](double &E, double &px, double &py, double &pz){
          double bp = bx*px + by*py + bz*pz;
          double gamma2 = (b2>0 ? (gamma-1)/b2 : 0);
          px += gamma2*bp*bx + gamma*bx*E;
          py += gamma2*bp*by + gamma*by*E;
          pz += gamma2*bp*bz + gamma*bz*E;
          E  = gamma*(E + bp);
        };
    
        // sample (m1, m2, pT1) from your inverse-PDF grid
        double mass1, mass2, pt1, u, w;
        do {
          u     = randomEngine().flat();
          pt1   = (fMaxPt - fMinPt)*randomEngine().flat() + fMinPt;
          mass1 = (fMaxMass - fMinMass)*randomEngine().flat() + fMinMass;
          mass2 = (fMaxMass - fMinMass)*randomEngine().flat() + fMinMass;
    
          w = lookup_invpdf3d(
            mass1, mass2, pt1,
            m0_bins_, pT_bins_,
            invpdf3d_
          );
        } while(u > w);
    
        // after the loop:
        const double M2   = mH*mH;
        const double sum  = mass1 + mass2;
        const double diff = mass1 - mass2;
        const double arg  = (M2 - sum*sum)*(M2 - diff*diff);
        if(arg < 0) return false;        // no real two-body solution
        const double pstar = std::sqrt(arg)/(2.0*mH);

        // if your sampled pt1 is kinematically too big, bail out
        if(pt1 > pstar) return false;

        // pick a random phi around the z-axis
        double phi1 = 2*M_PI * randomEngine().flat();

        // build px1,py1 in rest frame
        double px1 = pt1 * std::cos(phi1);
        double py1 = pt1 * std::sin(phi1);

        // the longitudinal piece is fixed by |p|=pstar
        double pz_abs = std::sqrt(pstar*pstar - pt1*pt1);
        // choose + or - with 50/50 chance
        double pz1 = (randomEngine().flat()<0.5 ? +1 : -1) * pz_abs;

        double E1 = std::sqrt(pstar*pstar + mass1*mass1);
    
        // daughter 2 exactly opposite:
        double px2 = -px1, py2 = -py1, pz2 = -pz1;
        double E2  = std::sqrt(pstar*pstar + mass2*mass2);
    
        // boost both daughters into the lab frame
        boostToLab(E1, px1, py1, pz1);
        boostToLab(E2, px2, py2, pz2);
    
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
