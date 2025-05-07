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

    // 5) Lookup the 4-D inverse-PDF weight
    double lookup_invpdf4d(
        double m1, double m2,
        double pt1, double pt2,
        const std::vector<double>& mBins,
        const std::vector<double>& ptBins,
        const std::vector<double>& invpdf4d
    ) {
        const size_t nMass = mBins.size();
        const size_t nPt   = ptBins.size();

        size_t ib1 = find_bin(m1,    mBins);
        size_t ip1 = find_bin(pt1,   ptBins);
        size_t ib2 = find_bin(m2,    mBins);
        size_t ip2 = find_bin(pt2,   ptBins);

        // flatten (ib1,ip1,ib2,ip2) → single index
        size_t idx = ((ib1 * nPt + ip1) * nMass + ib2) * nPt + ip2;
        return invpdf4d[idx];
    }

    // 6) Load or initialize a flat 4-D histogram of ints
    std::vector<int> load_or_init_csv4d(
        const std::string& filename,
        size_t n1, size_t m1,
        size_t n2, size_t m2
    ) {
        size_t Ntot = n1 * m1 * n2 * m2;
        std::vector<int> occ(Ntot, 0);

        std::ifstream ifs(filename);
        if (!ifs) {
            // file not found: leave zeros
            return occ;
        }
        std::string line;
        size_t idx = 0;
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
            double  fMinEta;
            double  fMaxEta;
            double  fMinPt ;
            double  fMaxPt ;
            double  fPtRes ;
            double  fMinMass ;
            double  fMaxMass ;
            double  fMassRes ;
            bool    fAddAntiParticle;
            std::vector<double> m0_bins_;
            std::vector<double> pT_bins_;
            std::vector<double> invpdf4d_;
    };

    // implementation
    //
    Py8PtGunV4::Py8PtGunV4( edm::ParameterSet const& ps )
        : Py8GunBase(ps)
    {
        // ParameterSet defpset ;
        auto p  = ps.getParameter<edm::ParameterSet>("PGunParameters");
        fMinEta        = p.getParameter<double>("MinEta");
        fMaxEta        = p.getParameter<double>("MaxEta");
        fMinPt         = p.getParameter<double>("MinPt");
        fMaxPt         = p.getParameter<double>("MaxPt");
        fPtRes         = p.getParameter<double>("PtRes");
        fMinMass       = p.getParameter<double>("MinMass");
        fMaxMass       = p.getParameter<double>("MaxMass");
        fMassRes       = p.getParameter<double>("MassRes");
        fAddAntiParticle = p.getParameter<bool>("AddAntiParticle");
        m0_bins_ = arange(fMinMass, fMaxMass + fMassRes, fMassRes);
        pT_bins_ = arange(fMinPt,   fMaxPt   + fPtRes,  fPtRes);
        size_t nMass = m0_bins_.size(), nPt = pT_bins_.size();
        auto occ4d  = load_or_init_csv4d("occ4d.csv", nMass, nPt, nMass, nPt);
        invpdf4d_   = get_inverse_pdf(occ4d);
    }

    // Modified event generator: sample BOTH daughters at once
    bool Py8PtGunV4::generatePartonsAndHadronize()
    {
        fMasterGen->event.reset();

        // Draw (mass1, pT1, mass2, pT2) with 4-D acceptance/rejection
        double mass1, mass2, pt1, pt2, u, w;
        do {
            u   = rand() / double(RAND_MAX);
            pt1 = (fMaxPt - fMinPt) * randomEngine().flat() + fMinPt;
            pt2 = (fMaxPt - fMinPt) * randomEngine().flat() + fMinPt;
            mass1 = (fMaxMass - fMinMass) * randomEngine().flat() + fMinMass;
            mass2 = (fMaxMass - fMinMass) * randomEngine().flat() + fMinMass;

            w = lookup_invpdf4d(
                mass1, mass2,
                pt1,   pt2,
                m0_bins_, pT_bins_,
                invpdf4d_
              );
        } while (u > w);

        // Build four-vector for daughter 1
        double phi1 = (fMaxPhi - fMinPhi)*randomEngine().flat() + fMinPhi;
        double eta1 = (fMaxEta - fMinEta)*randomEngine().flat() + fMinEta;
        double theta1 = 2.0 * std::atan(std::exp(-eta1));
        double pp1  = pt1 / std::sin(theta1);
        double E1   = std::sqrt(pp1*pp1 + mass1*mass1);
        double px1  = pt1 * std::cos(phi1);
        double py1  = pt1 * std::sin(phi1);
        double pz1  = pp1 * std::cos(theta1);

        // Build four-vector for daughter 2
        double phi2 = (fMaxPhi - fMinPhi)*randomEngine().flat() + fMinPhi;
        double eta2 = (fMaxEta - fMinEta)*randomEngine().flat() + fMinEta;
        double theta2 = 2.0 * std::atan(std::exp(-eta2));
        double pp2  = pt2 / std::sin(theta2);
        double E2   = std::sqrt(pp2*pp2 + mass2*mass2);
        double px2  = pt2 * std::cos(phi2);
        double py2  = pt2 * std::sin(phi2);
        double pz2  = pp2 * std::cos(theta2);

        // Append daughter 1
        int pid1 = fPartIDs[0];
        if (!(fMasterGen->particleData).isParticle(pid1)) pid1 = std::abs(pid1);
        fMasterGen->event.append(pid1, 23, 101, 0, px1, py1, pz1, E1, mass1);
        if (fAddAntiParticle) {
            int apid1 = -pid1;
            if (!(fMasterGen->particleData).isParticle(apid1)) apid1 = pid1;
            fMasterGen->event.append(apid1, 23, 0, 101, -px1, -py1, -pz1, E1, mass1);
        }

        // Append daughter 2
        int pid2 = fPartIDs[1];
        if (!(fMasterGen->particleData).isParticle(pid2)) pid2 = std::abs(pid2);
        fMasterGen->event.append(pid2, 23, 101, 0, px2, py2, pz2, E2, mass2);
        if (fAddAntiParticle) {
            int apid2 = -pid2;
            if (!(fMasterGen->particleData).isParticle(apid2)) apid2 = pid2;
            fMasterGen->event.append(apid2, 23, 0, 101, -px2, -py2, -pz2, E2, mass2);
        }

        // Finalize the event
        if (!fMasterGen->next()) return false;
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
