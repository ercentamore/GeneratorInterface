#include "GeneratorInterface/Pythia8Interface/interface/CustomHook.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Pythia8/Pythia.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

using namespace Pythia8;

// Helper utilities

static int sum(const std::vector<int>& dist) {
    return std::accumulate(dist.begin(), dist.end(), 0);
}

static double max_element(const std::vector<double>& dist) {
    return *std::max_element(dist.begin(), dist.end());
}

static std::vector<double> get_inverse_pdf(const std::vector<int>& dist) {
    std::vector<double> invpdf(dist.size());
    double tot = static_cast<double>(sum(dist));
    for (size_t i = 0; i < dist.size(); ++i) {
        invpdf[i] = (dist[i] > 0) ? tot / static_cast<double>(dist[i]) : 1.0;
    }
    double mx = max_element(invpdf);
    for (double &v : invpdf) v /= mx;
    return invpdf;
}

static size_t find_bin(double x, const std::vector<double> &bins) {
    auto it = std::upper_bound(bins.begin(), bins.end(), x);
    if (it == bins.begin()) return 0;
    size_t idx = std::distance(bins.begin(), it) - 1;
    return std::min(idx, bins.size() - 1);
}

static double lookup_invpdf1d(double m1,
                            const std::vector<double> &mBins,
                            const std::vector<double> &invpdf1d) {
    size_t ib1 = find_bin(m1, mBins);
    return invpdf1d[ib1];
}

static std::vector<int> load_or_init_csv1d(const std::string &filename,
                                           size_t n) {
  std::vector<int> occ(n, 100); // flat default
  if (filename.empty()) return occ;

  std::ifstream ifs(filename);
  if (!ifs) return occ;

  std::string line, cell;
  size_t idx = 0;
  std::getline(ifs, line);
  std::stringstream ss(line);
    while (idx < occ.size() && std::getline(ss, cell, ',')){
        occ[idx++] = std::stoi(cell);
  }
  return occ;
}

static std::vector<double> arange(double start, double stop, double step){
    std::vector<double> v;
    for (double x = start; x < stop - 0.5 * step; x += step)
        v.push_back(x);
    return v;
}

static inline double sqr(double x){
    return x * x;
}

// The actual Userhook
class HtoAAHook : public UserHooks {
public:
    explicit HtoAAHook(const edm::ParameterSet &ps){
                minM_(ps.getParameter<double>("minMass")),
                maxM_(ps.getParameter<double>("maxMass")),
                res_(ps.getParameter<double>("massRes")),
                csv_(ps.getUntrackedParameter<std::string>("csvFile", "")),
                idM_(ps.getParameter<int>("motherID")),
                idD1_(ps.getParameter<int>("daughter1ID")),
                idD2_(ps.getParameter<int>("daughter2ID")),
                rng_(0.0, 1.0),
            // Build mass grid and corresponding inverse-PDF table.
            bins_ = arange(minM_, maxM_ + res_, res_);
            auto occ1d = load_or_init_csv1d(csv_, bins_.size());
            invpdf1d_ = get_inverse_pdf(occ1d);
            // Seed RNG with nondeterministic value.
            engine_.seed(std::random_device{}());
    }

    // Operate at the process level after all standard decays.
    // bool canVetoProcessLevel() override { return true; }
    bool canVetoProcessLevel() override { return true; }

    void fixupTwoBodyDecays(Event& event, int iMom, int iDau1, int iDau2){
        double MM = event[iMom].m();
        double m1 = event[iDau1].m();
        double m2 = event[iDau2].m();

        const double p2cm = (sqr(MM) - sqr(m1 + m2)) * (sqr(MM) - sqr(m1 - m2));
        if (p2cm <= 0.0) return;
        const double p = 0.5 * std::sqrt(p2cm) / MM;
        const double cosTh = 2.0 * rndmPtr->flat() - 1.0;
        const double sinTh = std::sqrt(std::max(0.0, 1.0 - cosTh * cosTh));
        const double phi = 2.0 * M_PI * rndmPtr->flat();

        Vec4 p1RF(p * sinTh * std::cos(phi), p * sinTh * std::sin(phi),
        p * cosTh, std::sqrt(p * p + m1 * m1));
        Vec4 p2RF(-p1RF.px(), -p1RF.py(), -p1RF.pz(),
        std::sqrt(p * p + m2 * m2));
        p1RF.bst(event[iMom].p());
        p2RF.bst(event[iMom].p());
        event[iDau1].p( p1RF );
        event[iDau2].p( p2RF );
        return;
    }

    // bool doVetoProcessLevel(Event &event) override {
    bool doVetoProcessLevel(Event &event) override {
        // Loop over all Higgs bosons, irrespective of status code. We will
        // overwrite any pre-existing decay chain.
        for (int iH = 0; iH < event.size(); ++iH) {
            if (event[iH].id() != idM_) continue;

            // Require exactly two daughters to repurpose; otherwise skip.
            int d1 = event[iH].daughter1();
            int d2 = event[iH].daughter2();
            if (d1 <= 0 || d2 <= 0 || d2 < d1 || d2 != d1 + 1) {
                continue;
            }

            // Cache Higgs four-vector and mass before we modify anything.
            const double MM = event[iH].m();

            // Draw (m1, m2) according to user-supplied 2-D inverse PDF grid.
            double m1 = 0.0, m2 = 0.0, weight = 0.0, u = 0.0;
            do {
                m1 = m2 = rndmPtr->flat() * (maxM_ - minM_) + minM_;
                if (m1 + m2 >= MM) continue; // kinematic veto
                weight = lookup_invpdf1d(m1, bins_, invpdf1d_);
                u = rndmPtr->flat();
            } while (u > weight);
            
            if( m1 + m2 < MM ) {
                event[d1].m(m1);
                event[d2].m(m2);
                fixupTwoBodyDecays(event, iH, d1, d2);
            }

            int iDau1 = event[d1].daughter1();
            int iDau2 = event[d1].daughter2();
            if( event[iDau1].m() + event[iDau2].m() > m1 ) return true;
            fixupTwoBodyDecays(event, d1, iDau1, iDau2);
            iDau1 = event[d2].daughter1();
            iDau2 = event[d2].daughter2();
            if( event[iDau1].m() + event[iDau2].m() > m2 ) return true;
            fixupTwoBodyDecays(event, d2, iDau1, iDau2);
        
        }

        // Never veto - we have modified in-place.
        return false;
    }

private:
    // User-configurable parameters.
    double minM_, maxM_, res_;
    std::string csv_;

    // Pre-computed tables.
    std::vector<double> bins_, invpdf1d_;

    // PDG IDs.
    int idM_, idD1_, idD2_;

    // Random number machinery.
    std::mt19937 engine_;
    std::uniform_real_distribution<> rng_;
};

REGISTER_USERHOOK(HtoAAHook);