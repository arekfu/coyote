/* Copyright CEA (2024)
 * Contributors: Théophile Bonnet, Davide Mancusi
 * tlfab2@cam.ac.uk, davide.mancusi@cea.fr
 * 
 * This software is a computer program whose purpose is to investigate
 * correlations in neutron transport calculations in very simplified
 * conditions.
 * 
 * This software is governed by the CeCILL  license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef __static_scores__
#define __static_scores__

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <optional>
#include <utility>
#include <tuple>
#include "ressources.hpp"
#include "tally.hpp"
#include "particle_class.hpp"
#include "vector_class.hpp"
#include "ndarray.hpp"
#include "geometry_class.hpp"
#include "settings.hpp"
#include "random_class.hpp"

class StaticScores
{
    public:
        StaticScores();
        StaticScores(const StaticScores& scores) = default;
        
        void scoreGeneration(const double& k, std::vector<Particle>& _buffer,
                             const int& _gen);
        void scoreFissionSitesInner(const int& _gen);
        void scoreReplica();
        void scoreCol(const Particle& obj);
        void scoreTrace(const Particle& obj, const double& t_dist);
        void scoreEntropy(const std::vector<Particle>& n_buffer, const size_t& _gen);
        void scoreFamilyTrees(const std::vector<Particle>& n_buffer, const size_t& _gen);
        void scorePairDistance(std::vector<Particle>& _buffer, const size_t& _gen);
        void scoreLeakage(const double& wgt);
        void scoreFeynmanMomentGenTmp();
        void scoreFissionSitesInnerTmp(const std::vector<Particle>& obj);
        void scoreFissionSitesRepTmp(const std::vector<Particle>& obj);
        void computeFeynmanMomentRep();
        void scoreEmissionDensityTmp(const Particle& obj);
        void saveWeights(const std::vector<Particle>& _buffer, const size_t _gen);
        void nextReplica();
        void normalize_replicas();
        void normalize_generations();
        void reinitialize();
        void save();

        Tally1D const & FluxCol() const;
        Tally1D const & FluxTrace() const;
        Tally1D const & Entropy() const;
        Tally1D const & FamilyTrees() const;
        Tally1D const & N2() const;
        Tally1D const & PairDistance() const;
        Tally1D const & FissionSitesInner() const;
        std::pair<Tally1D, Tally1D> const & FissionSitesOuter() const;
        Tally1D const & FissionSitesRep() const;
        Tally1D const & FeynmanMomentGen() const;
        Tally1D const & FeynmanMomentRep() const;
        Tally1D const & EmissionDensityInner() const;
        Tally1D const & TotEmissionDensityInner() const;
        std::pair<Tally1D, Tally1D> const & EmissionDensityOuter() const;
        std::pair<Tally1D, Tally1D> const & TotEmissionDensityOuter() const;
        Tally1D const & EmissionDensityRep() const;
        Tally1D const & TotEmissionDensityRep() const;
        Tally1D const & CollisionDensityInner() const;
        Tally1D const & TotCollisionDensityInner() const;
        std::pair<Tally1D, Tally1D> const & CollisionDensityOuter() const;
        std::pair<Tally1D, Tally1D> const & TotCollisionDensityOuter() const;
        Tally1D const & CollisionDensityRep() const;
        Tally1D const & TotCollisionDensityRep() const;
        Tally0D const & Kstep() const;
        Tally0D const & Leakage() const; 
        Tally0D const & MigArea() const;
        void scoreMigAreaTmp(const Particle& obj);

        inline StaticScores& operator+=(const StaticScores& rhs);

    private:
        std::unique_ptr<Geometry> geometry_ptr;
        RNG rng;
        int m;
        Tally1D mFluxCol;
        Tally1D mFluxTrace;
        Tally1D mEntropy, mFamilyTrees;
        Tally1D mPairDistance;
        Tally1D mN2;

        Tally1D mFissionSitesInner; //only through inner generations!
        Tally1D mEmissionDensityInner; //only through inner generations!
        Tally1D mTotEmissionDensityInner; //only through inner generations!
        Tally1D mCollisionDensityInner; //only through inner generations!
        Tally1D mTotCollisionDensityInner; //only through inner generations!
        std::pair<Tally1D, Tally1D> mFissionSitesOuter;
        std::pair<Tally1D, Tally1D> mEmissionDensityOuter;
        std::pair<Tally1D, Tally1D> mTotEmissionDensityOuter;
        std::pair<Tally1D, Tally1D> mCollisionDensityOuter;
        std::pair<Tally1D, Tally1D> mTotCollisionDensityOuter;

        Tally1D mFissionSitesRep;
        Tally1D mFeynmanMomentGen;
        Tally1D mFeynmanMomentRep;

        Tally1D mEmissionDensityRep;
        Tally1D mTotEmissionDensityRep;

        Tally1D mCollisionDensityRep;
        Tally1D mTotCollisionDensityRep;

        Tally0D mKstep;
        Tally0D mLeakage;
        Tally0D mMigArea;
};

inline StaticScores& StaticScores::operator+=(const StaticScores& rhs) {
    if (settings::pi_flux_col)
        mFluxCol += rhs.FluxCol();
    if (settings::pi_flux_trace)
        mFluxTrace += rhs.FluxTrace();
    if (settings::pi_entropy) {
        mEntropy += rhs.Entropy();
        mFamilyTrees += rhs.FamilyTrees();
    }
    if (settings::pi_pair_distance) {
        mPairDistance += rhs.PairDistance();
        mN2 += rhs.N2();
    }
    if (settings::feynman_moment_gen) {
        mFeynmanMomentGen += rhs.FeynmanMomentGen();
        mFissionSitesOuter.first += rhs.FissionSitesOuter().first;
        mFissionSitesOuter.second += rhs.FissionSitesOuter().second;
    }
    if (settings::feynman_moment_rep) {
        mFeynmanMomentRep += rhs.FeynmanMomentRep();
        mFissionSitesRep += rhs.FissionSitesRep();
    }
    if (settings::pi_emission_density_gen) {
        mEmissionDensityOuter.first += rhs.EmissionDensityOuter().first;
        mEmissionDensityOuter.second += rhs.EmissionDensityOuter().second;
        mTotEmissionDensityOuter.first += rhs.TotEmissionDensityOuter().first;
        mTotEmissionDensityOuter.second += rhs.TotEmissionDensityOuter().second;
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep += rhs.EmissionDensityRep();
        mTotEmissionDensityRep += rhs.TotEmissionDensityRep();
    }
    if (settings::pi_collision_density_gen) {
        mCollisionDensityOuter.first += rhs.CollisionDensityOuter().first;
        mCollisionDensityOuter.second += rhs.CollisionDensityOuter().second;
        mTotCollisionDensityOuter.first += rhs.TotCollisionDensityOuter().first;
        mTotCollisionDensityOuter.second += rhs.TotCollisionDensityOuter().second;
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep += rhs.CollisionDensityRep();
        mTotCollisionDensityRep += rhs.TotCollisionDensityRep();
    }
    mKstep += rhs.Kstep();
    mLeakage += rhs.Leakage();
    mMigArea += rhs.MigArea();

    return *this;
}

#endif
