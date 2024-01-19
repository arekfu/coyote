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

#ifndef __scores__
#define __scores__

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <optional>
#include "ressources.hpp"
#include "tally.hpp"
#include "particle_class.hpp"
#include "vector_class.hpp"
#include "ndarray.hpp"
#include "geometry_class.hpp"
#include "random_class.hpp"

/* scores for tarapia tapioco simulations */
class Scores
{
    public:
    /* constructor*/
    Scores();
    Scores(const Scores& scores) = default;

    /* compute photo score contribution of obj to time bin time */
    void scoreTime(size_t time, std::vector<Particle>& buffer);

    void scoreKinSource(const std::vector<Particle>& _buffer);
    void scoreCombedKinSource(const std::vector<Particle>& _buffer);

    /* temp scores of Tallies corresponding to collision estimators */
    void scoreCol(const Particle& obj, size_t current_reactivity);

    /* temp scores of Tallies correponding to trace estimators */
    void scoreTrace(const Particle& obj, const double& t_dist);

    /* temp scores of time dependent entropy */
    void scoreEntropy(size_t time, const std::vector<Particle>& buffer);

    void scorePairDistance(size_t, std::vector<Particle>& buffer);

    void scoreEmissionDensityTmp(const Particle& obj, const size_t t_bin);

    /* compute photo related estimators */
    void computeBasicScores(); // compute ct, S, empty cell, h from temp_mat

    /* normalize all scores */
    void normalizeBasicScores(); // normalize wrt dxg and M
    /* reinitialize temporary scores from Tallies */
    void reinitializeTempScores();

    /* function corresponding to operator += */
    void addScores(const Scores& other);

    /* m++ */
    void nextReplica();

    /* getter */
    Tally1D const & getPPairDist() const;
    Tally1D const & getNPairDist() const;
    Tally1D const & getN2() const;
    Tally1D const & getP2() const;
    Tally1D const & getFamilyTrees() const;
    Tally2D const & getParticleTotal() const;
    Tally2D const & getNeutronDensity() const;
    Tally2D const & getExpNeutronDensity() const;
    Tally2D const & getPrecuDensity() const;
    Tally2D const & getResTimeCol() const;
    Tally2D const & getResTimeTrace() const;
    Tally2D const & getFluxCol() const;
    Tally2D const & getFluxTrace() const;
    Tally4D const & getNeutronCorrelations() const;
    Tally4D const & getColCorrelations() const;
    Tally4D const & getTraceCorrelations() const;
    Tally2D const & getPEntropy() const;
    Tally2D const & getNEntropy() const;
    Tally1D const & getEntropy() const;
    Tally4D const & getEntropyCov() const;
    NDArray<double> const & getSourceEntropy() const;
    Tally2D const & getEmissionDensity() const;
    Tally2D const & getCollisionDensity() const;
    Tally2D const & getTotEmissionDensity() const;
    Tally2D const & getTotCollisionDensity() const;
    Tally1D const & getCombedKinSourceNeutrons() const;
    Tally1D const & getKinSourceNeutrons() const;
    Tally1D const & getCombedKinSourcePrecu() const;
    Tally1D const & getKinSourcePrecu() const;

    void setSourceEntropy(const NDArray<double> & obj);

    Scores& operator+=(const Scores& rhs);

    /* save scores with path provided in the function itself */
    void saveScores(const std::string& output);

    virtual ~Scores();

    protected:
    std::unique_ptr<Geometry> geometry_ptr;
    int m; //current replica and number of energy groups
    Tally1D P_pair_dist, N_pair_dist;
    Tally1D N2, P2;
    Tally1D family_trees;
    Tally2D P_entropy, N_entropy; 
    Tally1D entropy;
    Tally4D entropy_cov;
    NDArray<double> source_entropy;
    Tally2D particle_total;
    Tally2D neutron_density, precu_density;
    Tally2D res_time_col, res_time_trace;
    Tally2D flux_col, flux_trace;
    Tally4D neutron_correlations;
    Tally4D col_correlations, trace_correlations;
    std::map<int, int> family_trees_map;

    Tally2D mEmissionDensity;
    Tally2D mCollisionDensity;
    Tally2D mTotEmissionDensity;
    Tally2D mTotCollisionDensity;
    Tally1D mCombedKinSourceNeutrons;
    Tally1D mKinSourceNeutrons;
    Tally1D mCombedKinSourcePrecu;
    Tally1D mKinSourcePrecu;
    RNG rng;
};

#endif
