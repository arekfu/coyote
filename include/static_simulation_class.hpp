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

#ifndef __new_pi__
#define __new_pi__

#include <algorithm>
#include <memory>
#include <random>
#include <vector>
#include <utility>
#include <tuple>
#include "pcg_random.hpp"
#include "ndarray.hpp"
#include "particle_class.hpp"
#include "geometry_class.hpp"
#include "random_class.hpp"
#include "population_control_class.hpp"
#include "source_class.hpp"
#include "static_collision_kernel_class.hpp"
#include "static_scores_class.hpp"
#include "static_population_control.hpp"

class StaticSimulation
{
    public:
        StaticSimulation();
        void source();
        //void run(const int& m);
        void antani(const int& k);
        //void antaniLast();
        void update(const int& k);

        bool transportNeutron(const int& n, const int& _gen);
        //void applyPopulationControl();

        //void updateTarapiaTapiocoSource(const int& n, const double& prob);
        //void computeTarapiaTapiocoSource(const int& cycles);

        void reinitialize();
        void normalize();
        void makeScores();
        void pi_combing();
        void takePhoto(const int& rep, const std::string& path) const;
        void computeTarapiaTapiocoBuffer(size_t cycles);
        void updateTarapiaTapiocoBuffer(int n, double prob);
        void antaniLast();
        std::vector<Particle> const & CurrentBuffer() const;
        std::vector<Particle> const & TarapiaTapiocoBuffer() const;

        void computeFeynmanMomentRep();
        StaticScores const & Scores() const;
        ~StaticSimulation();

    protected:
        RNG rng;
        std::unique_ptr<Geometry> geometry_ptr;
        std::unique_ptr<StaticControl> population_control;
        std::unique_ptr<StaticCollision> collision_kernel;
//        std::vector<std::shared_ptr<Source>> sources_ptr;
        std::vector<Particle> mCurrentBuffer;
        std::vector<Particle> mTarapiaTapiocoBuffer;
        std::vector<Particle> mCombedBuffer;
        std::vector<std::tuple<Vector,double, int>> mFissionSites;
        StaticScores mScores;
        bool mSourceSimulation;
        double mOldWgt;
};

#endif
