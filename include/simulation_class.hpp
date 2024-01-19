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

#ifndef __simulation__
#define __simulation__

#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include "pcg_random.hpp"
#include "ressources.hpp"
#include "particle_class.hpp"
#include "vector_class.hpp"
#include "scores_class.hpp"
#include "geometry_class.hpp"
#include "population_control_class.hpp"
#include "collision_kernel_class.hpp"
#include "decay_kernel_class.hpp"
#include "source_class.hpp"
#include "static_simulation_class.hpp"
#include "settings.hpp"
#include "random_class.hpp"

class Simulation
{
    public:
    Simulation();

    /* transport neutron of index n in buffer to its next collision */
    bool transportNeutron(const int& n);

    /* apply population control at the end a time step */
    void applyPopulationControl();

    /* getters */
    const Scores& getScores() const;
    const std::vector<Particle>& getBuffer();
    int neutronNum() const; // get number of particles
    int precursorNum() const; // get number of precursors

    /* clear buffer */
    void reinitialize();

    /* normalize wrt to the weight of the partial score */
    void normalize();

    /* arbitrary source with N_d precu, N1 neutrons in group 1 and N2 in group 2*/
    void setSource();

    /* removed dead particles from buffer */
    void cleanBuffer();

    void takePhoto(const int& rep, const std::string& path) const;

    /* advance time until tp */
    void advanceTime(const double& tp);

    /* calls computeBasicScores() of scores */
    void computeScores();
    void scoreEntropy(size_t time);
    void scorePairDist(size_t time);
    
    /* change current reactivity (next set of capture x section in material) */
    void changeReactivity();

    virtual ~Simulation();

    protected:
    RNG rng;
    std::vector<Particle> buffer;
    double t_bin;
    size_t current_reactivity;
    std::unique_ptr<Geometry> geometry_ptr;
    std::unique_ptr<PopulationControl> population_control; 
    std::unique_ptr<CollisionKernel> collision_kernel; 
    std::unique_ptr<DecayKernel> decay_kernel;
    std::vector<std::shared_ptr<Source>> sources_ptr;
    Scores scores; 
    StaticSimulation source_computation;
};

#endif
