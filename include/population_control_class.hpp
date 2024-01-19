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

#ifndef __population_control__
#define __population_control__

#include <random>
#include <vector>
#include <cassert>
#include "particle_class.hpp"
#include "vector_class.hpp"
#include "pcg_random.hpp"

class PopulationControl
{
    public:
        PopulationControl();
        
        /* Apply population control at the end of each time step, if enable */
        virtual void applyPopulationControl([[maybe_unused]] std::vector<Particle>& _buffer) = 0;
        virtual ~PopulationControl();
    protected:
        pcg64 rng;
        double t;
        std::uniform_real_distribution<double> R;
};

/* Analog case */
class NoControl: public PopulationControl
{
    public:
        NoControl();
        virtual void applyPopulationControl([[maybe_unused]] std::vector<Particle>& _buffer);
        virtual ~NoControl();
};

/* Apply splitting and roulette (not on collision, on end of time step) */
class SplittingRoulette: public PopulationControl
{
    public:
        SplittingRoulette();
        virtual void applyPopulationControl(std::vector<Particle>& _buffer);
        void splitting(std::vector<Particle>& _buffer);
        void russianRoulette(std::vector<Particle>& _buffer);
        virtual ~SplittingRoulette();
    private:
        double russian_threshold, splitting_threshold, base_weight;
};

/* weight-based combing, bringing the number of particles to _N while keeping
 * the weight a constant */
class WeightCombing: public PopulationControl
{
    public:
        WeightCombing();
        virtual void applyPopulationControl(std::vector<Particle>& _buffer);
        virtual ~WeightCombing();
    private:
        std::vector<Particle> n_buffer;
        std::vector<Particle> p_buffer;
};

#endif
