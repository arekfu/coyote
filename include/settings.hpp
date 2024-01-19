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

#ifndef __SETTINGS__
#define __SETTINGS__

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

namespace settings
{
    enum class CollisionMode {
        ANALOG,
        BRANCHING,
        BRANCHLESS,
        PROMPT_BRANCHLESS,
        ALT_BRANCHLESS };
    enum class ControlMode {
        ANALOG,
        ROULETTE,
        WGT_COMBING,
        WOR };
    enum class DecayMode {
        ANALOG,
        FORCED_DECAY };
    enum class BoundaryType {
        REFLECTIVE,
        PERIODIC,
        VACUUM };
    enum class SourceType {
        NONEQUILIBRIUM,
        EQUILIBRIUM};
    enum class ControlStrat{
        REGULAR,
        SINGULAR,
        SOURCE };
    enum class SimulationMode {
        TARAPIATAPIOCO,
        STATIC };

    extern int x_bin;
    extern int t_bin;
    extern double x_size;
    extern double t_size;
    extern double dx, dt;
    extern double IFn, IFp, IR;
    extern double roulette_value;
    extern int splitting_value;
    extern bool weighted_sites;
    extern BoundaryType boundary_type;
    extern CollisionMode collision_mode;
    extern CollisionMode static_collision;
    extern ControlMode control_mode;
    extern ControlMode static_control_mode;
    extern DecayMode decay_mode;
    extern SourceType source_type;
    extern ControlStrat control_strat;
    extern SimulationMode simulation_mode;
    extern std::vector<int> control_steps;
    extern bool collapsed_precu;
    extern bool control_on_collision;
    extern int n_particles;
    extern int m_replica;
    extern size_t groups, families;
    extern int regions;
    extern size_t pi_seed;
    extern int passive_gen, active_gen, tot_gen;
    extern double p;
    extern int decoupling_gen;
    extern std::vector<double> v;

    extern bool reactivity_change;
    extern std::vector<size_t> reactivity_time;

    extern bool col_correlations;
    extern bool trace_correlations;
    extern bool photo_correlations;
    extern bool col_flux;
    extern bool trace_flux;
    extern bool col_res_time;
    extern bool trace_res_time;
    extern bool entropy;
    extern bool source_entropy;
    extern bool pair_distance;
    extern bool family_trees;
    extern bool feynman_moment_gen;
    extern bool feynman_moment_rep;
    extern bool emission_density;
    extern bool collision_density;
    extern bool kin_source;
    
    extern bool pi_entropy;
    extern bool pi_flux_col;
    extern bool pi_flux_trace;
    extern bool pi_pair_distance;
    extern bool pi_emission_density_gen;
    extern bool pi_emission_density_rep;
    extern bool pi_collision_density_gen;
    extern bool pi_collision_density_rep;

    extern std::string input;
    extern std::string output;
}

#endif
