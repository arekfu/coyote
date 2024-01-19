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

#include "settings.hpp"

namespace settings 
{
    int x_bin = 100;
    int t_bin = 25;
    double x_size = 1;
    double t_size = 1;
    double dx = x_size*2/x_bin;
    double dt = t_size/t_bin;
    double IFn = 1, IFp = 1, IR = 1;
    double roulette_value = 0.8;
    int splitting_value = 2;
    bool weighted_sites = true;

    BoundaryType boundary_type = BoundaryType::REFLECTIVE;
    CollisionMode collision_mode = CollisionMode::ANALOG;
    CollisionMode static_collision = CollisionMode::BRANCHLESS;
    ControlMode control_mode = ControlMode::ANALOG;
    ControlMode static_control_mode = ControlMode::WGT_COMBING;
    DecayMode decay_mode = DecayMode::ANALOG;
    SourceType  source_type = SourceType::EQUILIBRIUM;
    ControlStrat control_strat = ControlStrat::REGULAR;
    SimulationMode simulation_mode = SimulationMode::STATIC;

    std::vector<int> control_steps = std::vector<int>();
    bool control_on_collision = false;
    bool collapsed_precu = false;
    int n_particles = 10000;
    int m_replica = 100;
    std::size_t groups = 3, families = 6;
    int regions = 3;
    size_t pi_seed = 2125529455;
    int passive_gen = 100, active_gen = 100, tot_gen = 200;
    double p = 1;
    int decoupling_gen = 1;
    std::vector<double> v(groups, 0);

    bool reactivity_change = false;
    std::vector<size_t> reactivity_time = std::vector<size_t>(1,0);

    bool col_flux = true;
    bool trace_flux = true;
    bool col_correlations = false;
    bool trace_correlations = false;
    bool photo_correlations = false;
    bool col_res_time = false;
    bool trace_res_time = false;
    bool entropy = false;
    bool source_entropy = false;
    bool pair_distance = false;
    bool family_trees = false;
    bool emission_density = false;
    bool collision_density = false;

    bool pi_entropy = false;
    bool pi_flux_col = false;
    bool pi_flux_trace = false;
    bool pi_pair_distance = false;
    bool feynman_moment_gen = false;
    bool feynman_moment_rep = false;
    bool pi_emission_density_gen = false;
    bool pi_emission_density_rep = false;
    bool pi_collision_density_gen = false;
    bool pi_collision_density_rep = false;
    bool kin_source = false;

    std::string input = "/home/catB/tb264612/Documents/thèse/code/coyote/input/";
    std::string output = "/home/catB/tb264612/Documents/thèse/code/coyote/data/saved/";
}
