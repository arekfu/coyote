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

#include "parser.hpp"
#include "error.hpp"
#include "settings.hpp"
#include "material_class.hpp"
#include "yaml-cpp/node/node.h"
#include "yaml-cpp/node/parse.h"
#include <string>
#include <iostream>
#include <memory>
#include <filesystem>

using namespace std;

vector<shared_ptr<Material>> materials;
shared_ptr<Geometry> geometry;
vector<shared_ptr<Source>> sources;

void parseInput()
{
    YAML::Node read_input = YAML::LoadFile(settings::input);
    if (read_input["settings"])
        makeSettings(read_input["settings"]);
    else
        fatal_error("no settings field", __FILE__, __LINE__);

    if (read_input["materials"] and read_input["materials"].IsSequence())
        makeMaterials(read_input["materials"]);
    else
        fatal_error("no materials field", __FILE__, __LINE__);

    if (read_input["regions"] and read_input["regions"].IsSequence())
        makeGeometry(read_input["regions"]);
    else
        fatal_error("no regions field", __FILE__, __LINE__);

    if (read_input["tallies"])
        makeTallies(read_input["tallies"]);
    else
        fatal_error("no tallies field", __FILE__, __LINE__);

    if (read_input["sources"] and read_input["sources"].IsSequence())
        makeSources(read_input["sources"]);
    else
        cout << "No sources provided: switching to equilibrium source." << endl;
    std::ofstream fout(settings::output + "config.yaml");
    fout << read_input;
}

void makeSettings(const YAML::Node& setnode)
{
    cout << "\nReading simulation configuration." << endl;

    if (setnode["simulation-mode"])
    {
        string mode = setnode["simulation-mode"].as<string>();
        if (mode == "tarapiaTapioco")
            settings::simulation_mode = settings::SimulationMode::TARAPIATAPIOCO;
        else if (mode == "static")
            settings::simulation_mode = settings::SimulationMode::STATIC;
        else
            fatal_error("invalid simulation mode", __FILE__, __LINE__);
    }
    else
        fatal_error("no simulation-mode field", __FILE__, __LINE__);

    if (setnode["x-bin"])
        settings::x_bin = setnode["x-bin"].as<int>();
    else
        fatal_error("no x-bin field", __FILE__, __LINE__);
    
    if (setnode["t-bin"])
        settings::t_bin = setnode["t-bin"].as<int>();
    else
        fatal_error("no t-bin field", __FILE__, __LINE__);

    if (setnode["x-size"])
        settings::x_size = setnode["x-size"].as<double>();
    else
        fatal_error("no x-size field", __FILE__, __LINE__);
    
    if (setnode["t-size"])
        settings::t_size = setnode["t-size"].as<double>();
    else
        fatal_error("no t-size field", __FILE__, __LINE__);

    settings::dx = 2*settings::x_size/settings::x_bin;
    settings::dt = settings::t_size/settings::t_bin;

    if (setnode["neutron-importance"])
        settings::IFn = setnode["neutron-importance"].as<double>();

    if (setnode["precu-importance"])
        settings::IFp = setnode["precu-importance"].as<double>();
    settings::IR = settings::IFn / settings::IFp;

    if (setnode["boundaries"])
    {
        string bound = setnode["boundaries"].as<string>();
        if (bound == "reflective")
            settings::boundary_type = settings::BoundaryType::REFLECTIVE;
        else if (bound == "periodic")
            settings::boundary_type = settings::BoundaryType::PERIODIC;
        else if (bound == "vacuum")
            settings::boundary_type = settings::BoundaryType::VACUUM;
        else
            fatal_error("invalid boundary condition", __FILE__, __LINE__);
    }
    else
        fatal_error("no boundaries field", __FILE__, __LINE__);

    if (setnode["collision"])
    {
        string col = setnode["collision"].as<string>();
        if (col == "analog")
            settings::collision_mode = settings::CollisionMode::ANALOG;
        else if (col == "branching")
            settings::collision_mode = settings::CollisionMode::BRANCHING;
        else if (col == "branchless")
            settings::collision_mode = settings::CollisionMode::BRANCHLESS;
        else if (col == "prompt-branchless")
            settings::collision_mode = 
                settings::CollisionMode::PROMPT_BRANCHLESS;
        else if (col == "alt-branchless")
            settings::collision_mode = 
                settings::CollisionMode::ALT_BRANCHLESS;
        else
            fatal_error("invalid collision mode", __FILE__, __LINE__);
    }
    else 
        fatal_error("no collision field", __FILE__, __LINE__);

    if (setnode["control"])
    {
        string control = setnode["control"].as<string>();
        if (control == "analog")
            settings::control_mode = settings::ControlMode::ANALOG;
        else if (control == "roulette")
            settings::control_mode = settings::ControlMode::ROULETTE;
        else if (control == "combing")
            settings::control_mode = settings::ControlMode::WGT_COMBING;
        else
            fatal_error("invalid control mode", __FILE__, __LINE__);
    }
    else
        fatal_error("no control field", __FILE__, __LINE__);

    if (setnode["static-control"]) {
        string control = setnode["static-control"].as<string>();
        if (control == "analog")
            settings::static_control_mode = settings::ControlMode::ANALOG;
        else if (control == "combing")
            settings::static_control_mode = settings::ControlMode::WGT_COMBING;
        else if (control == "wor")
            settings::static_control_mode = settings::ControlMode::WOR;
        else
            fatal_error("invalid static control field", __FILE__, __LINE__);
    }
    if (setnode["weighted-sites"]) 
        settings::weighted_sites = setnode["weighted-sites"].as<bool>();

    if (setnode["control-strategy"]) {
        string strat = setnode["control-strategy"].as<string>();
        if (strat == "regular")
            settings::control_strat = settings::ControlStrat::REGULAR;
        else if (strat == "singular") {
            settings::control_strat = settings::ControlStrat::SINGULAR;
            if (setnode["control-steps"] and setnode["control-steps"].IsSequence())
                settings::control_steps = setnode["control-steps"].as<vector<int>>();
            else
                fatal_error("no control-steps field provided", __FILE__, __LINE__);
        }
        else if (strat == "source")
            settings::control_strat = settings::ControlStrat::SOURCE;
    }
    else {
        cout << "Invalid control-strategy field: switching to regular control" 
             << endl;
    }

    if (setnode["decay"])
    {
        string dc = setnode["decay"].as<string>();
        if (dc == "analog")
            settings::decay_mode = settings::DecayMode::ANALOG;
        else if (dc == "forced")
            settings::decay_mode = settings::DecayMode::FORCED_DECAY;
        else
            fatal_error("invalid decay mode", __FILE__, __LINE__);
    }
    else
        fatal_error("no decay field", __FILE__, __LINE__);

    if (setnode["control-on-collision"])
        settings::control_on_collision = 
            setnode["control-on-collision"].as<bool>();
    else
        fatal_error("no control-on-collision field", __FILE__, __LINE__);
    

    if (setnode["collapsed-precu"]) {
        if (setnode["collapsed-precu"].as<bool>()) {
            if (settings::decay_mode == settings::DecayMode::FORCED_DECAY) 
                settings::collapsed_precu = setnode["collapsed-precu"].as<bool>();
            else
                fatal_error("cannot use collapsed precu with analogue decay", __FILE__, __LINE__);
        }
        else 
            cout << "Using analogue precursor particle" << endl;
        }
    else 
        cout << "Using analogue precursor particle" << endl;

    if (setnode["roulette-value"])
        settings::roulette_value = setnode["roulette-value"].as<double>();
    if (setnode["splitting-value"])
        settings::splitting_value = setnode["splitting-value"].as<int>();

    if (setnode["nparticles"])
        settings::n_particles = setnode["nparticles"].as<int>();
    else
        fatal_error("no nparticles field", __FILE__, __LINE__);
    
    if (setnode["groups"])
        settings::groups = setnode["groups"].as<size_t>();
    else
        fatal_error("no groups field", __FILE__, __LINE__);

    if (setnode["families"])
        settings::families = setnode["families"].as<size_t>();
    else
        fatal_error("no families field", __FILE__, __LINE__);

    if (setnode["regions"])
        settings::regions = setnode["regions"].as<int>();
    else
        fatal_error("no regions field", __FILE__, __LINE__);

    if (setnode["speed"] and setnode["speed"].IsSequence())
        settings::v = setnode["speed"].as<vector<double>>();
    else
        fatal_error("no speed field", __FILE__, __LINE__);

    if (setnode["static-collision"]) {
        string static_col = setnode["static-collision"].as<string>();
        if (static_col == "branching")
            settings::static_collision = settings::CollisionMode::BRANCHING;
        else if (static_col == "analog")
            settings::static_collision = settings::CollisionMode::ANALOG;
        else if (static_col == "branchless")
            settings::static_collision = settings::CollisionMode::BRANCHLESS;
        else
            fatal_error("invalid static collision field", __FILE__, __LINE__);
    }
    else
        fatal_error("no static-collision field", __FILE__, __LINE__);

    if (setnode["passive-gen"])
        settings::passive_gen = setnode["passive-gen"].as<int>();
    else fatal_error("no passive-gen field", __FILE__, __LINE__);

    if (setnode["active-gen"])
        settings::active_gen = setnode["active-gen"].as<int>();
    else
        fatal_error("no active-gen field", __FILE__, __LINE__);
    settings::tot_gen = settings::active_gen + settings::passive_gen;
    if (setnode["sampling-prob"])
        settings::p = setnode["sampling-prob"].as<double>();
    if (setnode["decoupling-gen"])
        settings::decoupling_gen = setnode["decoupling-gen"].as<int>();
    if (setnode["reactivity-change"])
        settings::reactivity_change = setnode["reactivity-change"].as<int>();
    if (setnode["reactivity-time"] and setnode["reactivity-time"].IsSequence())
        settings::reactivity_time = setnode["reactivity-time"].as<vector<size_t>>();
}

void makeMaterials(const YAML::Node& matsnode)
{
    cout << "Constructing materials." << endl;
    for (size_t i = 0; i < matsnode.size(); i++)
        materials.push_back(make_shared<Material>(matsnode[i]));
}

void makeGeometry(const YAML::Node& regionsnode)
{
    cout << "Constructing geometry." << endl;

    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry = 
            make_shared<GeometryRef>(regionsnode);
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry = 
            make_shared<GeometryPer>(regionsnode); 
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry = 
            make_shared<GeometryLeak>(regionsnode);
    else
        fatal_error("Invalid choice of boundary condition", __FILE__, __LINE__);
}       

void makeTallies(const YAML::Node& talliesnode)
{
    cout << "Reading enabled tallies." << endl;
    
    if (talliesnode["collision-flux"])
        settings::col_flux = talliesnode["collision-flux"].as<bool>();
    if (talliesnode["trace-flux"])
        settings::trace_flux = talliesnode["trace-flux"].as<bool>();
    if (talliesnode["collision-correlations"])
        settings::col_correlations = talliesnode["collision-correlations"].as<bool>();
    if (talliesnode["trace-correlations"])
        settings::trace_correlations = talliesnode["trace-correlations"].as<bool>();
    if (talliesnode["collision-residence"])
        settings::col_res_time = talliesnode["collision-residence"].as<bool>();
    if (talliesnode["trace-residence"])
        settings::trace_res_time = talliesnode["trace-residence"].as<bool>();
    if (talliesnode["position-correlations"])
        settings::photo_correlations = talliesnode["position-correlations"].as<bool>();
    if (talliesnode["entropy"])
        settings::entropy = talliesnode["entropy"].as<bool>();
    if (talliesnode["pair-distance"])
        settings::pair_distance = talliesnode["pair-distance"].as<bool>();
    if (talliesnode["family-trees"])
        settings::family_trees = talliesnode["family-trees"].as<bool>();
    if (talliesnode["emission-density"])
        settings::emission_density = talliesnode["emission-density"].as<bool>();
    if (talliesnode["collision-density"])
        settings::collision_density = talliesnode["collision-density"].as<bool>();
    if (talliesnode["kin-source"])
        settings::kin_source = talliesnode["kin-source"].as<bool>();
    if (talliesnode["pi-entropy"])
        settings::pi_entropy = talliesnode["pi-entropy"].as<bool>();
    if (talliesnode["pi-flux-col"])
        settings::pi_flux_col = talliesnode["pi-flux-col"].as<bool>();
    if (talliesnode["pi-flux-trace"])
        settings::pi_flux_trace = talliesnode["pi-flux-trace"].as<bool>();
    if (talliesnode["pi-pair-distance"])
        settings::pi_pair_distance = talliesnode["pi-pair-distance"].as<bool>();
    if (talliesnode["pi-feynman-moment-gen"])
        settings::feynman_moment_gen = talliesnode["pi-feynman-moment-gen"].as<bool>();
    if (talliesnode["pi-feynman-moment-rep"])
        settings::feynman_moment_rep = talliesnode["pi-feynman-moment-rep"].as<bool>();
    if (talliesnode["pi-emission-density-gen"])
        settings::pi_emission_density_gen = 
            talliesnode["pi-emission-density-gen"].as<bool>();
    if (talliesnode["pi-emission-density-rep"])
        settings::pi_emission_density_rep = 
            talliesnode["pi-emission-density-rep"].as<bool>();
    if (talliesnode["pi-collision-density-gen"])
        settings::pi_collision_density_gen = 
            talliesnode["pi-collision-density-gen"].as<bool>();
    if (talliesnode["pi-collision-density-rep"])
        settings::pi_collision_density_rep = 
            talliesnode["pi-collision-density-rep"].as<bool>();
    
    for (size_t i = 0; i < max(settings::groups, settings::families); i++)
        filesystem::create_directories(settings::output + "energy_group" + to_string(i) + "/");
}

void makeSources(const YAML::Node& sourcenode) {
    cout << "Setting up sources." << endl;

    settings::source_type = settings::SourceType::NONEQUILIBRIUM;
    int tmp_number = 0;

    for (size_t i = 0; i < sourcenode.size(); i++) {
        if (sourcenode[i]["type"].as<string>() == "delta")
            sources.push_back(make_shared<DeltaSource>(sourcenode[i]));
        else if (sourcenode[i]["type"].as<string>() == "uniform")
            sources.push_back(make_shared<UniformSource>(sourcenode[i]));
        else
            fatal_error("invalid choice of non-equilibrium source", __FILE__, __LINE__);
        tmp_number += sources.back()->number();
    }
    
    if (tmp_number != settings::n_particles)
        fatal_error("please check n_particles field", __FILE__, __LINE__);
}

const YAML::Node parseGeometry()
{
    YAML::Node read_input = YAML::LoadFile(settings::input);
    if (read_input["regions"] and read_input["regions"].IsSequence())
        return read_input["regions"];
    else
        fatal_error("no regions field provided", __FILE__, __LINE__);
    
    exit(1);
}
