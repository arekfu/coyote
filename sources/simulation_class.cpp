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

#include <iostream>
#include <stdexcept>
#include <vector>
#include <random>
#include <cmath>
#include <utility>
#include <memory>
#include "omp.h"
#include "random_class.hpp"
#include "simulation_class.hpp"
#include "population_control_class.hpp"
#include "source_class.hpp"
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

Simulation::Simulation() 
    : rng(1121055555 + omp_get_thread_num()), t_bin(0), current_reactivity(0), scores(), source_computation()
{
    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry_ptr = make_unique<GeometryRef>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry_ptr = make_unique<GeometryPer>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry_ptr = make_unique<GeometryLeak>(parseGeometry()); 
    else
        fatal_error("Invalid choice of boundary condition", __FILE__, __LINE__);


    if (settings::control_mode == settings::ControlMode::ANALOG)
        population_control = make_unique<NoControl>();
    else if (settings::control_mode == settings::ControlMode::ROULETTE)
        population_control = make_unique<SplittingRoulette>();
    else if (settings::control_mode == settings::ControlMode::WGT_COMBING)
        population_control = make_unique<WeightCombing>();
    else
        fatal_error("Invalid choice of population control", __FILE__, __LINE__);

    if (settings::collision_mode == settings::CollisionMode::ANALOG)
        collision_kernel = make_unique<AnalogCollisionKernel>();
    else if (settings::collision_mode == settings::CollisionMode::BRANCHING)
        collision_kernel = make_unique<BasicCollisionKernel>();
    else if (settings::collision_mode == settings::CollisionMode::BRANCHLESS)
        collision_kernel = make_unique<BranchlessCollisionKernel>();
    else if (settings::collision_mode == settings::CollisionMode::PROMPT_BRANCHLESS)
        collision_kernel = make_unique<PromptBranchlessCollisionKernel>();
    else if (settings::collision_mode == settings::CollisionMode::ALT_BRANCHLESS)
        collision_kernel = make_unique<AltBranchless>();
    else
        fatal_error("Invalid choice of collision mode", __FILE__, __LINE__);

    if (settings::decay_mode == settings::DecayMode::ANALOG)
        decay_kernel = make_unique<AnalogDecayKernel>();
    else if (settings::decay_mode == settings::DecayMode::FORCED_DECAY)
        decay_kernel = make_unique<ForcedDecayKernel>();
    else
        fatal_error("Invalid choice of decay mode", __FILE__, __LINE__);

    if (settings::source_type == settings::SourceType::NONEQUILIBRIUM)
        for (const auto& s: sources)
            sources_ptr.emplace_back(s);
    else {
        source_computation.source();
        for (int i = 1; i <= settings::passive_gen; i++) {
            source_computation.antani(i);
            source_computation.update(i);
        }
    }

    buffer.reserve(settings::n_particles);
}
const vector<Particle>& Simulation::getBuffer() {return buffer;}

int Simulation::neutronNum() const
{
    double tmp = 0;
    for (const auto& object: buffer)
        if (not object.is_precu())
            tmp += object.wgt();
    return tmp;
}

int Simulation::precursorNum() const
{
    double tmp = 0;
    for (const auto& object: buffer)
        if (object.is_precu())
            tmp += object.wgt();
    return tmp;
}

void Simulation::reinitialize() 
{
    buffer.clear();
    scores.reinitializeTempScores();
    current_reactivity = 0;
}

void Simulation::normalize() {scores.normalizeBasicScores();}

void Simulation::cleanBuffer()
{
    vector<Particle> tmp;
    tmp.reserve(buffer.size());
    for (const auto& obj: buffer)
        if (obj.is_alive())
            tmp.emplace_back(obj);
    swap(tmp, buffer);
    tmp.clear();
}

void Simulation::takePhoto(const int& rep, const string& path) const
{
    string name = path + "t" + to_string(rep) + ".txt";
    ofstream photo(name);
    if (photo)
    {
        for (const auto& obj: buffer)
            photo << obj.r().x() << " " << obj.r().y() << " " << obj.r().z() << " " << obj.wgt() << endl;
        photo.close();
    }
}

void Simulation::advanceTime(const double& tp)
{
    t_bin = tp*settings::dt; //set the final time for current time step
    for (unsigned int i = 0; i < buffer.size(); i++)
    {
        if (buffer[i].is_precu() and buffer[i].is_alive())
            decay_kernel->decayPrecursor(buffer, i, tp);
        else
        {
            while (buffer[i].t() < t_bin and buffer[i].is_alive())
            {
                if (settings::emission_density)
                    scores.scoreEmissionDensityTmp(buffer[i], tp-1);
                bool collision = transportNeutron(i);
                if (collision) 
                {
                    scores.scoreCol(buffer[i], current_reactivity);
                    collision_kernel->collisionNeutron(buffer,i,current_reactivity);
                }
                if (settings::control_on_collision) {
                    buffer[i].roulette(rng.R());
                    splitting(buffer, i);
                }
            }
        }
    }
    cleanBuffer();
    scores.scoreTime(tp, buffer);
}

void Simulation::scoreEntropy(size_t time) {scores.scoreEntropy(time, buffer);}
void Simulation::scorePairDist(size_t time) {scores.scorePairDistance(time, buffer);}

void Simulation::setSource()
{
    if (settings::source_type == settings::SourceType::NONEQUILIBRIUM)
        for (auto& ptr: sources_ptr) {
            for (int i = 0; i < ptr->number(); i++) {
                buffer.push_back(ptr->sampleSourceParticle());
                buffer.back().setIn_region(geometry_ptr->getRegionIndex(buffer.back()));
                buffer.back().setTree(i+1); 
            }
        }

    else {
        double total_weight = 0, w_factor;
        source_computation.computeTarapiaTapiocoBuffer(settings::decoupling_gen);
        buffer = source_computation.TarapiaTapiocoBuffer();

        for (const auto& obj: buffer)
            total_weight += obj.wgt();

        w_factor = settings::n_particles/total_weight;
        for (auto& obj: buffer)
        {
            //renormalizing total weight to wanted total population
            obj.multiplyWeight(w_factor);
            if (obj.is_precu()) {
                obj.multiplyExpWeight(w_factor/settings::IFp);
                obj.multiplyWeight(1./settings::IFp);//importance sampling adjustment
            }
            else
                obj.multiplyWeight(1./settings::IFn);//importance sampling adjustment
        }
    }
    
    if (settings::kin_source)
        scores.scoreKinSource(buffer);
    applyPopulationControl();
    cleanBuffer();
    if (settings::kin_source)
        scores.scoreCombedKinSource(buffer);
    
    scores.scoreTime(0,buffer);
    if (settings::entropy)
        scores.scoreEntropy(0,buffer);
    if (settings::pair_distance)
        scores.scorePairDistance(0,buffer);
}

void Simulation::computeScores() {
    scores.computeBasicScores();
    scores.nextReplica();
}

bool Simulation::transportNeutron(const int& n)
{
    bool collision = false, end_time = false, leakage = false;
    int g = buffer[n].Eg();

    while (not collision and not end_time and not leakage) //until next collision or next time_bin
    {
        //update tarapiaTapioco variables for n-th neutron
        double t = buffer[n].t();
        const Region& current = geometry_ptr->getRegion(buffer[n].in_region());
        
        //get distance to next border/boundary and corresponding time
        double max_displacement = geometry_ptr->distToFirstSurface(buffer[n]);
        double max_t_time = max_displacement/buffer[n].v();

        //sample transport distance
        double t_dist = -log(1-rng.R())/current.material().E_t(current_reactivity)[g];
        double t_time = t_dist/buffer[n].v();

        //if we go outside the time bin, then reduce the transport distance
        if (t + t_time >= t_bin)
        {
            t_time = t_bin - t;
            t_dist = t_time*buffer[n].v();
            if (t_dist > max_displacement)
            {
                scores.scoreTrace(buffer[n], max_displacement);
                geometry_ptr->crossSurface(buffer[n], max_displacement);
                leakage = buffer[n].has_leaked();
                if (t + max_t_time > t_bin)
                    end_time = true;
            }
            else
            {
                scores.scoreTrace(buffer[n], t_dist); //scoring trace estimator
                geometry_ptr->noCrossing(buffer[n], t_dist, t_time);
                end_time = true;
            }
        }
        else
        {
            if (t_dist > max_displacement)
            {
                scores.scoreTrace(buffer[n], max_displacement);
                geometry_ptr->crossSurface(buffer[n], max_displacement);
                leakage = buffer[n].has_leaked();
            }
            else
            {
                scores.scoreTrace(buffer[n], t_dist); //scoring trace estimator
                geometry_ptr->noCrossing(buffer[n], t_dist, t_time);
                collision = true;
            } 
        }
    }
    return collision;
}

void Simulation::applyPopulationControl() { population_control->applyPopulationControl(buffer); }

void Simulation::changeReactivity() {current_reactivity++;}

const Scores& Simulation::getScores() const {return scores;}
Simulation::~Simulation() {}


