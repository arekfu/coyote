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

#include "static_simulation_class.hpp"
#include "omp.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

StaticSimulation::StaticSimulation() :
    rng(settings::pi_seed+omp_get_thread_num()), mScores(), mOldWgt(settings::n_particles) 
{
    if (settings::simulation_mode == settings::SimulationMode::STATIC)
        mSourceSimulation = false;
    else
        mSourceSimulation = true;

    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry_ptr = make_unique<GeometryRef>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry_ptr = make_unique<GeometryPer>(parseGeometry());
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry_ptr = make_unique<GeometryLeak>(parseGeometry());
    else
        fatal_error("invalid choice of boundary conditions in static sim", __FILE__, __LINE__);

    if (settings::static_collision == settings::CollisionMode::ANALOG)
        collision_kernel = make_unique<StaticAnalog>();
    else if (settings::static_collision == settings::CollisionMode::BRANCHING)
        collision_kernel = make_unique<StaticBranching>();
    else if (settings::static_collision == settings::CollisionMode::BRANCHLESS)
        collision_kernel = make_unique<StaticBranchless>();
    else
        fatal_error("Invalid choice of collision mode in static sim", __FILE__, __LINE__);

    if (settings::static_control_mode == settings::ControlMode::ANALOG)
        population_control = make_unique<StaticNoControl>();
    else if (settings::static_control_mode == settings::ControlMode::WGT_COMBING)
        population_control = make_unique<StaticCombing>();
    else if (settings::static_control_mode == settings::ControlMode::WOR)
        population_control = make_unique<StaticWOR>();
    else
        fatal_error("Invalid choice of pop control mode in static sim", __FILE__, __LINE__);
    mCurrentBuffer.reserve(settings::n_particles);
    mFissionSites.reserve(settings::n_particles);
} 

void StaticSimulation::takePhoto(const int& rep, const string& path) const
{
    string name = path + "t" + to_string(rep) + ".txt";
    ofstream photo(name);
    if (photo)
    {
        for (const auto& obj: mCurrentBuffer)
            photo << obj.r().x() << " " << obj.r().y() << " " << obj.r().z() << " " << obj.wgt() << endl;
        photo.close();
    }
}


void StaticSimulation::source()
{
    uniform_real_distribution<double> Rl(-settings::x_size, settings::x_size);
    mCurrentBuffer.clear();
    for (int i = 0; i < settings::n_particles; i++)
    {
        double mu = rng.Rmu();
        double phi = rng.Rphi();
        double sin_theta = sqrt(max(0.,1.-mu*mu));
        mCurrentBuffer.emplace_back(Vector(rng.Rl(),0,0), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                              0, settings::v[0], 1.);
//        mCurrentBuffer.emplace_back(Vector(0,0,0), 
//                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
//                              0, settings::v[0], 1.);
        mCurrentBuffer.back().setIn_region(geometry_ptr->getRegionIndex(mCurrentBuffer.back()));
        mCurrentBuffer.back().setTree(i);
    }
    if (settings::pi_entropy and not mSourceSimulation) {
        mScores.scoreEntropy(mCurrentBuffer, 0);
        mScores.scoreFamilyTrees(mCurrentBuffer, 0);
    }
}

void StaticSimulation::antani(const int& k) {
    collision_kernel->clear();
    for (size_t n = 0; n < mCurrentBuffer.size(); n++)
    {
        bool collision = true;
        while (mCurrentBuffer[n].is_alive()) {
            if (k > settings::passive_gen and not mSourceSimulation and collision)
                mScores.scoreEmissionDensityTmp(mCurrentBuffer[n]);
            collision = transportNeutron(n, k);
            if (collision) {
                if (k > settings::passive_gen and not mSourceSimulation)
                    mScores.scoreCol(mCurrentBuffer[n]);
                collision_kernel->collisionNeutron(mCurrentBuffer, n);
            }
            else if (mCurrentBuffer[n].has_leaked() 
                     and not mSourceSimulation
                     and k > settings::passive_gen) {
                mScores.scoreLeakage(mCurrentBuffer[n].wgt());
                mScores.scoreMigAreaTmp(mCurrentBuffer[n]);
            }

            mCurrentBuffer[n].roulette(rng.R());
            if (mCurrentBuffer[n].is_alive()) {
                splitting(mCurrentBuffer, n);
            }
        }
    }
}

void StaticSimulation::update(const int& k) {
    double sum_weight = 0, w_factor = 0;
    if (settings::static_control_mode == settings::ControlMode::WOR)
        mCurrentBuffer = collision_kernel->fissionSites();
    else
        mCurrentBuffer = collision_kernel->fissionBank();

    for (const auto& obj: mCurrentBuffer)
        sum_weight += obj.wgt();
//    double k_eff = sum_weight/mOldWgt;
//    mOldWgt = sum_weight;

    double k_eff = sum_weight/settings::n_particles;
    if (sum_weight == 0)
        fatal_error("no particle left in generation " + to_string(k), 
                    __FILE__, __LINE__);

    population_control->applyPopulationControl(mCurrentBuffer);    

    sum_weight = 0;
    for (const auto& obj: mCurrentBuffer)
        sum_weight += obj.wgt();

    w_factor = sum_weight/settings::n_particles;
    for (auto& obj: mCurrentBuffer)
        obj.multiplyWeight(1./w_factor);

    mScores.scoreGeneration(k_eff, mCurrentBuffer, k);
}

bool StaticSimulation::transportNeutron(const int& n, const int& _gen)
{
    int g = mCurrentBuffer[n].Eg();
    double max_displacement = geometry_ptr->distToFirstSurface(mCurrentBuffer[n]);
    const Region& current = geometry_ptr->getRegion(mCurrentBuffer[n].in_region());
    double t_dist = -log(1-rng.R())/current.material().E_t()[g];

    if (t_dist > max_displacement) {
        if (_gen > settings::passive_gen)
            mScores.scoreTrace(mCurrentBuffer[n], max_displacement);
        geometry_ptr->crossSurface(mCurrentBuffer[n], max_displacement);
    }
    else {
        if (_gen > settings::passive_gen)
            mScores.scoreTrace(mCurrentBuffer[n], t_dist);
        geometry_ptr->noCrossing(mCurrentBuffer[n], t_dist, 0);
        return true;
    }
    return false;
}

void StaticSimulation::reinitialize() {
    mCurrentBuffer.clear();
    mTarapiaTapiocoBuffer.clear();
    mCombedBuffer.clear();
    mScores.reinitialize(); 
}

void StaticSimulation::normalize() {mScores.normalize_replicas();}
void StaticSimulation::makeScores() {
    mScores.normalize_generations();
    mScores.scoreReplica();
    mScores.nextReplica();
}

vector<Particle> const & StaticSimulation::CurrentBuffer() const {
    return mCurrentBuffer;}
vector<Particle> const & StaticSimulation::TarapiaTapiocoBuffer() const {
    return mTarapiaTapiocoBuffer;}
StaticScores const & StaticSimulation::Scores() const {return mScores;}

void StaticSimulation::computeTarapiaTapiocoBuffer(size_t cycles) {
    for (size_t k = 1; k <= cycles; k++) {
        antani(k);
        update(k);
    }
    antaniLast();
    update(cycles+1);
}

void StaticSimulation::updateTarapiaTapiocoBuffer(int n, double prob)
{
    const Particle& tmp = mCurrentBuffer[n];
    int g = tmp.Eg();
    const Material& mat = 
        geometry_ptr->getRegion(tmp.in_region()).material();

    mTarapiaTapiocoBuffer.emplace_back(tmp.r(), tmp.dir(), g, mat.v()[g], 0, 
                                tmp.wgt()/(mat.v()[g]*mat.E_t()[g]*prob));
    mTarapiaTapiocoBuffer.back().setIn_region(tmp.in_region());
    
    if (settings::collapsed_precu and mat.lam_bar() > 0) {
        double w_factor = 
            mat.beta_tot()*mat.nu()[g]*mat.E_f()[g]/(mat.E_t()[g]*mat.lam_bar()*prob);
        mTarapiaTapiocoBuffer.emplace_back(tmp.r(), 0, 0, tmp.wgt()*w_factor, 0, true);
        mTarapiaTapiocoBuffer.back().setIn_region(tmp.in_region());
    }
    else
    {
    for (size_t f = 0; f < settings::families; f++)
        if (mat.beta()[f] > 0 and mat.lambda()[f] > 0)
        {
            double w_factor = (mat.beta()[f]*mat.nu()[g]*mat.E_f()[g])/
                              (mat.lambda()[f]*prob*mat.E_t()[g]);
            double d_time = -log(1-rng.R())/mat.lambda()[f];
            mTarapiaTapiocoBuffer.emplace_back(tmp.r(), 0, d_time, tmp.wgt()*w_factor, f, true);
            mTarapiaTapiocoBuffer.back().setIn_region(tmp.in_region());

        }
    }
}

void StaticSimulation::antaniLast()
{
    mTarapiaTapiocoBuffer.clear();
    collision_kernel->clear();
    for (size_t n = 0; n < mCurrentBuffer.size(); n++)
        while (mCurrentBuffer[n].is_alive())
        {
            bool collision = transportNeutron(n, -1);
            if (collision)
            {
                if (rng.R() < settings::p)
                    updateTarapiaTapiocoBuffer(n, settings::p);
                collision_kernel->collisionNeutron(mCurrentBuffer, n);
            }
            mCurrentBuffer[n].roulette(rng.R());
            if (mCurrentBuffer[n].is_alive())
                splitting(mCurrentBuffer, n);
        }
}

void StaticSimulation::computeFeynmanMomentRep() {mScores.computeFeynmanMomentRep();}

StaticSimulation::~StaticSimulation() {}
