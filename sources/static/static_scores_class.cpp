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

#include "static_scores_class.hpp"
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"
#include "material_class.hpp"
#include "omp.h"
#include <map>

using namespace std;

StaticScores::StaticScores(): rng(omp_get_num_threads() + 1215), m(1), mKstep("kstep"), 
                              mLeakage("leakage"), mMigArea("mig_area") 
{
    size_t tot_gen = static_cast<size_t>(settings::tot_gen);
    if (settings::pi_entropy) {
        mEntropy = Tally1D({tot_gen}, "pi_entropy", 1);
        mFamilyTrees = Tally1D({tot_gen}, "pi_family_trees", 1);
    }
    if (settings::pi_flux_col)
        mFluxCol = Tally1D({(size_t) settings::x_bin}, "pi_flux_col", settings::groups);
    if (settings::pi_flux_trace)
        mFluxTrace = Tally1D({(size_t) settings::x_bin}, "pi_flux_trace", settings::groups);
    if (settings::pi_pair_distance) {
        mPairDistance = Tally1D({tot_gen}, "pi_pair_distance", 1);
        mN2 = Tally1D({tot_gen}, "pi_n2", 1);
    }
    if (settings::feynman_moment_gen) {
        mFissionSitesInner = Tally1D({(size_t) settings::x_bin}, "pi_fission_sites_inner", 1);
        mFissionSitesOuter = 
            make_pair(Tally1D({(size_t) settings::x_bin}, "pi_fission_sites_av_outer", 1), 
                      Tally1D({(size_t) settings::x_bin}, "pi_fission_sites_var_outer", 1));
        mFeynmanMomentGen = Tally1D({(size_t) settings::x_bin}, "pi_feynman_moment_gen", 1);
    }
    if (settings::feynman_moment_rep) {
        mFissionSitesRep = Tally1D({(size_t) settings::x_bin}, "pi_fission_sites_rep", 1);
        mFeynmanMomentRep = Tally1D({(size_t) settings::x_bin}, "pi_feynman_moment_rep", 1);
    }

    if (settings::pi_emission_density_gen) {
        mEmissionDensityInner = Tally1D({(size_t) settings::x_bin}, "pi_emission_density_inner", settings::groups);
        mEmissionDensityOuter =
            make_pair(Tally1D({(size_t) settings::x_bin}, "pi_emission_density_av_outer", settings::groups),
                      Tally1D({(size_t) settings::x_bin}, "pi_emission_density_var_outer", settings::groups));
        mTotEmissionDensityInner = Tally1D({(size_t) settings::x_bin}, "pi_tot_emission_density_inner", 1);
        mTotEmissionDensityOuter =
            make_pair(Tally1D({(size_t) settings::x_bin}, "pi_tot_emission_density_av_outer", 1),
                      Tally1D({(size_t) settings::x_bin}, "pi_tot_emission_density_var_outer", 1));
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep = Tally1D({(size_t) settings::x_bin}, "pi_emission_density_rep",settings::groups);
        mTotEmissionDensityRep = Tally1D({(size_t) settings::x_bin}, "pi_tot_emission_density_rep", 1);
    }

    if (settings::pi_collision_density_gen) {
        mCollisionDensityInner = Tally1D({(size_t) settings::x_bin}, "pi_collision_density_inner", settings::groups);
        mCollisionDensityOuter =
            make_pair(Tally1D({(size_t) settings::x_bin}, "pi_collision_density_av_outer", settings::groups),
                      Tally1D({(size_t) settings::x_bin}, "pi_collision_density_var_outer", settings::groups));
        mTotCollisionDensityInner = Tally1D({(size_t) settings::x_bin}, "pi_tot_collision_density_inner", 1);
        mTotCollisionDensityOuter =
            make_pair(Tally1D({(size_t) settings::x_bin}, "pi_tot_collision_density_av_outer", 1),
                      Tally1D({(size_t) settings::x_bin}, "pi_tot_collision_density_var_outer", 1));
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep = Tally1D({(size_t) settings::x_bin}, "pi_collision_density_rep",settings::groups);
        mTotCollisionDensityRep = Tally1D({(size_t) settings::x_bin}, "pi_tot_collision_density_rep", 1);
    }

    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry_ptr = make_unique<GeometryRef>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry_ptr = make_unique<GeometryPer>(parseGeometry());
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry_ptr = make_unique<GeometryLeak>(parseGeometry());
    else
        fatal_error("invalid choice of boundary conditions", __FILE__, __LINE__);
}

void StaticScores::scoreGeneration(const double& _k, vector<Particle>& _buffer, 
                                   const int& _gen) {
    if (_gen > settings::passive_gen 
            and settings::simulation_mode == settings::SimulationMode::STATIC) {
        mKstep.scoreTmp({0}, 0, _k);
        if (settings::pi_emission_density_gen) {
            mEmissionDensityInner.score(_gen - settings::passive_gen);
            mEmissionDensityInner.reinitializeTmp();
            mTotEmissionDensityInner.score(_gen - settings::passive_gen);
            mTotEmissionDensityInner.reinitializeTmp();
        }
        if (settings::pi_collision_density_gen) {
            mCollisionDensityInner.score(_gen - settings::passive_gen);
            mCollisionDensityInner.reinitializeTmp();
            mTotCollisionDensityInner.score(_gen - settings::passive_gen);
            mTotCollisionDensityInner.reinitializeTmp();
        }

        if (settings::feynman_moment_gen) {
            scoreFissionSitesInnerTmp(_buffer);
            mFissionSitesInner.score(_gen - settings::passive_gen);
        }
        if (_gen == settings::tot_gen and settings::feynman_moment_rep)
            scoreFissionSitesRepTmp(_buffer);
    }
    if (_gen < settings::tot_gen 
             and settings::simulation_mode == settings::SimulationMode::STATIC) {
        if (settings::pi_entropy) {
            scoreEntropy(_buffer, _gen);
            scoreFamilyTrees(_buffer, _gen);
        }
        if (settings::pi_pair_distance)
            scorePairDistance(_buffer, _gen);
    }
//    if (_gen > settings::passive_gen)
//        saveWeights(_buffer, _gen);
}

void StaticScores::saveWeights(const vector<Particle>& _buffer, const size_t _gen) {
    string name = settings::output + to_string(_gen);
    ofstream file(name);

    for (const auto& obj: _buffer)
        file << obj.wgt() << " ";
}



void StaticScores::scoreFissionSitesInner(const int& _gen) {
    if (settings::feynman_moment_gen) 
        mFissionSitesInner.score(_gen - settings::passive_gen);
}

void StaticScores::scoreReplica() {
    mKstep.score(m);
    mLeakage.score(m);
    mMigArea.score(m);

    if (settings::pi_flux_col)
        mFluxCol.score(m);
    if (settings::pi_flux_trace)
        mFluxTrace.score(m);
    if (settings::pi_entropy) {
        mFamilyTrees.score(m);
        mEntropy.score(m);
    }

    if (settings::pi_pair_distance) {
        mN2.score(m);
        mPairDistance.score(m);
    }

    if (settings::feynman_moment_gen) {
        mFissionSitesOuter.first.scoreTmp(mFissionSitesInner.getAverage()[0], 0);
        mFissionSitesOuter.second.scoreTmp(mFissionSitesInner.getVariance()[0], 0);
        mFissionSitesOuter.first.score(m);
        mFissionSitesOuter.second.score(m);
        scoreFeynmanMomentGenTmp();
        mFeynmanMomentGen.score(m);
    }

    if (settings::feynman_moment_rep) {
        mFissionSitesRep.score(m);
    }

    if (settings::pi_emission_density_gen) {
        for (size_t g = 0; g < settings::groups; g++) {
            mEmissionDensityOuter.first.scoreTmp(mEmissionDensityInner.getAverage()[g], g);
            mEmissionDensityOuter.second.scoreTmp(mEmissionDensityInner.getVariance()[g], g);
        }
        mEmissionDensityOuter.first.score(m);
        mEmissionDensityOuter.second.score(m);
        mTotEmissionDensityOuter.first.scoreTmp(mTotEmissionDensityInner.getAverage()[0], 0);
        mTotEmissionDensityOuter.second.scoreTmp(mTotEmissionDensityInner.getVariance()[0], 0);
        mTotEmissionDensityOuter.first.score(m);
        mTotEmissionDensityOuter.second.score(m);
    }

    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep.score(m);
        mTotEmissionDensityRep.score(m);
    }

    if (settings::pi_collision_density_gen) {
        for (size_t g = 0; g < settings::groups; g++) {
            mCollisionDensityOuter.first.scoreTmp(mCollisionDensityInner.getAverage()[g], g);
            mCollisionDensityOuter.second.scoreTmp(mCollisionDensityInner.getVariance()[g], g);
        }
        mCollisionDensityOuter.first.score(m);
        mCollisionDensityOuter.second.score(m);
        mTotCollisionDensityOuter.first.scoreTmp(mTotCollisionDensityInner.getAverage()[0], 0);
        mTotCollisionDensityOuter.second.scoreTmp(mTotCollisionDensityInner.getVariance()[0], 0);
        mTotCollisionDensityOuter.first.score(m);
        mTotCollisionDensityOuter.second.score(m);
    }

    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep.score(m);
        mTotCollisionDensityRep.score(m);
    }
}

void StaticScores::scoreCol(const Particle& obj)
{
    const Material& mat = geometry->getRegion(obj.in_region()).material();
    int g = obj.Eg(); 
    double w = obj.wgt();
    int bin = (obj.r().x() + settings::x_size)/settings::dx;
    checkIndex(bin);

    mMigArea.scoreTmp({0},0,settings::IFn*w*(mat.E_c()[g] + mat.E_f()[g])/mat.E_t()[g]*
            pow((obj.r() - obj.b_r()).norm(),2)); 

    if (settings::pi_flux_col)
        mFluxCol.scoreTmp({(size_t) bin}, g, w/mat.E_t()[g]);

    if (settings::pi_collision_density_gen) {
        mCollisionDensityInner.scoreTmp({static_cast<size_t>(bin)}, obj.Eg(), w);
        mTotCollisionDensityInner.scoreTmp({static_cast<size_t>(bin)}, 0, w);
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep.scoreTmp({static_cast<size_t>(bin)}, obj.Eg(), w);
        mTotCollisionDensityRep.scoreTmp({static_cast<size_t>(bin)}, 0, w);
    }
//    precu_density.scoreTmp({0, (size_t) bin}, 0, (w*c_p.nud[g]*c_p.sigma_f[g])/(c_p.sigma_tot[g]*p.lambda*obj.dir().v()));
}

void StaticScores::scoreEntropy(const std::vector<Particle>& n_buffer, const size_t& _gen) {
    vector<double> tmp(settings::x_bin, 0); 
    double tot_wgt = 0;
    for (const auto& obj : n_buffer) {
        int bin = (obj.r().x() + settings::x_size)/settings::dx;
        if (obj.is_alive())
            tmp[bin] += obj.wgt();
        tot_wgt += obj.wgt();
    }

    double tmp_entropy = 0;
    for (int i = 0; i < settings::x_bin; i++) {
        if (tmp[i] > 0)
            tmp_entropy -= tmp[i]/tot_wgt* log2(tmp[i]/tot_wgt);
    }

    mEntropy.scoreTmp({_gen}, 0, tmp_entropy);
}

void StaticScores::scoreFamilyTrees(const std::vector<Particle>& n_buffer, const size_t& _gen) {
    std::map<int, int> family_map;
    for (const auto& obj : n_buffer) {
        family_map.emplace(make_pair(obj.tree(), 1));
    }
    mFamilyTrees.scoreTmp({_gen}, 0, family_map.size());
//    cout << family_map.size() << endl;
}

void StaticScores::scorePairDistance(std::vector<Particle>& _buffer, const size_t& _gen)
{
    double tmp_x2 = 0, tot_wgt = 0;
    size_t part_to_sample = min((int) _buffer.size(), 1000);
    shuffle(_buffer.begin(), _buffer.end(), rng.getGenerator());

    for (size_t i = 0; i < part_to_sample; i++) {
        for (size_t j = 0; j < i; j++)
            tmp_x2 += 2*_buffer[i].wgt()*_buffer[j].wgt()*pow(_buffer[i].r().x() - _buffer[j].r().x(), 2);
        tot_wgt += _buffer[i].wgt();
    }
    
    mPairDistance.scoreTmp({_gen}, 0, tmp_x2);
    mN2.scoreTmp({_gen}, 0, tot_wgt*tot_wgt);
}

void StaticScores::scoreLeakage(const double& wgt) {mLeakage.scoreTmp({0}, 0, wgt);}

void StaticScores::nextReplica() {m++;}

void StaticScores::normalize_replicas() {
    
    double th_w = (double) settings::m_replica / (double) (m-1);
    mLeakage.normalize(settings::n_particles*th_w, 
            pow(settings::n_particles,2)*settings::active_gen*th_w);
    mKstep.normalize(th_w,settings::active_gen*th_w);
    mMigArea.normalize(settings::n_particles*th_w, 
            pow(settings::n_particles,2)*settings::active_gen*th_w);
    if (settings::pi_entropy) {
        mEntropy.normalize(th_w,th_w*settings::m_replica);
        mFamilyTrees.normalize(th_w, th_w*settings::m_replica);
    }
    if (settings::pi_pair_distance) {
        NDArray<double> tmp = mN2.getAverage()[0];
        tmp *= mN2.getAverage()[0];
        vector<NDArray<double>> tmp_v(1, tmp);
        mPairDistance.normalize(mN2.getAverage(), tmp_v);
        mPairDistance.normalize(th_w, th_w*settings::m_replica);
    }
    if (settings::pi_flux_trace)
        mFluxCol.normalize(settings::dx*settings::n_particles*th_w,
                settings::n_particles*settings::active_gen*pow(settings::dx,2)*th_w*settings::m_replica);
    if (settings::pi_flux_trace)
        mFluxTrace.normalize(settings::dx*settings::n_particles*th_w,
                settings::n_particles*settings::active_gen*pow(settings::dx,2)*th_w*settings::m_replica);
    if (settings::feynman_moment_rep) {
        mFissionSitesRep.normalize(th_w, th_w);
        mFeynmanMomentRep.normalize(th_w, th_w*settings::m_replica);
    }
    if (settings::feynman_moment_gen) {
        mFissionSitesOuter.first.normalize(th_w, th_w);
        mFissionSitesOuter.second.normalize(th_w, th_w);
        mFeynmanMomentGen.normalize(th_w, th_w*settings::m_replica);
    }
    if (settings::pi_emission_density_gen) { 
        mEmissionDensityOuter.first.normalize(th_w, th_w);
        mEmissionDensityOuter.second.normalize(th_w, th_w);
        mTotEmissionDensityOuter.first.normalize(th_w, th_w);
        mTotEmissionDensityOuter.second.normalize(th_w, th_w);
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep.normalize(th_w, th_w);
        mTotEmissionDensityRep.normalize(th_w, th_w);
    }
    if (settings::pi_collision_density_gen) {
        mCollisionDensityOuter.first.normalize(th_w, th_w);
        mCollisionDensityOuter.second.normalize(th_w, th_w);
        mTotCollisionDensityOuter.first.normalize(th_w, th_w);
        mTotCollisionDensityOuter.second.normalize(th_w, th_w);
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep.normalize(th_w, th_w);
        mTotCollisionDensityRep.normalize(th_w, th_w);
    }
}

void StaticScores::normalize_generations() {
    mKstep.normalizeTmp(settings::active_gen);
    mLeakage.normalizeTmp(settings::active_gen);
    mMigArea.normalizeTmp(settings::active_gen);
    if (settings::pi_flux_col)
        mFluxCol.normalizeTmp(settings::active_gen);
    if (settings::pi_flux_trace)
        mFluxTrace.normalizeTmp(settings::active_gen);
//    if (settings::feynman_moment)
//        mFissionSites.normalize(settings::dx, settings::dx*settings::dx);
}
    
void StaticScores::reinitialize() {
    mKstep.reinitializeTmp();
    mLeakage.reinitializeTmp();
    mMigArea.reinitializeTmp();
    if (settings::pi_entropy) {
        mEntropy.reinitializeTmp();
        mFamilyTrees.reinitializeTmp();
    }
    if (settings::pi_flux_col)
        mFluxCol.reinitializeTmp();
    if (settings::pi_flux_trace)
        mFluxTrace.reinitializeTmp();
    if (settings::pi_pair_distance) {
        mPairDistance.reinitializeTmp();
        mN2.reinitializeTmp();
    }
    if (settings::feynman_moment_gen) {
        mFissionSitesInner.reinitialize();
        mFissionSitesOuter.first.reinitializeTmp();
        mFissionSitesOuter.second.reinitializeTmp();
        mFeynmanMomentGen.reinitializeTmp();
    }
    if (settings::feynman_moment_rep)  
        mFissionSitesRep.reinitializeTmp();
    if (settings::pi_emission_density_gen) {
        mEmissionDensityInner.reinitialize();
        mEmissionDensityOuter.first.reinitializeTmp();
        mEmissionDensityOuter.second.reinitializeTmp();
        mTotEmissionDensityInner.reinitialize();
        mTotEmissionDensityOuter.first.reinitializeTmp();
        mTotEmissionDensityOuter.second.reinitializeTmp();
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep.reinitializeTmp();
        mTotEmissionDensityRep.reinitializeTmp();
    }
    if (settings::pi_collision_density_gen) {
        mCollisionDensityInner.reinitialize();
        mCollisionDensityOuter.first.reinitializeTmp();
        mCollisionDensityOuter.second.reinitializeTmp();
        mTotCollisionDensityInner.reinitialize();
        mTotCollisionDensityOuter.first.reinitializeTmp();
        mTotCollisionDensityOuter.second.reinitializeTmp();
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep.reinitializeTmp();
        mTotCollisionDensityRep.reinitializeTmp();
    }
}

void StaticScores::save() {
    if (settings::pi_entropy) {
        mEntropy.save(settings::output);
        mFamilyTrees.save(settings::output);
    }
    if (settings::pi_flux_col)
        mFluxCol.save(settings::output);
    if (settings::pi_flux_trace)
        mFluxTrace.save(settings::output);
    if (settings::pi_pair_distance)
        mPairDistance.save(settings::output);
    if (settings::feynman_moment_gen) {
        mFissionSitesOuter.first.save(settings::output);
        mFissionSitesOuter.second.save(settings::output);
        mFeynmanMomentGen.save(settings::output);
    }
    if (settings::feynman_moment_rep) {
        mFissionSitesRep.save(settings::output);
        mFeynmanMomentRep.save(settings::output);
    }
    if (settings::pi_emission_density_gen) {
        mEmissionDensityOuter.first.save(settings::output);
        mEmissionDensityOuter.second.save(settings::output);
        mTotEmissionDensityOuter.first.save(settings::output);
        mTotEmissionDensityOuter.second.save(settings::output);
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep.save(settings::output);
        mTotEmissionDensityRep.save(settings::output);
    }
    if (settings::pi_collision_density_gen) { 
        mCollisionDensityOuter.first.save(settings::output);
        mCollisionDensityOuter.second.save(settings::output);
        mTotCollisionDensityOuter.first.save(settings::output);
        mTotCollisionDensityOuter.second.save(settings::output);
    }
    if (settings::pi_collision_density_rep) {
        mCollisionDensityRep.save(settings::output);
        mTotCollisionDensityRep.save(settings::output);
    }
}


Tally1D const& StaticScores::FluxCol() const {return mFluxCol;}
Tally1D const& StaticScores::FluxTrace() const {return mFluxTrace;}
Tally1D const& StaticScores::Entropy() const {return mEntropy;}
Tally1D const& StaticScores::FamilyTrees() const {return mFamilyTrees;}
Tally1D const& StaticScores::PairDistance() const {return mPairDistance;}
Tally1D const& StaticScores::N2() const {return mN2;}
Tally1D const& StaticScores::FissionSitesInner() const {return mFissionSitesInner;}
pair<Tally1D, Tally1D> const& StaticScores::FissionSitesOuter() const {return mFissionSitesOuter;}
Tally1D const& StaticScores::FeynmanMomentGen() const {return mFeynmanMomentGen;}
Tally1D const& StaticScores::FissionSitesRep() const {return mFissionSitesRep;}
Tally1D const& StaticScores::FeynmanMomentRep() const {return mFeynmanMomentRep;}
Tally1D const& StaticScores::EmissionDensityInner() const {return mEmissionDensityInner;}
pair<Tally1D, Tally1D> const& StaticScores::EmissionDensityOuter() const {return mEmissionDensityOuter;}
Tally1D const& StaticScores::TotEmissionDensityInner() const {return mTotEmissionDensityInner;}
pair<Tally1D, Tally1D> const& StaticScores::TotEmissionDensityOuter() const {return mTotEmissionDensityOuter;}
Tally1D const& StaticScores::EmissionDensityRep() const {return mEmissionDensityRep;}
Tally1D const& StaticScores::TotEmissionDensityRep() const {return mTotEmissionDensityRep;}
Tally1D const& StaticScores::CollisionDensityInner() const {return mCollisionDensityInner;}
pair<Tally1D, Tally1D> const& StaticScores::CollisionDensityOuter() const {return mCollisionDensityOuter;}
Tally1D const& StaticScores::TotCollisionDensityInner() const {return mTotCollisionDensityInner;}
pair<Tally1D, Tally1D> const& StaticScores::TotCollisionDensityOuter() const {return mTotCollisionDensityOuter;}
Tally1D const& StaticScores::CollisionDensityRep() const {return mCollisionDensityRep;}
Tally1D const& StaticScores::TotCollisionDensityRep() const {return mTotCollisionDensityRep;}
Tally0D const& StaticScores::Kstep() const {return mKstep;}
Tally0D const& StaticScores::Leakage() const {return mLeakage;}
Tally0D const& StaticScores::MigArea() const {return mMigArea;}

void StaticScores::scoreMigAreaTmp(const Particle& obj) {
    double w = obj.wgt();
    mMigArea.scoreTmp({0},0,settings::IFn*w*pow((obj.r() - obj.b_r()).norm(),2)); 
}

void StaticScores::scoreTrace(const Particle& obj, const double& t_dist)
{
    const int& g = obj.Eg();
    const double& w = obj.wgt();
    double distance_left = t_dist;
    const double& mu = obj.dir().x();
    const double& x_i = obj.r().x();
    const double& x_f = obj.r().x() + t_dist*mu;

    int bin_i = (x_i + settings::x_size)/settings::dx; //select initial bin
    int bin_f = (x_f + settings::x_size)/settings::dx; //select final bin
    checkIndicies(bin_i, bin_f);

    if (mu > 0)
    {
        for (size_t i = bin_i; i < (size_t) min(settings::x_bin,bin_f+1) ; i++)
        {
            double diff = min(x_f, (i+1)*settings::dx - settings::x_size) - max(x_i,i*settings::dx - settings::x_size); 
            distance_left -= diff/mu;

            if (settings::pi_flux_trace)
                mFluxTrace.scoreTmp({(size_t) i}, g, settings::IFn*w*diff/mu);
        }
    }
    else if (mu < 0)
    {
        for (int i = bin_i; i >= max(0, bin_f); i--)
        {
            double diff = min(x_i, (i+1)*settings::dx - settings::x_size) - max(x_f,i*settings::dx - settings::x_size); 
            distance_left -= diff/abs(mu);
            if (settings::pi_flux_trace)
                mFluxTrace.scoreTmp({(size_t) i}, g, settings::IFn*w*diff/abs(mu));
            if (i == 0)
                break;
        }
    }
}

void StaticScores::scoreFeynmanMomentGenTmp() {
    NDArray<double> tmp = mFissionSitesInner.getVariance()[0];
    tmp /= mFissionSitesInner.getAverage()[0];
    mFeynmanMomentGen.scoreTmp(tmp, 0);
}

void StaticScores::scoreFissionSitesInnerTmp(const vector<Particle>& _sites) {
    mFissionSitesInner.reinitializeTmp();
    for (const auto& obj: _sites) { 
        int bin = (obj.r().x() + settings::x_size)/settings::dx;
        checkIndex(bin);
        mFissionSitesInner.scoreTmp({static_cast<size_t>(bin) }, 0, obj.wgt());
    }
}

void StaticScores::scoreFissionSitesRepTmp(const vector<Particle>& _sites) {
    for (const auto& obj: _sites) { 
        int bin = (obj.r().x() + settings::x_size)/settings::dx;
        checkIndex(bin);
        mFissionSitesRep.scoreTmp({static_cast<size_t>(bin) }, 0, obj.wgt());
    }
}

void StaticScores::computeFeynmanMomentRep() {
    if (settings::feynman_moment_rep) {
        NDArray<double> tmp = mFissionSitesRep.getVariance()[0];
        tmp /= mFissionSitesRep.getAverage()[0];
        mFeynmanMomentRep.scoreTmp(tmp, 0);
        mFeynmanMomentRep.score(1);
    }
}

void StaticScores::scoreEmissionDensityTmp(const Particle& obj) {
    int bin = (obj.r().x() + settings::x_size)/settings::dx;
    checkIndex(bin);
    if (settings::pi_emission_density_gen) {
        mEmissionDensityInner.scoreTmp({static_cast<size_t>(bin)}, obj.Eg(), obj.wgt());
        mTotEmissionDensityInner.scoreTmp({static_cast<size_t>(bin)}, 0, obj.wgt());
    }
    if (settings::pi_emission_density_rep) {
        mEmissionDensityRep.scoreTmp({static_cast<size_t>(bin)}, obj.Eg(), obj.wgt());
        mTotEmissionDensityRep.scoreTmp({static_cast<size_t>(bin)}, 0, obj.wgt());
    }
}
