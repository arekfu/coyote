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
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <limits>
#include <iomanip>
#include <filesystem>
#include "/usr/lib/gcc/x86_64-linux-gnu/9/include/omp.h"
#include "scores_class.hpp"
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

Scores::Scores() 
    : m(1),
    particle_total({(size_t) settings::t_bin + 1, 1}, "particle_total", 1),
    neutron_density({(size_t) settings::t_bin + 1, 
                     (size_t) settings::x_bin}, 
                     "ct", settings::groups),
    precu_density({(size_t) settings::t_bin + 1, 
                   (size_t) settings::x_bin}, 
                   "ct_d", settings::families),
    rng(2654654 + omp_get_thread_num())
{
    size_t Nt = (size_t) settings::t_bin;
    size_t Z = (size_t) settings::x_bin;
    
    if (settings::col_correlations)
        col_correlations = Tally4D({Nt,Z,Nt,Z}, "col_cor", settings::groups);
    if (settings::trace_correlations)
        trace_correlations = Tally4D({Nt,Z,Nt,Z}, "trace_cor", settings::groups);
    if (settings::photo_correlations)
        neutron_correlations = Tally4D({Nt+1, Z, Nt+1, Z}, "gt", settings::groups);
    if (settings::col_flux)
        flux_col = Tally2D({Nt,Z},"flux_col", settings::groups);
    if (settings::trace_flux)
        flux_trace = Tally2D({Nt,Z},"flux_trace", settings::groups);
    if (settings::col_res_time)
        res_time_col = Tally2D({1,Z},"res_time_col", settings::groups);
    if (settings::trace_res_time)
        res_time_trace = Tally2D({1,Z},"res_time_trace", settings::groups);
    if (settings::entropy) {
        P_entropy = Tally2D({Nt+1,1}, "p_entropy", 1);
        N_entropy = Tally2D({Nt+1,1}, "n_entropy", 1);
        entropy = Tally1D({Nt+1}, "entropy", 1);
    }
    if (settings::source_entropy)
        source_entropy = NDArray<double>({(size_t) settings::active_gen + settings::passive_gen, 1});
    if (settings::pair_distance) {
        N2 = Tally1D({Nt+1}, "N2", 1);
        P2 = Tally1D({Nt+1}, "P2", 1);
        P_pair_dist = Tally1D({Nt+1}, "p_pair_distance", 1);
        N_pair_dist = Tally1D({Nt+1}, "n_pair_distance", 1);
    }
    if (settings::family_trees)
        family_trees = Tally1D({Nt+1}, "family_trees", 1);
    if (settings::emission_density) {
        mTotEmissionDensity = Tally2D({Nt,Z}, "tot_emission_density", 1);
        mEmissionDensity = Tally2D({Nt,Z}, "emission_density", settings::groups);
    }
    if (settings::collision_density) {
        mTotCollisionDensity = Tally2D({Nt,Z}, "tot_collision_density", 1);
        mCollisionDensity = Tally2D({Nt,Z}, "collision_density", settings::groups);
    }
    if (settings::kin_source) {
        mKinSourceNeutrons = Tally1D({Z}, "kin_source_neutrons", settings::groups);
        mCombedKinSourceNeutrons = Tally1D({Z}, "combed_kin_source_neutrons", settings::groups);
        mKinSourcePrecu= Tally1D({Z}, "kin_source_precu", settings::families);
        mCombedKinSourcePrecu = Tally1D({Z}, "combed_kin_source_precu", settings::families);
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

void Scores::scoreKinSource(const vector<Particle>& _buffer) {
    for (const auto& obj: _buffer) {
        size_t bin = (obj.r().x() + settings::x_size)/settings::dx; //select initial bin
        bin = (bin <= 0 ? 0 : bin);
        if (obj.is_precu())
            mKinSourcePrecu.scoreTmp({bin}, obj.family(), obj.wgt());
        else
            mKinSourceNeutrons.scoreTmp({bin}, obj.Eg(), obj.wgt());
    }
}
                                                                 
void Scores::scoreCombedKinSource(const vector<Particle>& _buffer) {
    for (const auto& obj: _buffer) {
        size_t bin = (obj.r().x() + settings::x_size)/settings::dx; //select initial bin
        bin = (bin <= 0 ? 0 : bin);
        if (obj.is_precu())
            mCombedKinSourcePrecu.scoreTmp({bin}, obj.family(), obj.wgt());
        else
            mCombedKinSourceNeutrons.scoreTmp({bin}, obj.Eg(), obj.wgt());
    }
}

void Scores::scoreEntropy(size_t time, const vector<Particle>& buffer)
{
    vector<double> tmp_n(settings::x_bin,0);
    vector<double> tmp_p(settings::x_bin,0);
    double neutron_wgt = 0;
    double precu_wgt = 0;

    for (const auto& obj: buffer) {
        int bin = (obj.r().x() + settings::x_size)/settings::dx; //select initial bin
        double w = obj.wgt();
        checkIndex(bin);
        if (obj.is_alive()) {
            if (obj.is_precu()){
                precu_wgt += settings::IFp*w;
                tmp_p[bin] += settings::IFp*w;
            }
            else {
                neutron_wgt += settings::IFn*w;
                tmp_n[bin] += settings::IFn*w;
            }
        }
    }

    double tmp_entropy_n = 0;
    for (int i = 0; i < settings::x_bin; i++) {
        if (tmp_n[i] > 0)
            tmp_entropy_n -= tmp_n[i]/neutron_wgt * log2(tmp_n[i]/neutron_wgt);
    }

    double tmp_entropy_p = 0;
    for (int i = 0; i < settings::x_bin; i++) {
        if (tmp_p[i] > 0)
            tmp_entropy_p -= tmp_p[i]/precu_wgt * log2(tmp_p[i]/precu_wgt);
    }

    double tmp_entropy = 0;
    double tot_wgt = neutron_wgt + precu_wgt;

    for (int i = 0; i < settings::x_bin; i++) {
        double tmp = tmp_p[i] + tmp_n[i];
        if (tmp > 0)
            tmp_entropy -= tmp/tot_wgt * log2(tmp/tot_wgt);
    }
    P_entropy.scoreTmp({time,0}, 0, tmp_entropy_p);
    N_entropy.scoreTmp({time,0}, 0, tmp_entropy_n);
    entropy.scoreTmp({time}, 0, tmp_entropy);
}

void Scores::scorePairDistance(size_t time, vector<Particle>& buffer)
{
    double n_tmp = 0, p_tmp = 0;
    double n_tot = 0, p_tot = 0;

    //sample a given number of pairs to keep CPU time manageable

    size_t part_to_sample = min((int) buffer.size(), 100);
    shuffle(buffer.begin(), buffer.end(), rng.getGenerator());

    for (size_t i = 0; i < part_to_sample; i++) {
        for (size_t j = 0; j < i;j++) {
            if (buffer[i].is_precu() and buffer[j].is_precu())
                p_tmp += 2*buffer[i].wgt()*buffer[j].wgt()*pow(buffer[i].r().x() - buffer[j].r().x(),2);
            else if (not buffer[i].is_precu() and not buffer[j].is_precu())
                n_tmp += 2*buffer[i].wgt()*buffer[j].wgt()*pow(buffer[i].r().x() - buffer[j].r().x(),2);
        }
    }
    
    //compute the weight of sampled particles
    for (size_t i = 0; i < part_to_sample; i++) {
        if (buffer[i].is_precu())
            p_tot += buffer[i].wgt();
        else
            n_tot += buffer[i].wgt();
    }

    P_pair_dist.scoreTmp({time}, 0, p_tmp);
    N_pair_dist.scoreTmp({time}, 0, n_tmp);
    N2.scoreTmp({time}, 0, n_tot*n_tot);
    P2.scoreTmp({time}, 0, p_tot*p_tot);
}

void Scores::scoreTime(size_t time, vector<Particle>& buffer) //photos
{
    double t = time*settings::dt;
    family_trees_map.clear(); 
    for (auto& obj: buffer) {
        const Material& mat = 
            geometry_ptr->getRegion(obj.in_region()).material();
        int g = obj.Eg();
        double w = obj.wgt();
        int bin = (obj.r().x() + settings::x_size)/settings::dx; //select initial bin

        checkIndex(bin);
        if (obj.is_precu() and settings::collapsed_precu) {
            double sum = 0;
            double f_wgt = 0;
            if (obj.from_eq()) {
                for (size_t i = 0; i < settings::families; i++) {
                    f_wgt = obj.wgt()*mat.lam_bar()/(mat.beta_tot()*mat.lambda()[i])*
                        mat.beta()[i]*exp(-mat.lambda()[i]*(t - obj.b_t()));
                    sum += f_wgt;
                    precu_density.scoreTmp({time, (size_t) bin}, i, settings::IFp*f_wgt);
                }
                particle_total.scoreTmp({time,0}, 0, settings::IFp*sum);
                obj.setExpWgt(sum);
            }
            else {
                for (size_t i = 0; i < settings::families; i++) {
                    f_wgt = obj.wgt()*mat.beta()[i]/mat.beta_tot()*exp(-mat.lambda()[i]*(t-obj.b_t()));
                    sum += f_wgt;
                    precu_density.scoreTmp({time, (size_t) bin}, i, settings::IFp*f_wgt);
                }
                particle_total.scoreTmp({time,0}, 0, settings::IFp*sum);
                obj.setExpWgt(sum);
            }
        }
        else if (obj.is_precu()) {
            precu_density.scoreTmp({time, (size_t) bin}, obj.family(), settings::IFp*w);
            particle_total.scoreTmp({time,0}, 0, settings::IFp*w);
        }
        else
        {
            if (settings::photo_correlations)
                neutron_correlations.scoreTmp({time, (size_t) bin}, g, settings::IFn*w);
            neutron_density.scoreTmp({time, (size_t) bin}, g, settings::IFn*w);
            particle_total.scoreTmp({time,0}, 0, settings::IFn*w);
        }
        family_trees_map.emplace(make_pair(obj.tree(), 1));
    }
    if (settings::family_trees)
        family_trees.scoreTmp({time}, 0, family_trees_map.size());
}

void Scores::scoreCol(const Particle& obj, size_t current_reactivity)
{
    const Material& mat = 
        geometry_ptr->getRegion(obj.in_region()).material();
    double w = obj.wgt();
    size_t g = obj.Eg();
    size_t bin = (obj.r().x() + settings::x_size)/settings::dx;

    if (bin >= (size_t) settings::x_bin)
        bin = settings::x_bin - 1;

    size_t t_bin = floor(obj.t()/settings::dt);
    if (t_bin >= (size_t) settings::t_bin)
        t_bin = settings::t_bin-1;
    
    if (settings::col_correlations)
        col_correlations.scoreTmp({t_bin, bin}, g, settings::IFn*w/mat.E_t(current_reactivity)[g]);
    if (settings::col_flux)
        flux_col.scoreTmp({t_bin, bin}, g, settings::IFn*w/mat.E_t(current_reactivity)[g]);
    if (settings::col_res_time)
        res_time_col.scoreTmp({0,bin}, g, settings::IFn*w/mat.E_t(current_reactivity)[g]/obj.v());
    if (settings::collision_density) {
        mCollisionDensity.scoreTmp({t_bin,bin}, g, settings::IFn*w);
        mTotCollisionDensity.scoreTmp({t_bin,bin}, 0, settings::IFn*w);
    }
}


void Scores::computeBasicScores()
{
    particle_total.score(m);
    neutron_density.score(m);
    precu_density.score(m);
    
    if (settings::col_correlations)
        col_correlations.score(m);
    if (settings::trace_correlations)
        trace_correlations.score(m);
    if (settings::photo_correlations)
        neutron_correlations.score(m);
    if (settings::trace_flux)
        flux_trace.score(m);
    if (settings::col_flux)
        flux_col.score(m);
    if (settings::col_res_time)
        res_time_trace.score(m);
    if (settings::trace_res_time)
        res_time_col.score(m);
    if (settings::entropy) {
        P_entropy.score(m);
        N_entropy.score(m);
        entropy.score(m);
    }
    if (settings::pair_distance){
        N2.score(m);
        P2.score(m);
        P_pair_dist.score(m);
        N_pair_dist.score(m);
    }
    if (settings::family_trees)
        family_trees.score(m);
    if (settings::emission_density) {
        mTotEmissionDensity.score(m);
        mEmissionDensity.score(m);
    }
    if (settings::collision_density) {
        mTotCollisionDensity.score(m);
        mCollisionDensity.score(m);
    }
    if (settings::kin_source) {
        mKinSourceNeutrons.score(m);
        mCombedKinSourceNeutrons.score(m);
        mKinSourcePrecu.score(m);
        mCombedKinSourcePrecu.score(m);
    }
}

void Scores::normalizeBasicScores()
{
    double dx = settings::dx;
    double dt = settings::dt;
    double th_w = (double) settings::m_replica / (double) (m-1);
    neutron_density.normalize(th_w, th_w);
    precu_density.normalize(th_w, th_w);
//    neutron_density.normalize(dx*th_w, pow(dx,2)*settings::m_replica*th_w);
//    precu_density.normalize(dx*th_w, pow(dx,2)*settings::m_replica*th_w);
    particle_total.normalize(th_w, settings::m_replica*th_w);

    if (settings::col_correlations)
        col_correlations.normalize(pow(dx*dt*settings::n_particles,2)*th_w, pow(dx*dt*settings::n_particles,4)*settings::m_replica*th_w);
    if (settings::trace_correlations)
        trace_correlations.normalize(pow(dx*dt*settings::n_particles,2)*th_w, pow(dx*dt*settings::n_particles,4)*settings::m_replica*th_w);
    if (settings::photo_correlations)
        neutron_correlations.normalize(th_w*pow(dx,2), pow(dx,4)*settings::m_replica*th_w);
    if (settings::col_flux)
        flux_col.normalize(th_w*settings::n_particles*dt*dx, settings::m_replica*th_w*pow(dt*dx*settings::n_particles,2));
    if (settings::trace_flux)
        flux_trace.normalize(th_w*settings::n_particles*dt*dx, settings::m_replica*th_w*pow(dt*dx*settings::n_particles,2));
    if (settings::col_res_time)
        res_time_col.normalize(th_w*settings::n_particles, settings::n_particles*settings::m_replica*th_w);
    if (settings::trace_res_time)
        res_time_trace.normalize(th_w*settings::n_particles, settings::n_particles*settings::m_replica*th_w);
    if (settings::entropy) {
        P_entropy.normalize(th_w, settings::m_replica*th_w);
        N_entropy.normalize(th_w, settings::m_replica*th_w);
        entropy.normalize(th_w, settings::m_replica*th_w);
    }
    if (settings::pair_distance) {
        P_pair_dist.normalize(P2.getAverage(), P2.getVariance());
        N_pair_dist.normalize(N2.getAverage(), N2.getVariance());
        P_pair_dist.normalize(th_w, th_w*settings::m_replica);
        N_pair_dist.normalize(th_w, th_w*settings::m_replica);
    }
    if (settings::family_trees)
        family_trees.normalize(th_w, th_w*settings::m_replica);

    if (settings::emission_density) {
        mTotEmissionDensity.normalize(th_w*dt*dx, th_w*pow(dt*dx,2));
        mEmissionDensity.normalize(th_w*dt*dx, th_w*pow(dt*dx,2));
    }
    if (settings::collision_density) {
        mTotCollisionDensity.normalize(th_w, th_w);
        mCollisionDensity.normalize(th_w, th_w);
    }
    if (settings::kin_source) {
        mKinSourceNeutrons.normalize(th_w, th_w);
        mCombedKinSourceNeutrons.normalize(th_w, th_w);
        mKinSourcePrecu.normalize(th_w, th_w);
        mCombedKinSourcePrecu.normalize(th_w, th_w);
    }
}

Tally2D const & Scores::getParticleTotal() const {return particle_total;}
Tally2D const & Scores::getNeutronDensity() const {
    return neutron_density;
}
Tally2D const & Scores::getPrecuDensity() const {
    return precu_density;
}
Tally2D const & Scores::getResTimeCol() const {return res_time_col;}
Tally2D const & Scores::getResTimeTrace() const {return res_time_trace;}
Tally2D const & Scores::getFluxCol() const {return flux_col;}
Tally2D const & Scores::getFluxTrace() const {return flux_trace;}
Tally4D const & Scores::getNeutronCorrelations() const {return neutron_correlations;}
Tally4D const & Scores::getColCorrelations() const {return col_correlations;}
Tally4D const & Scores::getTraceCorrelations() const {return trace_correlations;}
Tally2D const & Scores::getPEntropy() const {return P_entropy;}
Tally2D const & Scores::getNEntropy() const {return N_entropy;}
Tally1D const & Scores::getEntropy() const {return entropy;}
NDArray<double> const & Scores::getSourceEntropy() const {return source_entropy;}
Tally1D const & Scores::getPPairDist() const {return P_pair_dist;}
Tally1D const & Scores::getNPairDist() const {return N_pair_dist;}
Tally1D const & Scores::getN2() const {return N2;}
Tally1D const & Scores::getP2() const {return P2;}
Tally1D const & Scores::getFamilyTrees() const {return family_trees;}
Tally2D const & Scores::getEmissionDensity() const {return mEmissionDensity;}
Tally2D const & Scores::getCollisionDensity() const {return mCollisionDensity;}
Tally2D const & Scores::getTotEmissionDensity() const {return mTotEmissionDensity;}
Tally2D const & Scores::getTotCollisionDensity() const {return mTotCollisionDensity;}
Tally1D const & Scores::getKinSourceNeutrons() const {return mKinSourceNeutrons;}
Tally1D const & Scores::getCombedKinSourceNeutrons() const {return mCombedKinSourceNeutrons;}
Tally1D const & Scores::getKinSourcePrecu() const {return mKinSourcePrecu;}
Tally1D const & Scores::getCombedKinSourcePrecu() const {return mCombedKinSourcePrecu;}

void Scores::setSourceEntropy(const NDArray<double> & obj) {source_entropy = obj;}

Scores& Scores::operator+=(const Scores& rhs)
{
    particle_total += rhs.getParticleTotal();
    neutron_density += rhs.getNeutronDensity();
    precu_density += rhs.getPrecuDensity();

    if (settings::col_correlations)
        col_correlations += rhs.getColCorrelations();
    if (settings::trace_correlations)
        trace_correlations += rhs.getTraceCorrelations();
    if (settings::photo_correlations)
        neutron_correlations += rhs.getNeutronCorrelations();
    if (settings::col_flux)
        flux_col += rhs.getFluxCol();
    if (settings::trace_flux)
        flux_trace += rhs.getFluxTrace();
    if (settings::col_res_time)
        res_time_col += rhs.getResTimeCol();
    if (settings::trace_res_time)
        res_time_trace += rhs.getResTimeTrace();
    if (settings::entropy) {
        P_entropy += rhs.getPEntropy();
        N_entropy += rhs.getNEntropy();
        entropy += rhs.getEntropy();
    }
    if (settings::pair_distance) {
        N2 += rhs.getN2();
        P2 += rhs.getP2();
        P_pair_dist += rhs.getPPairDist();
        N_pair_dist += rhs.getNPairDist();
    }
    if (settings::family_trees)
        family_trees += rhs.getFamilyTrees();
    if (settings::emission_density) {
        mTotEmissionDensity += rhs.getTotEmissionDensity();
        mEmissionDensity += rhs.getEmissionDensity();
    }
    if (settings::collision_density) {
        mTotCollisionDensity += rhs.getTotCollisionDensity();
        mCollisionDensity += rhs.getCollisionDensity();
    }
    if (settings::kin_source) {
        mKinSourceNeutrons += rhs.getKinSourceNeutrons();
        mCombedKinSourceNeutrons += rhs.getCombedKinSourceNeutrons();
        mKinSourcePrecu += rhs.getKinSourcePrecu();
        mCombedKinSourcePrecu += rhs.getCombedKinSourcePrecu();
    }
    return *this;
}

void Scores::addScores(const Scores& rhs)
{
    particle_total += rhs.getParticleTotal();
    neutron_density += rhs.getNeutronDensity();
    precu_density += rhs.getPrecuDensity();

    if (settings::col_correlations)
        col_correlations += rhs.getColCorrelations();
    if (settings::trace_correlations)
        trace_correlations += rhs.getTraceCorrelations();
    if (settings::photo_correlations)
        neutron_correlations += rhs.getNeutronCorrelations();
    if (settings::col_flux)
        flux_col += rhs.getFluxCol();
    if (settings::trace_flux)
        flux_trace += rhs.getFluxTrace();
    if (settings::col_res_time)
        res_time_col += rhs.getResTimeCol();
    if (settings::trace_res_time)
        res_time_trace += rhs.getResTimeTrace();
    if (settings::entropy) {
        P_entropy += rhs.getPEntropy();
        N_entropy += rhs.getNEntropy();
        entropy += rhs.getEntropy();
    }
    if (settings::pair_distance) {
        N2 += rhs.getN2();
        P2 += rhs.getP2();
        P_pair_dist += rhs.getPPairDist();
        N_pair_dist += rhs.getNPairDist();
    }
    if (settings::family_trees)
        family_trees += rhs.getFamilyTrees();
    if (settings::emission_density) {
        mTotEmissionDensity += rhs.getTotEmissionDensity();
        mEmissionDensity += rhs.getEmissionDensity();
    }
    if (settings::collision_density) {
        mTotCollisionDensity += rhs.getTotCollisionDensity();
        mCollisionDensity += rhs.getCollisionDensity();
    }
    if (settings::kin_source) {
        mKinSourceNeutrons += rhs.getKinSourceNeutrons();
        mCombedKinSourceNeutrons += rhs.getCombedKinSourceNeutrons();
        mKinSourcePrecu += rhs.getKinSourcePrecu();
        mCombedKinSourcePrecu += rhs.getCombedKinSourcePrecu();
    }
}

void Scores::saveScores(const std::string& path)
{
    particle_total.save(path);
    neutron_density.save(path);
    precu_density.save(path);
    //p.save("../data/saved/input/param.txt");
    //geometry_ptr->save("../data/saved/input/");
    
    if (settings::col_correlations)
        col_correlations.save(path);
    if (settings::trace_correlations)
        trace_correlations.save(path);
    if (settings::photo_correlations)
        neutron_correlations.save(path);
    if (settings::col_flux)
        flux_col.save(path);
    if (settings::trace_flux)
        flux_trace.save(path);
    if (settings::col_res_time)
        res_time_col.save(path);
    if (settings::trace_res_time)
        res_time_trace.save(path);
    if (settings::entropy) {
        P_entropy.save(path);
        N_entropy.save(path);
        entropy.save(path);
    }
    if (settings::pair_distance) {
        P_pair_dist.save(path);
        N_pair_dist.save(path);
    }
    if (settings::family_trees)
        family_trees.save(path);
    if (settings::emission_density) {
        mTotEmissionDensity.save(path);
        mEmissionDensity.save(path);
    }
    if (settings::collision_density) {
        mTotCollisionDensity.save(path);
        mCollisionDensity.save(path);
    }
    if (settings::kin_source) {
        mKinSourceNeutrons.save(path);
        mCombedKinSourceNeutrons.save(path);
        mKinSourcePrecu.save(path);
        mCombedKinSourcePrecu.save(path);
    }
}

void Scores::reinitializeTempScores() 
{
    particle_total.reinitializeTmp();
    neutron_density.reinitializeTmp();
    precu_density.reinitializeTmp();
    
    if (settings::col_correlations)
        col_correlations.reinitializeTmp();
    if (settings::trace_correlations)
        trace_correlations.reinitializeTmp();
    if (settings::photo_correlations)
        neutron_correlations.reinitializeTmp();
    if (settings::col_flux)
        flux_col.reinitializeTmp();
    if (settings::trace_flux)
        flux_trace.reinitializeTmp();
    if (settings::col_res_time)
        res_time_col.reinitializeTmp();
    if (settings::trace_res_time)
        res_time_trace.reinitializeTmp();
    if (settings::entropy) {
        P_entropy.reinitializeTmp();
        N_entropy.reinitializeTmp();
        entropy.reinitializeTmp();
    }
    if (settings::pair_distance) {
        N2.reinitializeTmp();
        P2.reinitializeTmp();
        P_pair_dist.reinitializeTmp();
        N_pair_dist.reinitializeTmp();
    }
    if (settings::family_trees)
        family_trees.reinitializeTmp();

    if (settings::emission_density) {
        mTotEmissionDensity.reinitializeTmp();
        mEmissionDensity.reinitializeTmp();
    }
    if (settings::collision_density) {
        mTotCollisionDensity.reinitializeTmp();
        mCollisionDensity.reinitializeTmp();
    }
    if (settings::kin_source) {
        mKinSourceNeutrons.reinitializeTmp();
        mCombedKinSourceNeutrons.reinitializeTmp();
        mKinSourcePrecu.reinitializeTmp();
        mCombedKinSourcePrecu.reinitializeTmp();
    }
}

void Scores::nextReplica() {m++;}

void Scores::scoreTrace(const Particle& obj, const double& t_dist)
{
    int g = obj.Eg();
    double w = obj.wgt();
    double distance_left = t_dist;
    double v = obj.v();
    double mu = obj.dir().x();
    double x_i = obj.r().x();
    double x_f = obj.r().x() + t_dist*mu;

    size_t t_bin = floor(obj.t()/settings::dt);
    if (t_bin >= (size_t) settings::t_bin)
        t_bin = settings::t_bin-1;

    int bin_i = (x_i + settings::x_size)/settings::dx; //select initial bin
    int bin_f = (x_f + settings::x_size)/settings::dx; //select final bin
    checkIndicies(bin_i, bin_f);

    if (mu > 0)
    {
        for (size_t i = bin_i; i < (size_t) min(settings::x_bin,bin_f+1) ; i++)
        {
            double diff = min(x_f, (i+1)*settings::dx - settings::x_size) - max(x_i,i*settings::dx - settings::x_size); 
            distance_left -= diff/mu;

            if (settings::trace_res_time)
                res_time_trace.scoreTmp({0,i}, g, settings::IFn*w*diff/v/mu);
            if (settings::trace_flux)
                flux_trace.scoreTmp({t_bin, i}, g, settings::IFn*w*diff/mu);
            if (settings::trace_correlations)
                trace_correlations.scoreTmp({t_bin, i}, g, w*diff/mu*settings::IFn);
        }
    }
    else if (mu < 0)
    {
        for (int i = bin_i; i >= max(0, bin_f); i--)
        {
            double diff = min(x_i, (i+1)*settings::dx - settings::x_size) - max(x_f,i*settings::dx - settings::x_size); 
            distance_left -= diff/abs(mu);
            if (settings::trace_res_time)
                res_time_trace.scoreTmp({0,(size_t) i}, g, settings::IFn*w*diff/v/abs(mu));
            if (settings::trace_flux)
                flux_trace.scoreTmp({t_bin, (size_t) i}, g, settings::IFn*w*diff/abs(mu));
            if (settings::trace_correlations)
                trace_correlations.scoreTmp({t_bin, (size_t) i}, g, settings::IFn*w*diff/abs(mu));
            if (i == 0)
                break;
        }
    }
    else
    {
        if (settings::trace_res_time)
            res_time_trace.scoreTmp({0, (size_t) bin_i}, g, settings::IFn*w*t_dist/v);
        if (settings::trace_correlations)
            trace_correlations.scoreTmp({t_bin, (size_t) bin_i}, g, settings::IFn*w*t_dist);
    }
}

void Scores::scoreEmissionDensityTmp(const Particle& obj, const size_t t_bin) {
    int bin = (obj.r().x() + settings::x_size)/settings::dx;
    checkIndex(bin);
    mEmissionDensity.scoreTmp({t_bin, static_cast<size_t>(bin)}, obj.Eg(), obj.wgt()*settings::IFn);
    mTotEmissionDensity.scoreTmp({t_bin, static_cast<size_t>(bin)}, 0, obj.wgt()*settings::IFn);
}

Scores::~Scores() {}
