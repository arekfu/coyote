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

#include "static_collision_kernel_class.hpp"
#include "omp.h"
#include <memory>
#include <numeric>
#include <algorithm>
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

StaticCollision::StaticCollision(): rng(8587845 + omp_get_thread_num()), mFissionBank(), mFissionSites() {
    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry_ptr = make_unique<GeometryRef>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry_ptr = make_unique<GeometryPer>(parseGeometry());
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry_ptr = make_unique<GeometryLeak>(parseGeometry());
    else
        fatal_error("invalid choice of boundary conditions", __FILE__, __LINE__);
}

void StaticCollision::applyScattering(Particle& obj) {
    int g = obj.Eg();
    const Material& mat = 
     geometry_ptr->getRegion(obj.in_region()).material();

    int new_g = integer_cumulative_law(mat.s_mat(), g, rng.R()*mat.E_s()[g]);
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));
    obj.setEg(new_g);
    obj.setV(mat.v()[new_g]);
    obj.setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
}

void StaticCollision::addFissionSite(const Particle& site, double wgt_f) {
    const Material& mat = 
        geometry_ptr->getRegion(site.in_region()).material();
    const int g = site.Eg();

    if (mat.nu()[g]*mat.E_f()[g] > 0) {
        if (settings::weighted_sites) {
            double mu = rng.Rmu();
            double phi = rng.Rphi(); 
            double sin_theta = sqrt(max(0.,1-mu*mu));
            mFissionSites.emplace_back(site);
            mFissionSites.back().multiplyWeight(wgt_f);
            mFissionSites.back().setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
        
            if(rng.R()*mat.nu()[g] < mat.nu_p()[g]) {
                int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
                mFissionSites.back().setEg(new_g);
            }
            else {
                int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
                int new_g = integer_cumulative_law(mat.xi_d(), fam, rng.R());
                mFissionSites.back().setEg(new_g);
            }
        }
        else {
            int nb_sites = floor(rng.R() + site.wgt()*wgt_f);
            for (int i = 0; i < nb_sites; i++) {
                double mu = rng.Rmu();
                double phi = rng.Rphi(); 
                double sin_theta = sqrt(max(0.,1-mu*mu));
                mFissionSites.emplace_back(site);
                mFissionSites.back().setWgt(1.);
                mFissionSites.back().setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
        
                if(rng.R()*mat.nu()[g] < mat.nu_p()[g]) {
                    int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
                    mFissionSites.back().setEg(new_g);
                }
                else {
                    int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
                    int new_g = integer_cumulative_law(mat.xi_d(), fam, rng.R());
                    mFissionSites.back().setEg(new_g);
                }
            }
        }
    }
}

void StaticCollision::clear() {
    mFissionBank.clear();
    mFissionSites.clear();
}
vector<Particle> const& StaticCollision::fissionBank() const {return mFissionBank;}
vector<Particle> const& StaticCollision::fissionSites() const {return mFissionSites;}
StaticCollision::~StaticCollision() {}

StaticAnalog::StaticAnalog(): StaticCollision() {}

void StaticAnalog::collisionNeutron(std::vector<Particle>& _buffer, const int& _n) {
    double rnd = rng.R();
    int g = _buffer[_n].Eg();
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();
    double wgt_f = mat.nu()[g]*mat.E_f()[g]/(mat.E_c()[g] + mat.E_f()[g]);


    if (rnd*mat.E_t()[g] < mat.E_c()[g]) {
        addFissionSite(_buffer[_n], wgt_f);
        _buffer[_n].kill();
    }
    else if (rnd*mat.E_t()[g] < (mat.E_f()[g] + mat.E_c()[g])) {
        addFissionSite(_buffer[_n], wgt_f);
        resolveFission(_buffer, _n);
    }
    else
        applyScattering(_buffer[_n]);
}


void StaticAnalog::resolveFission(vector<Particle>& _buffer, const int& _n)
{
    int g = _buffer[_n].Eg();
    const Material & mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();
    

    int yield = floor(rng.R() + mat.nu_p()[g]);
    for (int i = 0; i < yield; i++)
    {
        double mu = rng.Rmu();
        double phi = rng.Rphi(); 
        double sin_theta = sqrt(max(0.,1-mu*mu));
        int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              new_g, mat.v()[new_g], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }

    int delayed_yield = floor(rng.R() + mat.nu()[g] - mat.nu_p()[g]);
    for (int k = 0; k < delayed_yield; k++)
    {
        double mu = rng.Rmu();
        double phi = rng.Rphi(); 
        double sin_theta = sqrt(max(0.,1-mu*mu));
        int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
        int group = integer_cumulative_law(mat.xi_d(), fam, rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              group, mat.v()[group], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }

    _buffer[_n].kill();
}

StaticAnalog::~StaticAnalog() {}

StaticBranchless::StaticBranchless() : StaticCollision() {}

void StaticBranchless::collisionNeutron(vector<Particle>& _buffer, const int& _n)
{
    int g = _buffer[_n].Eg();
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();

    _buffer[_n].multiplyWeight((mat.E_s()[g] + mat.nu()[g]*mat.E_f()[g])/
                                mat.E_t()[g]);
    if (rng.R()*(mat.E_s()[g] + mat.nu()[g]*mat.E_f()[g]) < mat.E_s()[g])
        applyScattering(_buffer[_n]);
    else {
        addFissionSite(_buffer[_n], 1.);
        resolveFission(_buffer, _n);
    }

}

void StaticBranchless::resolveFission(vector<Particle>& _buffer, const int& _n)
{
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();
    size_t g = _buffer[_n].Eg();

    double mu = rng.Rmu();
    double phi = rng.Rphi(); 
    double sin_theta = sqrt(max(0.,1-mu*mu));

    if (rng.R()*mat.nu()[g] < mat.nu_p()[g])
    {
        int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              new_g, mat.v()[new_g], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }
    else
    {
        int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
        int group = integer_cumulative_law(mat.xi_d(), fam, rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              group, mat.v()[group], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }
    _buffer[_n].kill();
}

StaticBranchless::~StaticBranchless() {}

StaticBranching::StaticBranching() : StaticCollision() {}

void StaticBranching::collisionNeutron(vector<Particle>& _buffer, const int& _n)
{
    int g = _buffer[_n].Eg();

    const Material & mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();
    double wgt_f = mat.nu()[g]*mat.E_f()[g]/mat.E_t()[g];
    
    addFissionSite(_buffer[_n], wgt_f);
    resolveFission(_buffer,_n);
    applyScattering(_buffer[_n]);
    _buffer[_n].multiplyWeight(mat.E_s()[g]/mat.E_t()[g]);
}

void StaticBranching::resolveFission(vector<Particle>& _buffer, const int& _n)
{
    int g = _buffer[_n].Eg();
    const Material & mat = 
        geometry_ptr->getRegion(_buffer[_n].in_region()).material();
    int new_g;
    
    int yield = floor(rng.R() + mat.nu_p()[g]*mat.E_f()[g]/mat.E_t()[g]);
    for (int i = 0; i < yield; i++)
    {
        double mu = rng.Rmu();
        double phi = rng.Rphi(); 
        double sin_theta = sqrt(max(0.,1-mu*mu));
        new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              new_g, mat.v()[new_g], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }

    int delayed_yield = floor(rng.R()
            + (mat.nu()[g] - mat.nu_p()[g])*mat.E_f()[g]/mat.E_t()[g]);
    for (int k = 0; k < delayed_yield; k++)
    {
        double mu = rng.Rmu();
        double phi = rng.Rphi(); 
        double sin_theta = sqrt(max(0.,1-mu*mu));
        int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
        new_g = integer_cumulative_law(mat.xi_d(), fam, rng.R());
        mFissionBank.emplace_back(_buffer[_n].r(), 
                              Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                              new_g, mat.v()[new_g], _buffer[_n].wgt());
        mFissionBank.back().setIn_region(_buffer[_n].in_region());
        mFissionBank.back().setTree(_buffer[_n].tree());
    }
}

StaticBranching::~StaticBranching() {}
