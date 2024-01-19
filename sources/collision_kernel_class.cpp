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

#include "collision_kernel_class.hpp"
#include "omp.h"
#include <memory>
#include <numeric>
#include <algorithm>
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

CollisionKernel::CollisionKernel() : rng(888400 + omp_get_thread_num())
{
    if (settings::boundary_type == settings::BoundaryType::REFLECTIVE) 
        geometry_ptr = make_unique<GeometryRef>(parseGeometry()); 
    else if (settings::boundary_type == settings::BoundaryType::PERIODIC) 
        geometry_ptr = make_unique<GeometryPer>(parseGeometry());
    else if (settings::boundary_type == settings::BoundaryType::VACUUM) 
        geometry_ptr = make_unique<GeometryLeak>(parseGeometry());
    else
        fatal_error("invalid choice of boundary conditions", __FILE__, __LINE__);
}

void CollisionKernel::applyScattering(Particle& obj) 
{
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

CollisionKernel::~CollisionKernel() {}

AnalogCollisionKernel::AnalogCollisionKernel() : CollisionKernel() {}

void AnalogCollisionKernel::collisionNeutron(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    double rnd = rng.R();
    int g = _buffer[n].Eg();

    const Material& mat = geometry_ptr->getRegion(_buffer[n].in_region()).material();

    if (rnd < mat.E_c(current_reactivity)[g]/mat.E_t(current_reactivity)[g]) // capture a neutron
        _buffer[n].kill();
    else if (rnd < (mat.E_c(current_reactivity)[g] + mat.E_f()[g])/mat.E_t(current_reactivity)[g]) 
        resolveFission(_buffer,n,current_reactivity);
    else // isotropic and elastic scattering
        applyScattering(_buffer[n]);
}

void AnalogCollisionKernel::resolveFission(vector<Particle>& _buffer, int  n, size_t current_reactivity)
{
    size_t g = _buffer[n].Eg();
    
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    
    int yield = floor(rng.R() + mat.nu_p()[g]);
    for (int i = 0; i < yield; i++)
    {
        int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        double mu = rng.Rmu();
        double phi = rng.Rphi();
        double sin_theta = sqrt(max(0.,1.-mu*mu));
        _buffer.emplace_back(_buffer[n].r(), 
                             Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                             new_g, settings::v[new_g], _buffer[n].t(), _buffer[n].wgt());
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
    }

    int delayed_yield = floor(rng.R() + mat.beta_tot()*mat.nu()[g]);
    for (int i = 0; i < delayed_yield; i++)
    {
        vector<double> f;
        for (const auto& v: mat.beta())
            f.push_back(v/mat.beta_tot());

        int fam = integer_cumulative_law(f, rng.R());
        double d_time = -log(1-rng.R())/mat.lambda()[fam];

        _buffer.emplace_back(_buffer[n].r(), _buffer[n].t(), 
                             _buffer[n].t() + d_time, 
                             _buffer[n].wgt()*settings::IR, fam);
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
    }
    _buffer[n].kill();
} //not up-to-date

AnalogCollisionKernel::~AnalogCollisionKernel() {}

BasicCollisionKernel::BasicCollisionKernel() : CollisionKernel() {}

void BasicCollisionKernel::collisionNeutron(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    int g = _buffer[n].Eg();

    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();

    _buffer[n].multiplyWeight(mat.E_s()[g]/mat.E_t(current_reactivity)[g]);
    resolveFission(_buffer, n,current_reactivity);
    applyScattering(_buffer[n]);
}

void BasicCollisionKernel::resolveFission(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    int g = _buffer[n].Eg();

    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    
    int yield = floor(rng.R() + mat.nu_p()[g]*mat.E_f()[g]/mat.E_t(current_reactivity)[g]);
    for (int i = 0; i < yield; i++)
    {
        int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        double mu = rng.Rmu();
        double phi = rng.Rphi();
        double sin_theta = sqrt(max(0.,1.-mu*mu));

        _buffer.emplace_back(_buffer[n].r(), 
                             Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                             new_g, mat.v()[new_g], _buffer[n].t(), _buffer[n].wgt());
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
        _buffer.back().roulette(rng.R());
        splitting(_buffer, _buffer.size() - 1);
    }

    int delayed_yield = floor(rng.R()
            + (mat.nu()[g] - mat.nu_p()[g])*mat.E_f()[g]/mat.E_t(current_reactivity)[g]);
    for (int i = 0; i < delayed_yield; i++)
    {
        vector<double> f;
        for (const auto& v: mat.beta())
            f.push_back(v/mat.beta_tot());

        int fam = integer_cumulative_law(f, rng.R());
        double d_time = -log(1-rng.R())/mat.lambda()[fam];

        _buffer.emplace_back(_buffer[n].r(), _buffer[n].t(), 
                             _buffer[n].t() + d_time, 
                             _buffer[n].wgt()*settings::IR,fam);
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
    }
} 

BasicCollisionKernel::~BasicCollisionKernel() {}

BranchlessCollisionKernel::BranchlessCollisionKernel() : CollisionKernel() {}

void BranchlessCollisionKernel::collisionNeutron(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    int g = _buffer[n].Eg();

    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();

    _buffer[n].multiplyWeight((mat.E_s()[g] + mat.nu()[g]*mat.E_f()[g])/
                            mat.E_t(current_reactivity)[g]);
    if (rng.R()*(mat.E_s()[g] + mat.nu()[g]*mat.E_f()[g]) < mat.E_s()[g])
        applyScattering(_buffer[n]);
    else 
        resolveFission(_buffer, n,current_reactivity);
}

void BranchlessCollisionKernel::resolveFission(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    size_t g = _buffer[n].Eg();
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));

    if (rng.R()*mat.nu()[g] < mat.nu_p()[g])
    {
        int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
        _buffer[n].setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
        _buffer[n].setEg(new_g);
        _buffer[n].setV(mat.v()[new_g]);
    }
    else
    {
        int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());

        _buffer[n].kill();
        double d_time = -log(1-rng.R())/mat.lambda()[fam];
        _buffer.emplace_back(_buffer[n].r(), _buffer[n].t(), 
                             _buffer[n].t() + d_time, _buffer[n].wgt()*settings::IR,fam);
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
    }
}

BranchlessCollisionKernel::~BranchlessCollisionKernel() {}

PromptBranchlessCollisionKernel::PromptBranchlessCollisionKernel() : 
    CollisionKernel() {}

void PromptBranchlessCollisionKernel::collisionNeutron(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    int g = _buffer[n].Eg();

    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    double fission_yield = mat.nu_p()[g]; 

    int delayed_yield = floor(rng.R() + mat.beta_tot()*mat.nu()[g]*mat.E_f()[g]/mat.E_t(current_reactivity)[g]);

    _buffer[n].multiplyWeight((mat.E_s()[g] + fission_yield*mat.E_f()[g])/
                               mat.E_t(current_reactivity)[g]);

    for (int i = 0; i < delayed_yield; i++)
    {
        int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
        double d_time = -log(1-rng.R())/mat.lambda()[fam];

        
        _buffer.emplace_back(_buffer[n].r(), _buffer[n].t(), 
                             _buffer[n].t() + d_time, _buffer[n].wgt()*settings::IR, fam);
        _buffer.back().setIn_region(_buffer[n].in_region());
        _buffer.back().setTree(_buffer[n].tree());
    }

    if (rng.R()*(mat.E_s()[g] + fission_yield*mat.E_f()[g]) < mat.E_s()[g])
        applyScattering(_buffer[n]);
    else 
        resolveFission(_buffer, n,current_reactivity);

}

void PromptBranchlessCollisionKernel::resolveFission(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));


    _buffer[n].setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
    _buffer[n].setEg(new_g);
    _buffer[n].setV(mat.v()[new_g]);

}

PromptBranchlessCollisionKernel::~PromptBranchlessCollisionKernel() {}

AltBranchless::AltBranchless() :
    CollisionKernel() {}

void AltBranchless::collisionNeutron(vector<Particle>& _buffer, int n, size_t current_reactivity)
{
    int g = _buffer[n].Eg();
    double rd = rng.R();

    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();

    double pb_decay = 1;
    for (size_t i = 0; i < settings::families; i++)
        if (mat.beta_tot() > 0)
            pb_decay -= mat.beta()[i]/mat.beta_tot()*exp(-mat.lambda()[i]*(settings::t_size - _buffer[n].t()));
    double tot_yield = mat.E_s()[g] + (mat.nu_p()[g] + mat.beta_tot()*mat.nu()[g]*pb_decay)*mat.E_f()[g];
    
    _buffer[n].multiplyWeight(tot_yield/mat.E_t(current_reactivity)[g]);

    if (rd*tot_yield < mat.E_s()[g])
        applyScattering(_buffer[n]);
    else if (rd*tot_yield < mat.E_s()[g] + mat.nu_p()[g]*mat.E_f()[g])
        resolveFission(_buffer, n,current_reactivity);
    else
        createPrecu(_buffer, n);
}

void AltBranchless::resolveFission(vector<Particle>& _buffer, int n, size_t current_reactivity) 
{
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    int new_g = integer_cumulative_law(mat.xi_p(), rng.R());
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));

    _buffer[n].setDir(Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu));
    _buffer[n].setEg(new_g);
    _buffer[n].setV(mat.v()[new_g]);
}

void AltBranchless::createPrecu(vector<Particle>& _buffer, int n)
{
    const Material& mat = 
        geometry_ptr->getRegion(_buffer[n].in_region()).material();
    
    int fam = integer_cumulative_law(mat.beta(), rng.R()*mat.beta_tot());
    double pb_decay = 1-exp(-mat.lambda()[fam]*(settings::t_size - _buffer[n].t()));

    double max_time = settings::t_size - _buffer[n].t();
    double rd = rng.R();
    double d_time = settings::t_size 
        + log(rd + exp(-mat.lambda()[fam]*max_time)*(1-rd))/mat.lambda()[fam];

    _buffer.emplace_back(_buffer[n].r(), _buffer[n].t(), 
                     d_time, _buffer[n].wgt()*settings::IR*pb_decay, fam);
    _buffer.back().setIn_region(_buffer[n].in_region());
    _buffer.back().setTree(_buffer[n].tree());
    
    _buffer[n].kill();
}

AltBranchless::~AltBranchless() {}
