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

#include "decay_kernel_class.hpp"
#include "omp.h"
#include "settings.hpp"
#include "error.hpp"
#include "parser.hpp"
#include <algorithm>

using namespace std;

DecayKernel::DecayKernel() : rng(2140001 + omp_get_thread_num())
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

DecayKernel::~DecayKernel() {}

AnalogDecayKernel::AnalogDecayKernel() : DecayKernel() {}

void AnalogDecayKernel::decayPrecursor(vector<Particle>& _buffer, size_t precu, double tmax)
{
    tmax *= settings::dt;
    if (_buffer[precu].d_t() < tmax)
    {
        const Material& mat = 
            geometry_ptr->getRegion(_buffer[precu].in_region()).material();
        int g = 
            integer_cumulative_law(mat.xi_d(), _buffer[precu].family(), rng.R());
        double mu = rng.Rmu();
        double phi = rng.Rphi();
        double sin_theta = sqrt(max(0.,1.-mu*mu));

        _buffer.emplace_back(_buffer[precu].r(), 
                             Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu),
                             g, mat.v()[g], _buffer[precu].d_t(), 
                             _buffer[precu].wgt()/settings::IR);

        _buffer.back().setIn_region(_buffer[precu].in_region());
        _buffer[precu].kill();
    }
}

AnalogDecayKernel::~AnalogDecayKernel() {}

ForcedDecayKernel::ForcedDecayKernel() : DecayKernel() {}

void ForcedDecayKernel::decayPrecursor(vector<Particle>& buffer, 
        size_t precu, double time_step)
{
    const Material& mat = 
        geometry_ptr->getRegion(buffer[precu].in_region()).material();

    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));


    if (settings::collapsed_precu) {
        double precu_time = buffer[precu].b_t();
        double t_min = max(precu_time, (time_step-1)*settings::dt);
        double t_max = time_step*settings::dt;
        double t_interval = t_max - t_min;
        double decay_time = t_min + rng.R()*t_interval;
        int family = buffer[precu].sample_family(mat, decay_time, rng.R());
        int g = integer_cumulative_law(mat.xi_d(), family, rng.R());

        if (buffer[precu].from_eq()) {
            double sum = 0;
            for (size_t i = 0; i < settings::families; i++)
                sum += mat.lam_bar()*mat.beta()[i]/mat.beta_tot()*exp(-mat.lambda()[i]*(decay_time - precu_time));
            double delayed_weight = t_interval*sum*buffer[precu].wgt();

            buffer.emplace_back(buffer[precu].r(), 
                                Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                                g, mat.v()[g], decay_time, delayed_weight/settings::IR);
            buffer.back().setIn_region(buffer[precu].in_region());
            buffer.back().roulette(rng.R());
        }
        else {
            double sum = 0;
            for (size_t i = 0; i < settings::families; i++)
                sum += mat.lambda()[i]*mat.beta()[i]*exp(-mat.lambda()[i]*(decay_time-precu_time));
            double delayed_weight = t_interval*sum/mat.beta_tot()*buffer[precu].wgt();

            buffer.emplace_back(buffer[precu].r(), 
                                Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                                g, mat.v()[g], decay_time, delayed_weight/settings::IR);
            buffer.back().setIn_region(buffer[precu].in_region());
            buffer.back().roulette(rng.R());
        }
        buffer[precu].setT(t_max);
    }
    else {
        double precu_time = buffer[precu].t();
        double t_min = max(precu_time, (time_step-1)*settings::dt);
        double t_max = time_step*settings::dt;
        double t_interval = t_max - t_min;
        double decay_time = t_min + rng.R()*t_interval;
        int f = buffer[precu].family();
        int g = integer_cumulative_law(mat.xi_d(), f, rng.R());
        double delayed_weight = buffer[precu].wgt()*t_interval*mat.lambda()[f]*
            exp(-mat.lambda()[f]*(decay_time - precu_time));

        buffer[precu].setT(t_max);
        buffer[precu].multiplyWeight(exp(-mat.lambda()[f]*(t_max - precu_time)));

        buffer.emplace_back(buffer[precu].r(), 
                            Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                            g, mat.v()[g], decay_time, delayed_weight/settings::IR);
        buffer.back().setIn_region(buffer[precu].in_region());
        buffer.back().roulette(rng.R());
    }
}

ForcedDecayKernel::~ForcedDecayKernel() {}
