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

#include <algorithm>
#include "population_control_class.hpp"
#include "omp.h"
#include "settings.hpp"
#include "error.hpp"

using namespace std;

PopulationControl::PopulationControl() : 
    rng(), t(0)
{
    #ifdef _OPENMP
        rng.set_stream(omp_get_thread_num() + 1000);
    #else
        rng.set_stream(214564);
    #endif
    
    R = uniform_real_distribution<double>(0,1.);
}

PopulationControl::~PopulationControl() {}

NoControl::NoControl() : PopulationControl() {}
void NoControl::applyPopulationControl(vector<Particle>& _buffer) 
{
    [[maybe_unused]] bool b = _buffer.empty();
    assert(b); //just to get rid of unused variable warning...
}
NoControl::~NoControl() {}

SplittingRoulette::SplittingRoulette() : PopulationControl(), 
    russian_threshold(0.8), splitting_threshold(2.), base_weight(1.) {}

void SplittingRoulette::applyPopulationControl(vector<Particle>& _buffer)
{
    russianRoulette(_buffer);
    splitting(_buffer);
}

void SplittingRoulette::russianRoulette(vector<Particle>& _buffer)
{
    for (auto& obj: _buffer)
        if (obj.wgt() < russian_threshold) {
            if (R(rng)*base_weight < obj.wgt())
                obj.setWgt(base_weight);
            else
                obj.kill();
        }
    //without collapsed precu but forced fission
//    if (settings::decay_mode == settings::DecayMode::FORCED_DECAY) {
//        const Material& mat = 
//            geometry_ptr->getRegion(obj.in_region()).material();
//
//        double t_end = t + settings::dt;
//        for (auto& obj: _buffer) {
//            if (obj.is_precu()) {
//                size_t f = obj.family();
//                double delayed_wgt = obj.wgt()*exp(mat.lambda()[f]*obj.b_t())*(exp(-mat.lambda()[f]*t)
//                                        - exp(-mat.lambda()[f]*t_end));
//                if (delayed_wgt < russian_threshold) {
//                    if (R(rng)*base_weight < delayed_wgt)
//                        obj.multiplyWeight(base_weight/delayed_wgt);
//                    else
//                        obj.kill();
//                }
//            }
//            else
//                if (obj.wgt() < russian_threshold) {
//                    if (R(rng)*base_weight < obj.wgt())
//                        obj.setWgt(base_weight);
//                    else
//                        obj.kill();
//                }
//        }
//        t = t_end;
//    }
//    else  {
//        for (auto& obj: _buffer) {
//        if (obj.is_alive())
//            if (obj.wgt() < russian_threshold) {
//                if (R(rng)*base_weight < obj.wgt())
//                    obj.setWgt(base_weight);
//                else
//                    obj.kill();
//            }
//        }
//    }
}

void SplittingRoulette::splitting(vector<Particle>& _buffer)
{
    for (size_t i = 0; i < _buffer.size(); i++)
        if (_buffer[i].wgt() > splitting_threshold)
        {
            _buffer[i].multiplyWeight(0.5);
            _buffer.push_back(_buffer[i]);
        }
}

SplittingRoulette::~SplittingRoulette() {}

WeightCombing::WeightCombing() : 
    PopulationControl()
{
    n_buffer.reserve(settings::n_particles);
    p_buffer.reserve(settings::n_particles);
}

void WeightCombing::applyPopulationControl(vector<Particle>& _buffer)
{
    double wgt = 0;
    size_t num;
    double av_wgt, comb_pos, current_particle;

    for (const auto& obj: _buffer) {
        if (obj.is_precu())
            wgt += obj.exp_wgt();
        else 
            wgt += obj.wgt();
    }

    n_buffer.clear();
    shuffle(_buffer.begin(), _buffer.end(), rng);
    num = round(wgt);
    

    av_wgt = wgt/static_cast<double>(num);
    current_particle = 0;
    comb_pos = R(rng)*av_wgt;

    if (not _buffer.empty())
        for (const auto& obj: _buffer)
        {
            if (obj.is_precu()) {
                current_particle += obj.exp_wgt();
                while (comb_pos < current_particle)
                {
                    n_buffer.emplace_back(obj);
                    n_buffer.back().setWgt(av_wgt*obj.wgt()/obj.exp_wgt());
                    n_buffer.back().setExpWgt(av_wgt);
                    comb_pos += av_wgt;
                }
            }
            else {
                current_particle += obj.wgt();
                while (comb_pos < current_particle)
                {
                    n_buffer.emplace_back(obj);
                    n_buffer.back().setWgt(av_wgt);
                    comb_pos += av_wgt;
                }
            }
        }

    swap(n_buffer, _buffer);
    n_buffer.clear();
    n_buffer.reserve(_buffer.size());
}

WeightCombing::~WeightCombing() {}
