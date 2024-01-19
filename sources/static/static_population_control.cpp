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

#include "static_population_control.hpp"
#include "settings.hpp"
#include "ressources.hpp"
#include "omp.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <iterator>

using namespace std;

StaticControl::StaticControl() : rng(254000 + omp_get_num_threads()) {} 

StaticNoControl::StaticNoControl() : StaticControl() {}
void StaticNoControl::applyPopulationControl(vector<Particle>& _buffer) {
    [[maybe_unused]] bool b = _buffer.empty();
}

StaticCombing::StaticCombing() : StaticControl() {mCombedBuffer.reserve(settings::n_particles);}
void StaticCombing::applyPopulationControl(vector<Particle>& _buffer)
{
    uniform_real_distribution<double> R(0,1.);
    double total_weight = 0, av_weight;
    double comb_pos;
    double current_particle = 0;
    mCombedBuffer.clear();
    shuffle(_buffer.begin(), _buffer.end(), rng);

    for (const auto& obj: _buffer)
        if (obj.is_alive())
            total_weight += obj.wgt();

    av_weight = total_weight/settings::n_particles;
    comb_pos = R(rng)*av_weight;
    
    for (const auto& obj: _buffer)
    {
        if (obj.is_alive())
        {
            current_particle += obj.wgt();
            while (comb_pos < current_particle)
            {
                mCombedBuffer.push_back(obj);
                mCombedBuffer.back().setWgt(av_weight);
                comb_pos += av_weight;
            }
        }
    }
    swap(_buffer, mCombedBuffer);
}

StaticWOR::StaticWOR() : StaticControl() {mNextBuffer.reserve(settings::n_particles);}
void StaticWOR::applyPopulationControl(vector<Particle>& _buffer) {

    mNextBuffer.clear();
    mNextBuffer.reserve(settings::n_particles);
    shuffle(_buffer.begin(), _buffer.end(), rng);
    size_t to_sample = static_cast<size_t>(settings::n_particles);

    if (to_sample > _buffer.size()) {
        const size_t to_copy = static_cast<size_t>(
                floor(static_cast<double>(to_sample) / static_cast<double>(_buffer.size())));
        for (const auto& obj : _buffer)

            for (size_t i = 0; i < to_copy; i++) {
                mNextBuffer.emplace_back(obj);
//                mNextBuffer.back().setWgt(1.);
            }
        to_sample -= _buffer.size()*to_copy;
    }
    
//    goood ! 
    sample(_buffer.begin(), _buffer.end(), back_inserter(mNextBuffer), to_sample, rng);
    swap(_buffer, mNextBuffer);


//      baaaad !
//    vector<double> wgt(_buffer.size(), 0.);
//    for (size_t i = 0; i < _buffer.size(); i++)
//        wgt[i] = _buffer[i].wgt();
//
//    for (size_t i = 0; i < to_sample; i++) {
//        discrete_distribution<size_t> R(wgt.begin(), wgt.end());
//        size_t id = R(rng);
//        mNextBuffer.emplace_back(_buffer[id]);
//        mNextBuffer.back().setWgt(1.);
//        mNextBuffer.back().setTree(_buffer[id].tree());
//        mNextBuffer.back().setB_r(mNextBuffer.back().r());
//        wgt[id] = 0.;
//    }
//    swap(_buffer, mNextBuffer);
}
