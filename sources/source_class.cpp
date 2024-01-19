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

#include "source_class.hpp"
#include "settings.hpp"
#include "error.hpp"

using namespace std;

Source::Source() : m_type("none"), m_number(0),
                   rng(2149 + omp_get_thread_num()) {}
Source::Source(const YAML::Node& node): rng(2149 + omp_get_thread_num()) {
    if (node["type"])
        m_type = node["type"].as<string>();
    else
        fatal_error("no type field in source field", __FILE__, __LINE__);
    if (node["number"] and node["number"].IsScalar())
        m_number = node["number"].as<int>();
    else
        fatal_error("no number field in source field", __FILE__, __LINE__);
}

const string& Source::type() const {return m_type;}
const int& Source::number() const {return m_number;}

Source::~Source() {}

DeltaSource::DeltaSource() : Source(), m_pos(0,0,0) {}
DeltaSource::DeltaSource(const YAML::Node& node) : Source(node) {
    if (node["position"] and node["position"].IsSequence())
        m_pos = Vector(node["position"].as<vector<double>>());
    else
        fatal_error("no position field in source field", __FILE__, __LINE__);
}


Particle DeltaSource::sampleSourceParticle() {
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));
    Particle obj(m_pos, Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                 0, settings::v[0], 0, 1.);
    
    return obj;
}

const Vector& DeltaSource::pos() const {return m_pos;}
DeltaSource::~DeltaSource() {}

UniformSource::UniformSource() : Source(), m_segment(), m_length(0) {}
UniformSource::UniformSource(const YAML::Node& node) : Source(node), m_segment() {
    if (node["segment"] and node["segment"].IsSequence())
        m_segment = node["segment"].as<vector<double>>();
    else {
        cout << "No segment field provided for uniform source: taking whole space."
             << endl;
        m_segment.push_back(-settings::x_size);
        m_segment.push_back(+settings::x_size);
    }
    m_length = m_segment[1] - m_segment[0];

    if (node["group"] and node["group"].IsScalar())
        m_group = node["group"].as<int>();
    else
        m_group = 0;
}

Particle UniformSource::sampleSourceParticle() {
    double mu = rng.Rmu();
    double phi = rng.Rphi();
    double sin_theta = sqrt(max(0.,1.-mu*mu));

    double x = m_segment[0] + rng.R()*m_length;
    Particle obj(Vector(x,0,0), Vector(sin_theta*cos(phi), sin_theta*sin(phi), mu), 
                 m_group, settings::v[m_group], 0, 1.);
    
    return obj;
}

UniformSource::~UniformSource() {}
