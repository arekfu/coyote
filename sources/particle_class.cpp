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

#include "particle_class.hpp"

using namespace std;

Particle::Particle(const Vector& _r) : m_r(_r), m_b_r(_r), m_wgt(1.) {}

//for neutrons in PI
Particle::Particle(const Vector& _r, const Vector& _dir, 
                   int _Eg, double _v, double _wgt):
    m_r(_r), m_b_r(_r), m_dir(_dir), m_v(_v), m_wgt(_wgt), m_t(0), m_b_t(0), 
    m_d_t(0), m_is_precu(false), m_is_alive(true), m_Eg(_Eg), m_leak(false),
    m_boundary_cross(false), m_tree(0) {}

//for neutrons in tarapiaTapiocos
Particle::Particle(const Vector& _r, const Vector& _dir,
                   int _Eg, double _v, double _t, double _wgt) :
    m_r(_r), m_b_r(_r), m_dir(_dir), m_v(_v), m_wgt(_wgt), m_t(_t), m_b_t(_t), 
    m_d_t(0), m_is_precu(false), m_is_alive(true), m_Eg(_Eg), m_leak(false),
    m_boundary_cross(false), m_tree(0) {}

//for precusors
Particle::Particle(const Vector& _r, double _t, double _d_t, 
                   double _wgt, size_t _family, bool from_eq) :
    m_r(_r), m_b_r(_r), m_dir(), m_wgt(_wgt), m_exp_wgt(_wgt), m_t(_t), m_b_t(_t), m_d_t(_d_t),
    m_is_precu(true), m_is_alive(true), m_Eg(0), m_family(_family), m_leak(false),
    m_boundary_cross(false), m_tree(0), m_from_eq(from_eq) {}

const Vector& Particle::r() const {return m_r;}
const Vector& Particle::b_r() const {return m_b_r;}
const Vector& Particle::dir() const {return m_dir;}
double Particle::v() const {return m_v;}
double Particle::wgt() const {return m_wgt;}
double Particle::exp_wgt() const {return m_exp_wgt;}
double Particle::t() const {return m_t;}
double Particle::b_t() const {return m_b_t;}
double Particle::d_t() const {return m_d_t;}
bool Particle::is_precu() const {return m_is_precu;}
bool Particle::is_alive() const {return m_is_alive;}
bool Particle::has_leaked() const {return m_leak;}
int Particle::Eg() const {return m_Eg;}
int Particle::family() const {return m_family;}
int Particle::tree() const {return m_tree;}
int Particle::in_region() const {return m_in_region;}
int Particle::next_region() const {return m_next_region;}
bool Particle::boundary_cross() const {return m_boundary_cross;}
bool Particle::from_eq() const {return m_from_eq;}
void Particle::setBoundary_cross(bool b_cross) {m_boundary_cross = b_cross;}
void Particle::setR(const Vector& _r) {m_r = _r;}
void Particle::setB_r(const Vector& _b_r) {m_b_r = _b_r;}
void Particle::setT(double _t) {m_t = _t;}
void Particle::setB_t(const double _t) {m_b_t = _t;}
void Particle::setDir(const Vector& _dir) {m_dir = _dir;}
void Particle::setV(double _v) {m_v = _v;}
void Particle::setWgt(const double _wgt) {m_wgt = _wgt;}
void Particle::setExpWgt(const double _wgt) {m_exp_wgt = _wgt;}
void Particle::setEg(int _Eg) {m_Eg = _Eg;}
void Particle::setFamily(int _family) {m_family = _family;}
void Particle::setTree(int _tree) {m_tree = _tree;}
void Particle::multiplyWeight(double x) {m_wgt *= x;}
void Particle::multiplyExpWeight(double x) {m_exp_wgt *= x;}
void Particle::kill() {m_is_alive = false;}
void Particle::leakage() {m_leak = true;}
void Particle::make_precu() {m_is_precu = true;}
void Particle::make_neutron() {m_is_precu = false;}
void Particle::setIn_region(int _region) 
{
    m_in_region = _region;
    m_next_region = _region;
}
void Particle::setNext_region(int _next) {m_next_region = _next;}
void Particle::update_region() {m_in_region = m_next_region;}

int Particle::sample_family(const Material& mat, const double& t, const double& r) {
    if (m_from_eq) {
        double sum = 0;
        for (size_t i = 0; i < settings::families; i++)
            sum += mat.beta()[i]*exp(-mat.lambda()[i]*(t-m_b_t));
        double tmp_sample = 0;
        for (size_t i = 0; i < settings::families; i++) {
            tmp_sample += mat.beta()[i]*exp(-mat.lambda()[i]*(t - m_b_t));
            if (r*sum < tmp_sample)
                return i;
        }
    }
    else {
        double sum = 0;
        for (size_t i = 0; i < settings::families; i++)
            sum += mat.lambda()[i]*mat.beta()[i]*exp(-mat.lambda()[i]*(t-m_b_t));
        double tmp_sample = 0;
        for (size_t i = 0; i < settings::families; i++) {
            tmp_sample += mat.lambda()[i]*mat.beta()[i]*exp(-mat.lambda()[i]*(t-m_b_t));
            if (r*sum < tmp_sample)
                return i;
        }
    }

    fatal_error("wrong sampling of family when using collapsed precu", __FILE__, __LINE__);
    return 0;
}

void Particle::roulette(const double& r) {
    double russian_threshold = settings::roulette_value, base_weight = 1.;
    if (m_wgt < russian_threshold and m_is_alive)
    {
        if (r > m_wgt/base_weight)
            m_is_alive = false;
        else
            m_wgt = base_weight;
    }
}
