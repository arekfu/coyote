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

#ifndef __PARTICLE__
#define __PARTICLE__

#include "vector_class.hpp"
#include "material_class.hpp"
#include "settings.hpp"
#include "error.hpp"
#include <iostream>

class Particle
{
    public:
        Particle() = default;
        Particle(const Vector& _r);
        Particle(const Vector& _r, const Vector& _dir, 
                 int _Eg, double _v, double _wgt);
        /* constructor for PI particles */
        Particle(const Vector& _r, const Vector& _dir,
                 int _Eg, double _v, double _t, double _wgt);
        /* constructor for dynamic neutrons */
        Particle(const Vector& _r, double _t, double _d_t, 
                 double _wgt, size_t _family, bool from_eq = false);
        /* constructor for dynamic precu */
        void setDir(const Vector& _dir);
        void setV(double _v);
        void setR(const Vector& _r);
        void setB_r(const Vector& _b_r);
        void setT(double _t);
        void setB_t(const double _t);
        void setWgt(const double _wgt);
        void setExpWgt(const double _wgt);
        void setEg(int _Eg);
        void setFamily(int _family);
        void setIn_region(int _region);
        void setNext_region(int _next);
        void setTree(int _tree);
        void multiplyWeight(double x);
        void multiplyExpWeight(double x);
        void kill();
        void leakage();
        void make_precu();
        void make_neutron();
        void setBoundary_cross(bool b_cross);
        const Vector& r() const;
        const Vector& b_r() const;
        const Vector& dir() const;
        double v() const;
        double wgt() const;
        double exp_wgt() const;
        double t() const;
        double b_t() const;
        double d_t() const;
        bool is_precu() const;
        bool is_alive() const;
        bool has_leaked() const;
        int Eg() const;
        int family() const;
        int tree() const;
        bool from_eq() const;
        int in_region() const;
        int next_region() const;
        bool boundary_cross() const;
        void update_region();
        int sample_family(const Material& mat, const double& t, const double& r);
        void roulette(const double& r);

    private:
        Vector m_r, m_b_r; //currention pos and birth position
        Vector m_dir;
        double m_v, m_wgt, m_exp_wgt;
        double m_t, m_b_t, m_d_t; //current time, birth time and death time
        bool m_is_precu, m_is_alive;
        int m_Eg;
        int m_family;
        bool m_leak;
        bool m_boundary_cross;
        int m_in_region, m_next_region;
        int m_tree; //family tree, based on parent in the source (for NONEQUILIBRIUM source only)
        bool m_from_eq;
};

#endif
