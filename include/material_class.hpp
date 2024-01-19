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

#ifndef __MATERIAL_CLASS__
#define __MATERIAL_CLASS__

#include <vector>
#include <string>
#include "ndarray.hpp"
#include "yaml-cpp/yaml.h"

class Material
{
    public:
        Material() = default;
        Material(const YAML::Node& matnode);
        const std::vector<double>& E_f() const; 
        const std::vector<double>& E_c(const size_t time = 0) const; 
        const std::vector<double>& E_s() const; 
        const std::vector<double>& E_t(const size_t time = 0) const; 
        const std::vector<double>& beta() const; 
        const std::vector<double>& nu() const; 
        const std::vector<double>& nu_p() const; 
        const std::vector<double>& lambda() const; 
        const std::vector<double>& xi_p() const; 
        const std::vector<double>& v() const; 
        const std::vector<std::vector<double>>& s_mat() const;
        const std::vector<std::vector<double>>& xi_d() const;
        double beta_tot() const;
        double lam_bar() const;
        double a_par() const;
        double b_par() const;
        std::string name() const;

    private:
        std::string m_name;
        std::vector<double> m_E_f;
        std::vector<std::vector<double>> m_E_c;
        std::vector<double> m_E_s;
        std::vector<std::vector<double>> m_E_t;
        std::vector<double> m_beta;
        std::vector<double> m_nu;
        std::vector<double> m_nu_p;
        std::vector<double> m_lambda;
        std::vector<double> m_xi_p;
        std::vector<double> m_v;
        std::vector<std::vector<double>> m_s_mat;
        std::vector<std::vector<double>> m_xi_d;
        double m_a_par, m_b_par;
        double m_beta_tot;
        double m_lam_bar;
};

#endif
