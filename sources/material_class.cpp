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

#include "material_class.hpp"
#include "settings.hpp"
#include "ressources.hpp"
#include "error.hpp"
#include <iostream>

using namespace std;

Material::Material(const YAML::Node& matnode) : m_beta_tot(0), m_lam_bar(0)
{
    if (matnode["a-par"] and matnode["a-par"].IsScalar()) {
        m_a_par = matnode["a-par"].as<double>();
        if (settings::reactivity_change and m_a_par != 1.) 
            cout << "Caution: don't use a and b parameters with reactivity change" << endl;
    }
    else
        m_a_par = 1;
    
    if (matnode["speed"] and matnode["speed"].size() == settings::groups)
        m_v = matnode["speed"].as<vector<double>>();
    else
        fatal_error("no speed field provided", __FILE__, __LINE__);

    if (matnode["b-par"] and matnode["b-par"].IsScalar()) {
        m_b_par = matnode["b-par"].as<double>();
        if (settings::reactivity_change and m_b_par != 0.) 
            cout << "Caution: don't use a and b parameters with reactivity change" << endl;
    }
    else
        m_b_par = 0.;

    if (matnode["name"])
        m_name = matnode["name"].as<string>();
    
    if (matnode["nu"] and matnode["nu"].size() == settings::groups)
        m_nu = matnode["nu"].as<vector<double>>();
    else
        fatal_error("no nu field provided", __FILE__, __LINE__);

    if (matnode["fission"] and matnode["fission"].size() == settings::groups)
        m_E_f = matnode["fission"].as<vector<double>>();
    else 
        fatal_error("no fission field provided", __FILE__, __LINE__);
    
    if (not settings::reactivity_change) {
        m_E_c = vector<vector<double>>(1, vector<double>(settings::groups, 0));
        if (matnode["capture"] and matnode["capture"].size() == settings::groups)
            m_E_c[0] = matnode["capture"].as<vector<double>>();
        else
            fatal_error("no valid capture field provided", __FILE__, __LINE__);
    }
    else {
        m_E_c = 
            vector<vector<double>>(settings::reactivity_time.size()+1, vector<double>(settings::groups, 0));
        if (matnode["capture"].IsSequence() and matnode["capture"].size() == settings::reactivity_time.size()+1) {
            for (size_t i = 0; i < settings::reactivity_time.size()+1; i++){
                if (matnode["capture"][i].IsSequence() and matnode["capture"][i].size() == settings::groups)
                    m_E_c[i] = matnode["capture"][i].as<vector<double>>();
                else
                    fatal_error("invalid group number for capture cross sections", __FILE__, __LINE__);
            }
        }
        else
            fatal_error("invalid reactivity change number for capture cross sections", __FILE__, __LINE__);
    }

    for (size_t i = 0; i < m_E_c[0].size(); i++)
        m_E_c[0][i] = m_a_par*m_E_c[0][i] + m_b_par/m_v[i];

    if (matnode["scattering"] and matnode["scattering"].size() == settings::groups)
        m_E_s = matnode["scattering"].as<vector<double>>();
    else
        fatal_error("no scattering field provided", __FILE__, __LINE__);

    m_s_mat = 
        vector<vector<double>>(settings::groups, vector<double>(settings::groups,0));
    if (matnode["scattering-mat"] and matnode["scattering-mat"].IsSequence())
        for (size_t i = 0; i < settings::groups; i++)
        {
            if (matnode["scattering-mat"] and matnode["scattering-mat"].IsSequence())
                m_s_mat[i] = matnode["scattering-mat"][i].as<vector<double>>();
            else
                fatal_error("no scattering-mat field", __FILE__, __LINE__);
        }
    else
        fatal_error("no scattering-mat field provided", __FILE__, __LINE__);

    if (matnode["xi-p"] and matnode["xi-p"].size() == settings::groups)
        m_xi_p = matnode["xi-p"].as<vector<double>>();
    else 
        fatal_error("no xi-p field provided", __FILE__, __LINE__);


    if (matnode["lambda"] and matnode["lambda"].size() == settings::families)
        m_lambda = matnode["lambda"].as<vector<double>>();
    else
        fatal_error("no lambda field provided", __FILE__, __LINE__);

    if (matnode["beta"] and matnode["beta"].size() == settings::families)
        m_beta = matnode["beta"].as<vector<double>>();
    else
        fatal_error("no beta field provided", __FILE__, __LINE__);

    m_xi_d = vector<vector<double>>(settings::families, 
                                    vector<double>(settings::groups, 0));
    if (matnode["xi-d"] and matnode["xi-d"].size() == settings::families)
        for (size_t i = 0; i < settings::families; i++)
        {
            if (matnode["xi-d"][i] and matnode["xi-d"][i].size() == settings::groups)
                m_xi_d[i] = matnode["xi-d"][i].as<vector<double>>();
            else
                fatal_error("no xi-d field provided", __FILE__, __LINE__);
        }
    else
        fatal_error("no xi-d field provided", __FILE__, __LINE__);

    for (const double& value : m_beta)
        m_beta_tot += value;

    if (settings::reactivity_change) {
        m_E_t = 
            vector<vector<double>>(settings::reactivity_time.size()+1, vector<double>(settings::groups, 0));
        for (size_t t = 0; t < settings::reactivity_time.size()+1; t++)
            for (size_t i = 0; i < settings::groups; i++)
                m_E_t[t][i] = m_E_f[i] + m_E_c[t][i] + m_E_s[i];
    }
    else {
        m_E_t = 
            vector<vector<double>>(1, vector<double>(settings::groups, 0));
        for (size_t i = 0; i < settings::groups; i++)
            m_E_t[0][i] = m_E_f[i] + m_E_c[0][i] + m_E_s[i];
    }

    for (size_t i = 0; i < m_nu.size(); i++)
        m_nu_p.push_back(m_nu[i]*(1-m_beta_tot));

    double tmp = 0;
    for (size_t i = 0; i < settings::families; i++)
        tmp += m_beta[i]/m_lambda[i];
    if (tmp > 0)
        m_lam_bar = m_beta_tot/tmp;
}

const vector<double>& Material::E_f() const {return m_E_f;}
const vector<double>& Material::E_c(const size_t time) const {return m_E_c[time];}
const vector<double>& Material::E_s() const {return m_E_s;}
const vector<double>& Material::E_t(const size_t time) const {return m_E_t[time];}
const vector<double>& Material::beta() const {return m_beta;}
const vector<double>& Material::nu() const {return m_nu;}
const vector<double>& Material::nu_p() const {return m_nu_p;}
const vector<double>& Material::lambda() const {return m_lambda;}
const vector<double>& Material::xi_p() const {return m_xi_p;}
const vector<double>& Material::v() const {return m_v;}
const vector<vector<double>>& Material::s_mat() const {return m_s_mat;}
const vector<vector<double>>& Material::xi_d() const {return m_xi_d;}
double Material::beta_tot() const {return m_beta_tot;}
double Material::lam_bar() const {return m_lam_bar;}
double Material::a_par() const {return m_a_par;}
double Material::b_par() const {return m_b_par;}
string Material::name() const {return m_name;}
