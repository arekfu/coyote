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

#include "region_class.hpp"
#include "error.hpp"
#include "parser.hpp"

using namespace std;

Region::Region(const YAML::Node& regionnode)
{
    if (regionnode["start"] and regionnode["start"].IsScalar())
        m_start = regionnode["start"].as<double>();
    else
        fatal_error("no start field provided", __FILE__, __LINE__);

    if (regionnode["end"] and regionnode["end"].IsScalar())
        m_end = regionnode["end"].as<double>();
    else
        fatal_error("no end field provided", __FILE__, __LINE__);

    if (regionnode["boundary-left"] and regionnode["boundary-left"].IsScalar())
        m_l_bound = regionnode["boundary-left"].as<bool>();
    else
        fatal_error("no boundary-left field provided", __FILE__, __LINE__);

    if (regionnode["boundary-right"] and regionnode["boundary-right"].IsScalar())
        m_r_bound = regionnode["boundary-right"].as<bool>();
    else
        fatal_error("no boundary-right field provided", __FILE__, __LINE__);

    if (regionnode["id"] and regionnode["id"].IsScalar())
        m_region_num = regionnode["id"].as<int>();
    else
        fatal_error("no id field provided", __FILE__, __LINE__);

    if (regionnode["material"])
    {
        for (auto obj: materials)
            if (obj->name() == regionnode["material"].as<string>())
            {
                m_mat = *obj;
                break;
            }
    }
    else
        fatal_error("no material field provided", __FILE__, __LINE__);

}
Region::Region() : m_region_num(0) {}
int Region::getRegionNumber() const {return m_region_num;}
double Region::start() const {return m_start;}
double Region::end() const {return m_end;}
bool Region::l_bound() const {return m_l_bound;}
bool Region::r_bound() const {return m_r_bound;}
const Material& Region::material() const {return m_mat;}
Region::~Region() {}
