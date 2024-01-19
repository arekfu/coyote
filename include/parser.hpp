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

#ifndef __PARSER__
#define __PARSER__

#include <string>
#include <memory>
#include "material_class.hpp"
#include "region_class.hpp"
#include "settings.hpp"
#include "geometry_class.hpp"
#include "source_class.hpp"
#include <yaml-cpp/yaml.h>

extern std::vector<std::shared_ptr<Material>> materials;
extern std::shared_ptr<Geometry> geometry;
extern std::vector<std::shared_ptr<Source>> sources;

void parseInput();
void makeSettings(const YAML::Node& setnode);
void makeMaterials(const YAML::Node& matsnode);
void makeGeometry(const YAML::Node& regionsnode);
void makeTallies(const YAML::Node& talliesnode);
void makeSources(const YAML::Node& sourcenode);
const YAML::Node parseGeometry();

#endif
