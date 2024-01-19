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

#include <iostream>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "geometry_class.hpp"
#include "settings.hpp"
#include "error.hpp"

using namespace std;

Geometry::Geometry() : regions(), m_size(0), b_cross(0) {}

Geometry::Geometry(const YAML::Node& regionsnode)
{
    regions.clear();
    for (int i = 0; i < settings::regions; i++)
    {
        regions.emplace_back(regionsnode[i]);
    }
}

int Geometry::getRegionIndex(const Particle& obj) const
{
    return getRegionIndex(obj.r(), obj.dir());
}

int Geometry::getRegionIndex(const Vector& r, const Vector& dir) const
{
    int pos = 0;
    for (const auto& reg : regions)
    {
        if (r.x() >= reg.start() and r.x() <= reg.end())
            return pos;
        pos++;
    }
    cout << "r = " << r << endl;
    cout << "dir = " << dir << endl;
    throw std::invalid_argument("particle not in any region");
}

const Region& Geometry::getRegion(const int& _reg) const {return regions[_reg];}

void Geometry::noCrossing(Particle& obj, double t_dist, double t_time) const
{
    obj.setR(obj.r() + obj.dir()*t_dist);
    obj.setT(obj.t() + t_time); //synchronise neutrons to t_bin
}

double Geometry::distToFirstSurface(Particle& obj) const
{
    int reg = obj.in_region();
    if (obj.dir().x() < 0)
    {
        obj.setBoundary_cross(regions[reg].l_bound());
        return abs((regions[reg].start() - obj.r().x())/obj.dir().x());
    }
    else
    {
        obj.setBoundary_cross(regions[reg].r_bound());
        return abs((regions[reg].end() - obj.r().x())/obj.dir().x());
    }
}

void Geometry::crossSurface(Particle& obj, const double& t_dist) const
{
    double t_time = t_dist/obj.v();
    obj.setT(obj.t() + t_time);
    goToSurface(obj);

    if (obj.boundary_cross())
        applyBoundaryConditions(obj);
    obj.setBoundary_cross(false);
}


void Geometry::goToSurface(Particle& obj) const
{
    int reg = obj.in_region();
    double x_s = regions[reg].start(), x_e = regions[reg].end();
    if (obj.dir().signX() < 0)
    {
        obj.setR(Vector(x_s, 
                        obj.r().y() + (x_s-obj.r().x())*obj.dir().y()/obj.dir().x(),
                        obj.r().z() + (x_s-obj.r().x())*obj.dir().z()/obj.dir().x()));
        if (obj.in_region() - 1 >= 0)
            obj.setIn_region(obj.in_region() - 1);
    }
    else
    {
        obj.setR(Vector(x_e, 
                        obj.r().y() + (x_e-obj.r().x())*obj.dir().y()/obj.dir().x(),
                        obj.r().z() + (x_e-obj.r().x())*obj.dir().z()/obj.dir().x()));
        if (obj.in_region() + 1 < (int) regions.size())
            obj.setIn_region(obj.in_region() + 1);
    }
}

Geometry::~Geometry() {}

GeometryRef::GeometryRef(const YAML::Node& regionsnode) : Geometry(regionsnode) {}

void GeometryRef::applyBoundaryConditions(Particle& obj) const
{
    Vector n(0,0,0);
    int reg = obj.in_region();

    if (obj.dir().x() < 0)
    {
        obj.setB_r(Vector(2*regions[reg].start()-obj.b_r().x(), 
                            obj.b_r().y(), obj.b_r().z()));
        n = Vector(1., 0, 0);
    }
    else
    {
        obj.setB_r(Vector(2*regions[reg].end()-obj.b_r().x(), 
                            obj.b_r().y(), obj.b_r().z()));
        n = Vector(-1., 0, 0);
    }
    obj.setDir((obj.dir() - 2.*(obj.dir()*n)*n).unitary());
}
GeometryRef::~GeometryRef() {}

GeometryPer::GeometryPer(const YAML::Node& regionsnode) : Geometry(regionsnode) {}

void GeometryPer::applyBoundaryConditions(Particle& obj) const
{
    if (obj.in_region() == 0)
        obj.setNext_region(regions.size()-1);
    else
        obj.setNext_region(0);
    obj.setR(Vector(-obj.r().x(), obj.r().y(), obj.r().z()));
}
GeometryPer::~GeometryPer() {}

GeometryLeak::GeometryLeak(const YAML::Node& regionsnode) : Geometry(regionsnode) {}

void GeometryLeak::applyBoundaryConditions(Particle& obj) const
{ 
    obj.leakage();
    obj.kill();
}

GeometryLeak::~GeometryLeak() {}
