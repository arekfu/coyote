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

#include "vector_class.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Vector::Vector() : m_x(0), m_y(0), m_z(0) {}
Vector::Vector(double vx, double vy, double vz) : 
    m_x(vx), m_y(vy), m_z(vz) {}
Vector::Vector(double vx, double vy, double vz, double _norm) :
    m_x(vx), m_y(vy), m_z(vz)
{
    double bad_norm = norm();
    m_x *= _norm/bad_norm;
    m_y *= _norm/bad_norm;
    m_z *= _norm/bad_norm;
}

Vector::Vector(vector<double> _v) : 
    m_x(_v[0]), m_y(_v[1]), m_z(_v[2]) {}

void Vector::setX(double vx) {m_x = vx;}
void Vector::setY(double vy) {m_y = vy;}
void Vector::setZ(double vz) {m_z = vz;}

double Vector::x() const {return m_x;}
double Vector::y() const {return m_y;}
double Vector::z() const {return m_z;}

Vector Vector::Vx() const {return Vector(m_x,0,0); }
Vector Vector::Vy() const {return Vector(0,m_y,0); }
Vector Vector::Vz() const {return Vector(0,0,m_z); }

int Vector::signX() const {return copysign(1., m_x);}
int Vector::signY() const {return copysign(1., m_y);}
int Vector::signZ() const {return copysign(1., m_z);}

double Vector::norm() const {return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);}
Vector Vector::unitary() const
{
    double Norm = norm();
    return Vector(m_x/Norm, m_y/Norm, m_z/Norm);
}

void Vector::print() const 
{ 
    cout << "(" << m_x << "," << m_y << "," << m_z << ")" << endl;
}
