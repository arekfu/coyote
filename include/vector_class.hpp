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

#ifndef __VECTOR_CLASS__
#define __VECTOR_CLASS__

#include <iostream>
#include <vector>

class Vector
{
    public:
        Vector();
        Vector(double vx, double vy, double vz);
        Vector(double vx, double vy, double vz, double _norm);
        Vector(std::vector<double> _v);
        void setX(double vx);
        void setY(double vy);
        void setZ(double vz);
        double x() const;
        Vector Vx() const;
        double y() const;
        Vector Vy() const;
        double z() const;
        Vector Vz() const;
        int signX() const;
        int signY() const;
        int signZ() const;
        double norm() const;
        Vector unitary() const;
        void print() const;

    private:
        double m_x,m_y,m_z;
};

inline Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a.x() + b.x(), a.y() + b.y(), a.z() + b.z()); }

inline Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a.x() - b.x(), a.y() - b.y(), a.z() - b.z()); }

inline Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a.x()*b.x(), a.y()*b.y(), a.z()*b.z()); }

inline Vector operator*(const Vector& a, const double v) {
    return Vector(a.x()*v, a.y()*v, a.z()*v); }

inline Vector operator*(const double v, const Vector& a) {
    return Vector(a.x()*v, a.y()*v, a.z()*v); }

inline Vector operator*(const Vector& a, const int v) {
    return Vector(a.x()*v, a.y()*v, a.z()*v); }

inline std::ostream& operator<<(std::ostream& output, const Vector& a) {
    output << "(" << a.x() << "," << a.y() << "," << a.z() << ")";
    return output;
}
#endif
