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
#include <string>
#include <vector>
#include "tally.hpp"

using namespace std;

Tally::Tally() : name("nothing"), shape(), groups(2), 
    average(2, NDArray<double>()),
    variance(2, NDArray<double>()) {}
    

Tally::Tally(vector<size_t> _shape, const string& _name, const size_t _g) : 
    name(_name), shape(_shape), groups(_g),  
    average(_g, NDArray<double>(_shape)),
    variance(_g, NDArray<double>(_shape)) {}

const vector<NDArray<double>>& Tally::getAverage() const {return average;}
const vector<NDArray<double>>& Tally::getVariance() const {return variance;}
const vector<NDArray<double>>& Tally::getTmp() const {return tmp;}
const string& Tally::getName() const {return name;}
Tally& Tally::operator+=(const Tally& rhs)
{
    for (size_t i = 0; i < groups; i++)
    {
        average[i] += rhs.getAverage()[i];
        variance[i] += rhs.getVariance()[i];
        tmp[i] += rhs.getTmp()[i];
    }

    return *this;
}

void Tally::save(string path) const
{
    for (size_t i = 0; i < groups; i++)
    {
        average[i].save(path + "energy_group" + to_string(i) + "/" + name + ".npy");
        variance[i].save(path + "energy_group" + to_string(i) + "/" + name + "_var.npy");
    }
}

void Tally::reinitialize()
{
    for (size_t i = 0; i < groups; i++)
    {
        average[i].fill(0);
        variance[i].fill(0);
        tmp[i].fill(0);
    }
}

void Tally::reinitializeTmp()
{
    for (size_t i = 0; i < groups; i++)
        tmp[i].fill(0.);
}

void Tally::normalizeTmp(const double& _norm) {
    for (auto& obj: tmp)
        obj /= _norm;
}

void Tally::scoreTmp(const vector<size_t>& indices, size_t energy_group, double value) 
    {tmp[energy_group](indices) += value;}

void Tally::scoreTmp(const NDArray<double>& _temp, size_t _g) {
    tmp[_g] += _temp;
}

void Tally::normalize(double norm_av, double norm_var)
{
    for (size_t i = 0; i < groups; i++)
    {
        average[i] /= norm_av;
        variance[i] /= norm_var;
    }
}

void Tally::normalize(const vector<NDArray<double>>& norm_av,
                      const vector<NDArray<double>>& norm_var) {
    for (size_t i = 0; i < groups; i++) {
        average[i] /= norm_av[i];
        variance[i] /= norm_var[i];
    }
}

Tally::~Tally() {}

Tally0D::Tally0D(const string& _name) : Tally({1}, _name, 1) {
    tmp = vector<NDArray<double>>(1, NDArray<double>({1}));}

void Tally0D::score(const int& m) {
    double old_av = average[0][0];
    double old_var = variance[0][0];
    average[0][0] = old_av + (tmp[0][0] - old_av)/m;
    if (m > 1)
        variance[0][0] = old_var*(m-2)/(m-1) 
            + pow(tmp[0][0] - old_av,2)/m;
}

Tally0D::~Tally0D() {}

Tally1D::Tally1D() : Tally() {}
Tally1D::Tally1D(std::vector<size_t> _shape, const std::string& _name,
                 const size_t _g) : Tally(_shape, _name, _g)
{ tmp = vector<NDArray<double>>(_g, NDArray<double>(_shape));}

void Tally1D::score(const int& m)
{
    for (size_t i = 0; i < average[0].size(); i++)
        for (size_t e = 0; e < groups; e++)
        {
            double av_old = average[e][i];
            double var_old = variance[e][i];
            average[e][i] = av_old + (tmp[e][i] - av_old)/m;
            if (m > 1)
                variance[e][i] = var_old*(m-2)/(m-1) 
                    + pow(tmp[e][i] - av_old,2)/m;
        }
}

Tally1D::~Tally1D() {}

Tally2D::Tally2D() : Tally() {}
Tally2D::Tally2D(vector<size_t> _shape, const string& _name, const size_t _g) 
    : Tally(_shape, _name, _g)
{ tmp = vector<NDArray<double>>(_g, NDArray<double>(_shape));}

void Tally2D::score(const int& m)
{
    for (size_t i = 0; i < average[0].size(); i++)
        for (size_t e = 0; e < groups; e++)
        {
            double av_old = average[e][i];
            double var_old = variance[e][i];
            average[e][i] = av_old + (tmp[e][i] - av_old)/m;
            if (m > 1)
                variance[e][i] = var_old*(m-2)/(m-1) 
                    + pow(tmp[e][i] - av_old,2)/m;
        }
}

void Tally2D::score(const NDArray<double> & tmp, int m)
{
    for (size_t i = 0; i < average[0].size(); i++)
    {
        double av_old = average[0][i];
        double var_old = variance[0][i];
        average[0][i] = av_old + (tmp[i] - av_old)/m;
        if (m > 1)
            variance[0][i] = var_old*(m-2)/(m-1) 
                + pow(tmp[i] - av_old,2)/m;
    }
}

Tally2D::~Tally2D() {}

Tally4D::Tally4D() : Tally(){}

Tally4D::Tally4D(vector<size_t> _shape, const string& _name, const size_t _g)
    : Tally(_shape, _name, _g)
{ 
    _shape.pop_back();
    _shape.pop_back();
    tmp = vector<NDArray<double>>(_g, NDArray<double>(_shape)); 
}

void Tally4D::score(const int& m)
{
    for (size_t i = 0; i < shape[0]; i++)
        for (size_t j = 0; j < shape[1]; j++)
            for (size_t k = 0; k < shape[2]; k++)
                for (size_t l = 0; l < shape[3]; l++)
                    for (size_t e = 0; e < groups; e++)
                    {
                        double av_old = average[e](i,j,k,l);
                        double var_old = variance[e](i,j,k,l);
                        average[e](i,j,k,l) = av_old 
                            + (tmp[e](i,j)*tmp[e](k,l) - av_old)/m;
                        if (m > 1)
                            variance[e](i,j,k,l) = var_old*(m-2)/(m-1) 
                               + pow(tmp[e](i,j)*tmp[e](k,l) - av_old,2)/m;
                    }
}

Tally4D::~Tally4D() {}

