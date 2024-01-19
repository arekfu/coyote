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
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "ressources.hpp"
#include "settings.hpp"
#include "error.hpp"

using namespace std;

size_t integer_cumulative_law(const vector<double>& f, double r)
{
    double tmp = 0;
    for (size_t k = 0; k < f.size(); k++)
    {
        tmp += f[k];
        if (r < tmp)
            return k;
    }
    cout << "tmp = " << tmp << endl;
    throw invalid_argument("not a vector<double> cumulative law");
}

size_t integer_cumulative_law(const NDArray<double>& f, int ind, double r)
{
    double tmp = 0;
    for (size_t k = 0; k < f.shape()[1]; k++)
    {
        tmp += f(ind, k);
        if (r < tmp)
            return k;
    }
    throw invalid_argument("not a NDArray<double> cumulative law");
}

size_t integer_cumulative_law(const vector<vector<double>>& f, int ind, double r)
{
    double tmp = 0;
    for (size_t k = 0; k < f[ind].size(); k++)
    {
        tmp += f[ind][k];
        if (r < tmp)
            return k;
    }
    throw invalid_argument("not a vector<vector<double>> cumulative law");
}

double modulo(const double& a, const double& b)
{
    double q = 0;
    if (a == 0)
        return 0.;
    else if (a > 0)
        while ((q+1)*b < a)
            q++;
    else
        while (q*b > a)
            q--;

    return a - q*b;
}

vector<double> read_param(string name)
{
    ifstream par(name);
    vector<double> tmp;
    
    if (par)
    {
        string s;
        while(not par.eof())
        {
            par >> s;
            if (s.find_first_of("0123456789") != string::npos and not par.eof())
                tmp.push_back(stod(s));
        }
        par.close();
        return tmp;
    }
    else
        throw std::invalid_argument("wrong parameter name");
}

int factorial(const int& num)
{
    int temp = 1;
    if (num > 0)
    {
        for (int i = 1; i <= num; i++)
            temp *= i;
        return temp;
    }
    else if (num == 0)
        return 1;
    else
    {
        cout << "Factorial of a negative number does not exist" << endl;
        exit(0);
    }
}

void splitting(vector<Particle>& buffer, const int& n) {
    if (buffer[n].wgt() > settings::splitting_value and buffer[n].is_alive())
    {
        buffer[n].multiplyWeight(1./settings::splitting_value);
        for (int i = 0; i < settings::splitting_value-1; i++) 
            buffer.push_back(buffer[n]);
    }
}

void checkIndicies(int& bin_i, int& bin_f)
{
    if (bin_i >= settings::x_bin)
        bin_i = settings::x_bin - 1;
    else if (bin_i < 0)
        bin_i = 0;

    if (bin_f >= settings::x_bin)
        bin_f = settings::x_bin - 1;
    else if (bin_f < 0)
        bin_f = 0;
}

void checkIndex(int& bin)
{
    if (bin >= settings::x_bin)
        bin  = settings::x_bin - 1;
    else if (bin < 0)
        bin = 0;
}
