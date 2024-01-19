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

#ifndef __tally__
#define __tally__

#include <vector>
#include <string>
#include "ndarray.hpp"

class Tally
{
    public:
        Tally();
        Tally(std::vector<size_t> _shape, const std::string& _name, 
              const size_t _g);
        virtual void score(const int& m) = 0;
        void scoreTmp(const std::vector<size_t>& indices, 
                size_t energy_group, double value);
        void scoreTmp(const NDArray<double>& _temp, size_t _g);
        void normalizeTmp(const double& _norm);
        void normalize(double norm_av, double norm_var);
        void normalize(const std::vector<NDArray<double>>& norm_av,
                       const std::vector<NDArray<double>>& norm_var);
        const std::vector<NDArray<double>>& getAverage() const;
        const std::vector<NDArray<double>>& getVariance() const;
        const std::vector<NDArray<double>>& getTmp() const;
        const std::string& getName() const;
        void save(std::string path) const;
        void reinitialize();
        void reinitializeTmp();
        Tally& operator+=(const Tally& rhs);
        virtual ~Tally();

    protected:
        std::string name; //name of the observable
        std::vector<size_t> shape;
        size_t groups; //number of groups
        std::vector<NDArray<double>> average, variance, tmp; 
                        //average, error, and temporary data
};

class Tally0D: public Tally
{
    public:
        Tally0D(const std::string& _name);
        void score(const int& m);
        virtual ~Tally0D();
};

class Tally1D: public Tally
{
    public:
        Tally1D();
        Tally1D(std::vector<size_t> _shape, const std::string& _name,
                const size_t _g);
        void score(const int& m);
        ~Tally1D();
};

class Tally2D: public Tally
{
    public:
        Tally2D();
        Tally2D(std::vector<size_t> _shape, const std::string& _name,
                const size_t _g);
        void score(const int& m);
        void score(const NDArray<double> & tmp, int m);
        ~Tally2D();
};

class Tally4D: public Tally
{
    public:
        Tally4D();
        Tally4D(std::vector<size_t> _shape, const std::string& _name,
                const size_t _g);
        void score(const int& m);
        ~Tally4D();

    private:
        std::vector<NDArray<double>>* ptr_to_tmp;
};
#endif
