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

#include "front_screen_class.hpp"
#include "settings.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cmath>

using namespace std;

FrontScreen::FrontScreen() 
{
    start_time = chrono::high_resolution_clock::now();
}

void FrontScreen::displayName() const 
{
    puts(
"                                            _       \n"
"                                           | |      \n"
"                       ___ ___  _   _  ___ | |_ ___ \n"
"                      / __/ _ \\| | | |/ _ \\| __/ _ \\ \n"
"                     | (_| (_) | |_| | (_) | ||  __/\n"
"                      \\___\\___/ \\__, |\\___/ \\__\\___|\n"
"                                 __/ |              \n"
"                                |___/       \n\n"
"                          Just a coyote in need\n");
}

void FrontScreen::start(int num_threads) const 
{
    time_t now = time(0);
    char* dt = ctime(&now);
    cout << "Date/Time        : " << dt;
    cout << "OpenMP Threads   : " << num_threads << endl; 
    cout << "Input File       : " << settings::input << endl;
    cout << "Output Path      : " << settings::output << endl;
}

void FrontScreen::end(const double& k, const double& k_err, const double& leak, const double& leak_err) 
{
    end_time = chrono::high_resolution_clock::now();
    double t_elapsed = chrono::duration<double>(end_time - start_time).count();

    cout << std::fixed;
    cout << std::setprecision(5);
    cout << " ______________________________________";
    cout << "\n|                                 |" << endl;
    cout << "| kstep = " << k << " +/- " << k_err << " |" << endl;
    cout << "| leak  = " << leak << " +/- " << leak_err << " |" << endl;
    cout << "|______________________________________|" << endl;
    cout << "\nTotal Runtime: " << t_elapsed << " seconds." << endl;
}

void FrontScreen::end(const StaticScores& _scores)
{
    double k = _scores.Kstep().getAverage()[0][0];
    double k_err = sqrt(_scores.Kstep().getVariance()[0][0]/static_cast<double>(settings::m_replica));
    double leak = _scores.Leakage().getAverage()[0][0];
    double leak_err = sqrt(_scores.Leakage().getVariance()[0][0]/static_cast<double>(settings::m_replica));
    double mig_area = _scores.MigArea().getAverage()[0][0];
    double mig_area_err = sqrt(_scores.MigArea().getVariance()[0][0]/static_cast<double>(settings::m_replica));
    end_time = chrono::high_resolution_clock::now();
    double t_elapsed = chrono::duration<double>(end_time - start_time).count();

    cout << std::fixed;
    cout << std::setprecision(5);
    cout << " _____________________________";
    cout << "\n|                             |" << endl;
    cout << "| kstep = " << k << " +/- " << k_err << " |" << endl;
    cout << "| leak  = " << leak << " +/- " << leak_err << " |" << endl;
    cout << "|_____________________________|" << endl;
    cout << "\nMigration Area = " << mig_area << " +/- " << mig_area_err << endl;
    cout << "\nTotal Runtime: " << t_elapsed << " seconds." << endl;
}

void FrontScreen::end() {
    end_time = chrono::high_resolution_clock::now();
    double t_elapsed = chrono::duration<double>(end_time - start_time).count();
    cout << "Total Runtime: " << t_elapsed << " seconds." << endl;
}
