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
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <omp.h>
#include <memory>
#include <unistd.h>
#include <algorithm>
#include <filesystem>
#include "pcg_random.hpp"
#include "simulation_class.hpp"
#include "scores_class.hpp"
#include "ressources.hpp"
#include "front_screen_class.hpp"
#include "error.hpp"
#include "settings.hpp"
#include "parser.hpp"
#include "static_simulation_class.hpp"
#include "static_scores_class.hpp"

using namespace std;

const char* OPTIONS = "hr:i:o:t:s:";

void usage() 
{
    cout << "\n   -i <input folder>   : <input folder> path to input folder\n";
    cout << "   -o <output folder>  : <output folder> path to output folder\n";
    cout << "   -r <replicas>       : <replicas> number of tarapiaTapioco replicas, default to 0\n";
    cout << "   -t <num threads>    : <num_threads> number of omp threads, default to 48\n";
    cout << "   -s                  : Random seed for antani\n";
    cout << "   -h                  : Show help and exit\n\n";
}

int main(int argc, char* argv[]) // arg is M
{
    FrontScreen front_screen;
    int num_threads = omp_get_max_threads();

    int opt = 0;
    while ((opt = getopt(argc, argv, OPTIONS)) != -1)
        switch (opt)
        {
            case 'r':
                settings::m_replica = stoi(optarg);
                break;
            case 't':
                num_threads = stoi(optarg);
                omp_set_num_threads(num_threads);
                break;
            case 'i':
                settings::input = string(optarg);
                break;
            case 'o':
                settings::output = string(optarg);
                break;
            case 's': 
                settings::pi_seed = stoi(optarg);
                break;
            case 'h':
                usage();
                return 0;
                break;
            default:
                usage();
                return 1;
                break;
        }
    const int period = max(1, settings::m_replica/100);
    front_screen.displayName();
    front_screen.start(num_threads);

    parseInput();
    Scores global_scores;

    int PI_ended = false;

    vector<double> Tp(settings::t_bin+1,0);
    for (int i = 0; i < settings::t_bin+1; i++)
        Tp[i] = i*settings::dt;

    if (settings::simulation_mode == settings::SimulationMode::TARAPIATAPIOCO) {
        cout << "Generating source particles." << endl;
        #ifdef _OPENMP
            #pragma omp parallel
        #endif
        {
            unique_ptr<Simulation> sim_local = make_unique<Simulation>();

            if (not PI_ended)
            {
                #pragma omp atomic
                    PI_ended++;
                cout << "Starting tarapiaTapioco simulation." << endl;
            }

            #ifdef _OPENMP
                #pragma omp for schedule(dynamic)
            #endif
            for (int i = 1; i <= settings::m_replica; i++)
            {
               if (i%period == 0)
                {
                    printf("Computation at %d percents.\r", i/period);
                    fflush(stdout);
                }
               
                sim_local->reinitialize();
                sim_local->setSource();
                for (int t = 0; t < settings::t_bin; t++)
                {
                    if (count(settings::reactivity_time.begin(), settings::reactivity_time.end(), t) and settings::reactivity_change)
                        sim_local->changeReactivity();
                    sim_local->advanceTime(t+1);
                        //advancing time to t + dt

                    if (settings::control_strat == settings::ControlStrat::REGULAR) {
                        sim_local->applyPopulationControl();
                        sim_local->cleanBuffer();
                    }
                    else if (settings::control_strat == settings::ControlStrat::SINGULAR) {
                        if (count(settings::control_steps.begin(), 
                                  settings::control_steps.end(), t)) {
                            sim_local->applyPopulationControl();
                            sim_local->cleanBuffer();
                        }
                    }


                    if (settings::entropy)
                        sim_local->scoreEntropy(t+1);
                    if (settings::pair_distance)
                        sim_local->scorePairDist(t+1);
                } //end of time loop
                sim_local->computeScores();
                    //on the fly scores with welford formulaes
            }
            sim_local->normalize();
            #ifdef _OPENMP
                #pragma omp critical
            #endif
            {
                global_scores += sim_local->getScores();
            }
        }//end of pragma parallel
        cout << "\nSaving scores." << endl;
        global_scores.saveScores(settings::output);
        front_screen.end();
    }

    else {
        cout << "Now using new static simulation capabilities" << endl;

        StaticScores pi_scores;
        #pragma omp parallel
        {
            StaticSimulation pi;
            #pragma omp for schedule(dynamic)
            for (int i = 1; i <= settings::m_replica; i++)
            {
                if (i%period == 0)
                {
                    printf("Computation at %d percents.\r", i/period);
                    fflush(stdout);
                }
                pi.reinitialize();
                pi.source();
                for (int j = 1; j <= settings::passive_gen + settings::active_gen; j++) {
                    pi.antani(j);
                    pi.update(j);
                }
                pi.makeScores();
                if (i == settings::m_replica)
                    pi.takePhoto(i, settings::output); 
            }
            pi.computeFeynmanMomentRep();
            pi.normalize();
            #pragma omp critical
            {
                pi_scores += pi.Scores();
            }
        }
        cout << "\nSaving scores..." << endl;
        pi_scores.save();
        front_screen.end(pi_scores);
    }
    return 0;
}
