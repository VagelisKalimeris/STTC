/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: main.cpp                                                         *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include "common.hpp"
#include "p_p_null_dist.hpp"
#include "cond_null_dist.hpp"
#include "tuplets_STTC.hpp"
#include "triplets_STTC.hpp"

// neurons 183
// total_time_samples 11970
// circ_shifts_num 50

using namespace std;

/******************************************************************************
* FUNCTION NAME: main                                                         *
*                                                                             *
* ARGUMENTS: The total numbers of: neurons, Time Samples, Circular Shifts.    *
*             And the size of the tile Δt.                                    *
*                                                                             *
* PURPOSE: This main function is for testing purposes only.                   *
*                                                                             *
* RETURNS: 0 on program comletion.                                            *
*                                                                             *
* I/O: Opens/Reads a file containing the spike trains.                        *
*                                                                             *
******************************************************************************/
int main(int argc, char const *argv[])
{
// Prevent warning: unused parameter 'argc'
    (void) argc;
// Command Line Arguments. First give random sample size, then tile size. 
    const int circ_shifts_num = stoi(argv[1]), Dt = stoi(argv[2]);
// Shifted spike trains will be copied here
    vector<int> to_shift;
// STTC values of shifted spike trains
    double shifted_res_arr[circ_shifts_num];

// Caclulation variables
    int ttl_sgnfcnt_tuplets = 0, ttl_sgnfcnt_triplets = 0;
    double tupl_sttc, trip_sttc, mean, st_dev, threshold;

// Open File
    ifstream data;
    data.open("../psm_avalanche", ifstream::in);
    string line;

// Get total number of neurons from file
    getline(data, line);
    // For debugging i choose only 10 neurons
    const int neurons = line.length() - 1;
    data.seekg(0, data.beg);

// Our main data structure
    vector<int> spike_trains[neurons];

// Store each neuron's firing (1's) to the data structure
    int total_time_samples = 0;
    while (getline(data, line)) {
        for (int n = 0; n < neurons; n++) {
            if (line[n] == '1') {
                spike_trains[n].push_back(total_time_samples);
            }
        }
        total_time_samples++;
    }
// Start random sequence
    srand(time(NULL));
    
// Find all T
    double tAp[neurons];
    double tBm[neurons];
    for(int n = 0; n < neurons; ++n) {
        tAp[n] = T_A_plus(spike_trains[n], total_time_samples, Dt);
        tBm[n] = T_B_minus(spike_trains[n], total_time_samples, Dt);
    }

// Calculate per pair STTC
    double tAps;
    for (int a = 0; a < neurons; a++) { // Neuron A
        for (int b = 0; b < neurons; b++) { // Neuron B
            if (a == b) {continue;} // Skip same neurons
            tupl_sttc = STTC_A_B(spike_trains[a], spike_trains[b], 
                                        total_time_samples, Dt, tBm[b], tAp[a]);
            //cout<<"TUPLE STTC: "<<tupl_sttc<<endl;
            for (int shift = 0; shift < circ_shifts_num; shift++) {
                to_shift = spike_trains[a];
                unsigned int random = random_gen(total_time_samples);
                circular_shift(to_shift, random, total_time_samples);
                tAps = T_A_plus(to_shift, total_time_samples, Dt);
                shifted_res_arr[shift] = STTC_A_B(to_shift, spike_trains[b], 
                                        total_time_samples, Dt, tBm[b], tAps);
            }
            mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
            st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
            threshold = sign_thresh(mean, st_dev);
            if (tupl_sttc > threshold) {
                ttl_sgnfcnt_tuplets++;
            }
        }
    }
    cout<<"Number of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<endl; 
    cout<<endl;



// Calculate conditional STTC
    double tApt;
    for (int a = 0; a < neurons; a++) { // Neuron A
        for (int c = 0; c < neurons; c++) { // Neuron C
            if (a == c) {continue;} // Skip same neurons
            if (!sign_trpl_limit(spike_trains[a], spike_trains[c] ,Dt)) {
                continue; // Reduced A spike train has < 5 spikes
            }
            tApt = T_A_plus_tripl(spike_trains[a], spike_trains[c], 
                                                        total_time_samples, Dt);
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (b == a || b == c) {continue;} // Skip same neurons
                trip_sttc = STTC_AB_C(spike_trains[a], spike_trains[b], 
                                    spike_trains[c], total_time_samples, 
                                    Dt, tBm[b], tApt);
                //cout<<"TRIPLE STTC: "<<trip_sttc<<endl;
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    to_shift = spike_trains[c];
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(to_shift, random, total_time_samples);
                    tApt = T_A_plus_tripl(spike_trains[a], to_shift, 
                                                        total_time_samples, Dt);
                    shifted_res_arr[shift] = STTC_AB_C(spike_trains[a], 
                                spike_trains[b], to_shift, total_time_samples, 
                                Dt, tBm[b], tApt);
                }
                mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
                st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
                threshold = sign_thresh(mean, st_dev);
                if ( trip_sttc > threshold) {
                    ttl_sgnfcnt_triplets++;
                }
            }
        }
    }
    cout<<"Number of total significant triplets: "<<ttl_sgnfcnt_triplets<<endl; 


// Print the data structure and total number of firings in experiment
    int total_firings = 0;
    cout<<"\nThe data structure: "<<endl;
    for (int neur = 0; neur < neurons; neur++) {
        cout<<"No "<<neur + 1<<" neuron's spikes:\n";
        for (size_t fire = 0; fire < spike_trains[neur].size(); fire++) {
            cout<<spike_trains[neur][fire] + 1<<' ';
        total_firings++;
        }
        cout<<endl<<endl;
    }
    cout<<"\nTotal number of spikes: "<<total_firings<<endl;
    
    data.close();
    return 0;
}
