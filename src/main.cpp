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

#define NEURONS 183
// total_time_samples 11970
// circ_shifts_num 50

using namespace std;

/******************************************************************************
* FUNCTION NAME: main                                                         *
*                                                                             *
* ARGUMENTS: The total numbers of: Neurons, Time Samples, Circular Shifts.    *
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
// command line arguments
    const int total_time_samples = stoi(argv[2]),
                           circ_shifts_num = stoi(argv[3]), Dt = stoi(argv[3]);
// Shifted spike trains will be copied here
    vector<int> to_shift;
// STTC values of shifted spike trains
    double shifted_res_arr[circ_shifts_num];
// Our main data structure
    vector<int> spike_trains[NEURONS];
// Caclulation variables
    int ttl_sgnfcnt_tuplets = 0, ttl_sgnfcnt_triplets = 0;
    double tupl_sttc, trip_sttc, mean, st_dev, threshold;

// Open File
    ifstream data;
    data.open("../psm_avalanche", ifstream::in);
    string line;

// Store each neuron's firing (1's) to the data structure
    int count = 0;
    while (getline(data, line)) {
        for (int n = 0; n < NEURONS; n++) {
            if (line[n] == '1') {
                spike_trains[n].push_back(count);
            }
            count++;
        }
    }


// Calculate per pair STTC
    for (int i = 0; i < NEURONS; i++) { // Neuron A
        for (int j = 0; j < NEURONS; j++) { // Neuron B
            if (i == j) {continue;} // Skip same neurons
            tupl_sttc = STTC_A_B(spike_trains[i], spike_trains[j], Dt);

            for (int shift = 0; shift < circ_shifts_num; shift++) {
                to_shift = spike_trains[i];
                circular_shift(to_shift, circ_shifts_num);
                shifted_res_arr[shift] = STTC_A_B(spike_trains[i], 
                                                          spike_trains[j], Dt);
            }
            mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
            st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
            threshold = sign_thresh(mean, st_dev);
            if ( tupl_sttc > threshold) {
                    ttl_sgnfcnt_tuplets++;
                }
        }
    }
    cout<<"Number of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<endl; 


// Calculate conditional STTC
    for (int i = 0; i < NEURONS; i++) { // Neuron A
        for (int j = 0; j < NEURONS; j++) { // Neuron B
            if (i == j) {continue;} // Skip same neurons
            for (int k = 0; k < NEURONS; k++) { // Neuron C
                if (j == k || i == k) {continue;} // Skip same neurons
                if (!sign_trpl_limit(spike_trains[i], spike_trains[k] ,Dt)) {
                    continue; // Reduced A spike train has < 5 spikes
                }
                trip_sttc = STTC_AB_C(spike_trains[i], spike_trains[j]
                                                        , spike_trains[k], Dt);

                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    to_shift = spike_trains[k];
                    circular_shift(to_shift, circ_shifts_num);
                    shifted_res_arr[shift] = STTC_AB_C(spike_trains[i], 
                                                spike_trains[j], to_shift, Dt);
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
    for (int neur = 0; neur < NEURONS; neur++) {
        for (int fire = 0; fire < spike_trains[neur].size(); fire++) {
            cout<<spike_trains[neur][fire]<<' '<<endl;
        total_firings++;
        }
        cout<<endl;
    }
    cout<<endl<<total_firings<<endl;
    
    data.close();
    return 0;
}
