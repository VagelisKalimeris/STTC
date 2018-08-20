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
// Our data structure
    vector<int> spike_trains[NEURONS];
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
    for (int i = 0; i < NEURONS; i++) {
        for (int j = 0; j < NEURONS; j++) {
            if (i == j) {continue;}
            for (int shift = 0; shift < circ_shifts_num; shift++) {

            }
        }
    }

// Calculate conditional STTC
    int ttl_sgnfcnt_trplts = 0;
    vector<double> to_shift;
    double shifted_res_arr[circ_shifts_num];

    for (int i = 0; i < NEURONS; i++) { // Neuron A
        for (int j = 0; j < NEURONS; j++) { // Neuron B
            if (i == j) {continue;}
            for (int k = 0; k < NEURONS; k++) { // Neuron C
                if (j == k || i == k) {continue;}
                if (!sign_trpl_limit(spike_trains[i], spike_trains[k] ,Dt)) {
                    continue; // Reduced A spike train has < 5 spikes
                }
                double trip_sttc = STTC_AB_C(spike_trains[i], spike_trains[j]
                                                        , spike_trains[k], Dt);
                double shifted_res_arr[circ_shifts_num]; //stores shifted sttc's
                to_shift = spike_trains[k];

                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    circular_shift(to_shift, circ_shifts_num);
                    shifted_res_arr[k] = STTC_AB_C(time_line_A, time_line_B, 
                                                                 to_shift, Dt);
                }
                double mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
                double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
                double threshhold = sign_thresh(mean, st_dev);
                if ( trip_sttc > threshhold) {
                    ttl_sgnfcnt_trplts++;
                }
            }
        }
    } 

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
