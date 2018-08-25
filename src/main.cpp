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
#include <ctime>

#include "common.hpp"
#include "p_p_null_dist.hpp"
#include "cond_null_dist.hpp"
#include "tuplets_STTC.hpp"
#include "triplets_STTC.hpp"
#include "motif.hpp"

// neurons 183
// total_time_samples 11970
// circ_shifts_num 50

using namespace std;

/******************************************************************************
* FUNCTION NAME: main                                                         *
*                                                                             *
* ARGUMENTS: The total numbers Circular Shifts and the size of the tile Δt.   *                                    *
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
    auto main_func = [](){
// Prevent warning: unused parameter 'argc'
    (void) argc;
// Command Line Arguments. First give random sample size, then tile size. 
    const int circ_shifts_num = atoi(argv[1]), Dt = atoi(argv[2]);
// Shifted spike trains will be copied here
    vector<int> to_shift;
// STTC values of shifted spike trains
    //double shifted_res_arr[circ_shifts_num];

// Caclulation variables
    int ttl_sgnfcnt_tuplets = 0, ttl_sgnfcnt_triplets = 0;
    //double tupl_sttc, trip_sttc, mean, st_dev, threshold;

// Open File
    ifstream data;
    data.open("C:\\Users\\MANOS\\Source\\Repos\\STTC\\psm_avalanche", ifstream::in);
    string line;

// Get total number of neurons from file
    getline(data, line);
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
    data.close();
// Start random sequence
    srand((unsigned int)time(NULL));
    
// Find all T
    double tAp[neurons];
    double tBm[neurons];
    #pragma omp parallel for
    for(int n = 0; n < neurons; ++n) {
        tAp[n] = T_A_plus(spike_trains[n], total_time_samples, Dt);
        tBm[n] = T_B_minus(spike_trains[n], total_time_samples, Dt);
    }

// Significant tuplets
    bool sgnfcnt_tuplets[neurons][neurons] = {false};

// Calculate per pair STTC
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        double tAp_tmp = tAp[a];
        #pragma omp parallel for private(tAp_tmp)
        for (int b = 0; b < neurons; b++) { // Neuron B
            if (a == b) {continue;} // Skip same neurons
            vector<int> time_line_B = spike_trains[b];
            double tBm_tmp = tBm[b];
            double tupl_sttc = STTC_A_B(time_line_A, time_line_B, 
                                        total_time_samples, Dt, tBm_tmp, tAp_tmp);
            double shifted_res_arr[circ_shifts_num];
            for (int shift = 0; shift < circ_shifts_num; shift++) {
                vector<int> to_shift = time_line_A;
                unsigned int random = random_gen(total_time_samples);
                circular_shift(to_shift, random, total_time_samples);
                tAp_tmp = T_A_plus(to_shift, total_time_samples, Dt);
                shifted_res_arr[shift] = STTC_A_B(to_shift, time_line_B, 
                                        total_time_samples, Dt, tBm_tmp, tAp_tmp);
            }
            double mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
            double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
            double threshold = sign_thresh(mean, st_dev);
            if (tupl_sttc > threshold) {
                #pragma omp atomic
                ttl_sgnfcnt_tuplets++;
                sgnfcnt_tuplets[a][b] = true;
            }
        }
    }
    cout<<"\nNumber of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<endl;



// Motif arrays
    unsigned int motifs_triplets[8] = {0};
    unsigned int motifs_sgnfcnt[8] = {0};
    
// Calculate conditional STTC
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        #pragma omp parallel for
        for (int c = 0; c < neurons; c++) { // Neuron C
            if (a == c) {continue;} // Skip same neurons
            vector<int> time_line_C = spike_trains[c];
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (b == a || b == c) {continue;} // Skip same neurons
                categorization(sgnfcnt_tuplets[c][a], sgnfcnt_tuplets[c][b], 
                                        sgnfcnt_tuplets[a][b], motifs_triplets);
            }
            if (!sign_trpl_limit(time_line_A, time_line_C ,Dt)) {
                continue; // Reduced A spike train has < 5 spikes
            }
            double tApt = T_A_plus_tripl(time_line_A, time_line_C, 
                                                        total_time_samples, Dt);
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (b == a || b == c) {continue;} // Skip same neurons
                vector<int> time_line_B = spike_trains[b];
                double tBm_tmp = tBm[b];
                double trip_sttc = STTC_AB_C(time_line_A, time_line_B, 
                                    time_line_C, total_time_samples, 
                                    Dt, tBm_tmp, tApt);
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    vector<int> to_shift = time_line_C;
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(to_shift, random, total_time_samples);
                    tApt = T_A_plus_tripl(time_line_A, to_shift, 
                                                        total_time_samples, Dt);
                    shifted_res_arr[shift] = STTC_AB_C(time_line_A, 
                                time_line_B, to_shift, total_time_samples, 
                                Dt, tBm_tmp, tApt);
                }
                double mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
                double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
                double threshold = sign_thresh(mean, st_dev);
                if ( trip_sttc > threshold) {
                    #pragma omp atomic
                    ttl_sgnfcnt_triplets++;
                    categorization(sgnfcnt_tuplets[c][a], sgnfcnt_tuplets[c][b], 
                                         sgnfcnt_tuplets[a][b], motifs_sgnfcnt);
                }
            }
        }
    }
    cout<<"Number of total significant triplets: "<<ttl_sgnfcnt_triplets<<endl<<endl; 


// Print Motifs
    /*cout<<"\nMotif - Triplets - Significant"<<endl;
    for (int m = 0; m < 8; m++) {
        cout<<"  "<<m<<"   - "<<motifs_triplets[m]<<(
                                int(motifs_triplets[m]/10000000)?"":
                                int(motifs_triplets[m]/1000000)?" ":
                                int(motifs_triplets[m]/100000)?"  ":
                                int(motifs_triplets[m]/10000)?"   ":
                                int(motifs_triplets[m]/1000)?"    ":
                                int(motifs_triplets[m]/100)?"     ":
                                int(motifs_triplets[m]/10)?"      ":
                                "        ")<<" - "<<motifs_sgnfcnt[m]<<endl;
    }*/
    
// Print the data structure and total number of firings in experiment
    // print_all_spikes(spike_trains, neurons);
    }

    
    
    return 0;
}
