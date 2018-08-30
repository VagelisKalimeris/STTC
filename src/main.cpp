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
#include <iomanip>
#include <algorithm>

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
* ARGUMENTS: The total numbers Circular Shifts and the size of the tile Δt.   *
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
    const int circ_shifts_num = atoi(argv[1]), Dt = atoi(argv[2]);
    
// Caclulation variables
    int ttl_sgnfcnt_tuplets = 0, ttl_sgnfcnt_triplets = 0;
    
// Open File
    ifstream data, extra;
    data.open("../psm_avalanche", ifstream::in);
    extra.open("../astrocytes", ifstream::in);
    string line;
    
// Our main astrocytes structure
    vector<int> astrocytes;
    
// Get astrocytes
    while (getline(extra, line)) {
        astrocytes.push_back(stoi(line) - 1);
    }
    int astrocytes_size = astrocytes.size();
    
// Get total number of neurons from file
    getline(data, line);
    const int neurons = line.length() - 1;
    data.seekg(0, data.beg);
    const int neur_clean = neurons - astrocytes_size;
    
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
    int astro = 0;
    int astrocyte = astrocytes[astro];
    for(int n = 0; n < neurons; ++n) {
        if (astrocyte == n) {
            astrocyte = astrocytes[(++astro) % astrocytes_size];
            continue;
        }
        tAp[n] = T_A_plus(spike_trains[n], total_time_samples, Dt);
        tBm[n] = T_B_minus(spike_trains[n], total_time_samples, Dt);
    }
    
// Significant tuplets
    bool sgnfcnt_tuplets[neurons][neurons];
    
// Calculate per pair STTC
    // print_sgnfcnt_tuplet_begin();
    ofstream tuplets;
    tuplets.open("tuplets.csv");
    tuplets<<"NeuronA,NeuronB,STTC,Percentile"<<endl;
    int astro_a = 0;
    int astrocyte_a = astrocytes[astro_a];
    for (int a = 0; a < neurons; a++) { // Neuron A
        if (astrocyte_a == a) {
            astrocyte_a = astrocytes[(++astro_a) % astrocytes_size];
            continue;
        }
        vector<int> time_line_A = spike_trains[a];
        double tAp_tmp = tAp[a];
        int astro_b = 0;
        int astrocyte_b = astrocytes[astro_b];
        for (int b = 0; b < neurons; b++) { // Neuron B
            if (astrocyte_b == b) {
                astrocyte_b = astrocytes[(++astro_b) % astrocytes_size];
                continue;
            }
            sgnfcnt_tuplets[a][b] = false;
            if (a == b) {continue;} // Skip same neurons
            vector<int> time_line_B = spike_trains[b];
            double tBm_tmp = tBm[b];
            double tupl_sttc = STTC_A_B(time_line_A, time_line_B, 
                                                        Dt, tBm_tmp, tAp_tmp);
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            for (int shift = 0; shift < circ_shifts_num; shift++) {
            // Shifted spike trains will be copied here
                vector<int> to_shift = time_line_A;
                unsigned int random = random_gen(total_time_samples);
                circular_shift(to_shift, random, total_time_samples);
                tAp_tmp = T_A_plus(to_shift, total_time_samples, Dt);
                shifted_res_arr[shift] = STTC_A_B(to_shift, time_line_B, 
                                                        Dt, tBm_tmp, tAp_tmp);
            }
            double mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
            double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
            double threshold = sign_thresh(mean, st_dev);
            if (tupl_sttc > threshold) {
                ttl_sgnfcnt_tuplets++;
                sgnfcnt_tuplets[a][b] = true;
                sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                int pos = 0; 
                while (pos < circ_shifts_num && 
                                        shifted_res_arr[pos] <= tupl_sttc) {
                    ++pos;
                }
                // print_sgnfcnt_tuplet(a+1, b+1, tupl_sttc, 
                //                             pos/double(circ_shifts_num));
                tuplets<<a+1<<','<<b+1<<','<<tupl_sttc<<','
                                        <<pos/double(circ_shifts_num)<<endl;
            }
        }
    }
    // print_sgnfcnt_tuplet_end();
    tuplets.close();
    cout<<"\nNumber of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<" ( "
                            <<(ttl_sgnfcnt_tuplets * 100 / double(neur_clean * 
                            (neur_clean - 1)))<<"% )"<<endl;
    
    
    
// Motif arrays
    int motifs_triplets[8] = {0};
    int motifs_sgnfcnts[8] = {0};
    
// Calculate conditional STTC
    // print_sgnfcnt_triplet_begin();
    ofstream triplets;
    triplets.open("triplets.csv");
    triplets<<"NeuronA,NeuronB,NeuronC,STTC,Percentile"<<endl;
    astro_a = 0;
    astrocyte_a = astrocytes[astro_a];
    for (int a = 0; a < neurons; a++) { // Neuron A
        if (astrocyte_a == a) {
            astrocyte_a = astrocytes[(++astro_a) % astrocytes_size];
            continue;
        }
        vector<int> time_line_A = spike_trains[a];
        int astro_c = 0;
        int astrocyte_c = astrocytes[astro_c];
        for (int c = 0; c < neurons; c++) { // Neuron C
            if (astrocyte_c == c) {
                astrocyte_c = astrocytes[(++astro_c) % astrocytes_size];
                continue;
            }
            if (a == c) {continue;} // Skip same neurons
            vector<int> time_line_C = spike_trains[c];
            int astro_b = 0;
            int astrocyte_b = astrocytes[astro_b];
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (astrocyte_b == b) {
                    astrocyte_b = astrocytes[(++astro_b) % astrocytes_size];
                    continue;
                }
                if (b == a || b == c) {continue;} // Skip same neurons
                int pos = sgnfcnt_tuplets[c][a] * 4 + 
                        sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b] * 1;
                ++motifs_triplets[pos];
            }
            if (!sign_trpl_limit(time_line_A, time_line_C ,Dt)) {
                continue; // Reduced A spike train has < 5 spikes
            }
            double tApt = T_A_plus_tripl(time_line_A, time_line_C, 
                                                    total_time_samples, Dt);
            astro_b = 0;
            astrocyte_b = astrocytes[astro_b];
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (astrocyte_b == b) {
                    astrocyte_b = astrocytes[(++astro_b) % astrocytes_size];
                    continue;
                }
                if (b == a || b == c) {continue;} // Skip same neurons
                vector<int> time_line_B = spike_trains[b];
                double tBm_tmp = tBm[b];
                double trip_sttc = STTC_AB_C(time_line_A, time_line_B, 
                                            time_line_C, Dt, tBm_tmp, tApt);
            // STTC values of shifted spike trains
                double shifted_res_arr[circ_shifts_num];
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                // Shifted spike trains will be copied here
                    vector<int> to_shift = time_line_C;
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(to_shift, random, total_time_samples);
                    tApt = T_A_plus_tripl(time_line_A, to_shift, 
                                                    total_time_samples, Dt);
                    shifted_res_arr[shift] = STTC_AB_C(time_line_A, 
                                    time_line_B, to_shift, Dt, tBm_tmp, tApt);
                }
                double mean = mean_STTC_dir(shifted_res_arr, circ_shifts_num);
                double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num);
                double threshold = sign_thresh(mean, st_dev);
                if ( trip_sttc > threshold) {
                    ttl_sgnfcnt_triplets++;
                    int pos = sgnfcnt_tuplets[c][a] * 4 + 
                            sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b];
                    ++motifs_sgnfcnts[pos];
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    pos = 0; 
                    while (pos < circ_shifts_num && 
                                        shifted_res_arr[pos] <= trip_sttc) {
                        ++pos;
                    }
                    // print_sgnfcnt_triplet(a+1, b+1, c+1, trip_sttc, 
                    //                         pos/double(circ_shifts_num));
                    triplets<<a+1<<','<<b+1<<','<<c+1<<','<<trip_sttc<<','
                                        <<pos/double(circ_shifts_num)<<endl;
                }
            }
        }
    }
    // print_sgnfcnt_triplet_end();
    triplets.close();
    cout<<"\nNumber of total significant triplets: "<<ttl_sgnfcnt_triplets
                    <<" ( "<<(ttl_sgnfcnt_triplets * 100 / double(neur_clean * 
                    (neur_clean - 1) * (neur_clean - 2)))<<"% )"<<endl;
    
    
// Print Motifs
    print_motifs(motifs_triplets, motifs_sgnfcnts);
    
// Print the data structure and total number of firings in experiment
    // print_all_spikes(spike_trains, neurons, astrocytes);
    
    extra.close();
    data.close();
    return 0;
}
