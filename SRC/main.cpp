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

#include "../INCLUDE/common.hpp"
#include "../INCLUDE/p_p_null_dist.hpp"
#include "../INCLUDE/cond_null_dist.hpp"
#include "../INCLUDE/tuplets_STTC.hpp"
#include "../INCLUDE/triplets_STTC.hpp"
#include "../INCLUDE/motif.hpp"

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
    
// Open Files
    ifstream data, astros;
    data.open((string("DATASETS/") + argv[3]).c_str(), ifstream::in);
    astros.open((string("ASTROCYTES/") + argv[3]).c_str(), ifstream::in);
    if (!data.is_open() || !astros.is_open()) {
        cout<<"Error opening input file!"<<endl;
        return 0;
    }
    string line;
    
// Get astrocytes
    vector<int> astrocytes;
    while (getline(astros, line)) {
        astrocytes.push_back(atoi(line.c_str()) - 1);
    }
    const int astro_size = astrocytes.size();

// Get total number of neurons from file
    getline(data, line);
    const int neurons = line.length() - astro_size - 1;
    data.seekg(0, data.beg);
    
// Our main data structure and astrocyte list
    vector<int> spike_trains[neurons];
    
// Store each neuron's firing (1's) to the data structure
    int total_time_samples = 0;
    int astros_gone = 0;
    while (getline(data, line)) {
        for (int n = 0; n < neurons; n++) {
            if ((n+astros_gone) == astrocytes[astros_gone]) {astros_gone++; continue;}
            if (line[n+astros_gone] == '1') {
                spike_trains[n].push_back(total_time_samples);
            }
        }
        total_time_samples++;
    }

// Close files
    data.close();
    astros.close();
    
// Start random sequence
    srand(time(NULL));
    
// Find all T
    double tAp[neurons];
    double tBm[neurons];
    #pragma omp parallel for
    for(int n = 0; n < neurons; ++n) {
        tAp[n] = T_A_plus(spike_trains[n], total_time_samples, Dt);
        tBm[n] = T_B_minus(spike_trains[n], total_time_samples, Dt);
    }
    
// Significant tuplets
    bool sgnfcnt_tuplets[neurons][neurons];
    
// Calculate per pair STTC
    ofstream tuplets;
    tuplets.open(("RESULTS/" + string(argv[3]) + "_tuplets.csv").c_str());
    if (!tuplets.is_open()) {
        cout<<"Error opening results file!"<<endl;
        return 0;
    }
    tuplets<<"NeuronA,NeuronB,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        #pragma omp parallel
        {
            double tAp_tmp = tAp[a];
            #pragma omp for
            for (int b = 0; b < neurons; b++) { // Neuron B
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
                    #pragma omp atomic
                    ++ttl_sgnfcnt_tuplets;
                    sgnfcnt_tuplets[a][b] = true;
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    int pos = 0; 
                    while (pos < circ_shifts_num && 
                                            shifted_res_arr[pos] <= tupl_sttc) {
                        ++pos;
                    }
                    #pragma omp critical
                    tuplets<<a+1<<','<<b+1<<','<<tupl_sttc<<','
                                            <<pos/double(circ_shifts_num)<<'\n';
                }
            }
        }
    }
    // print_sgnfcnt_tuplet_end();
    tuplets.close();
    cout<<"\nNumber of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<" ( "
        <<(ttl_sgnfcnt_tuplets*100/double(neurons*(neurons-1)))<<"% )"<<endl;
    
    
    
// Motif arrays
    int motifs_triplets[8] = {0};
    int motifs_sgnfcnts[8] = {0};
    
// Calculate conditional STTC
    ofstream triplets;
    triplets.open(("RESULTS/" + string(argv[3]) + "_triplets.csv").c_str());
    if (!triplets.is_open()) {
        cout<<"Error opening results file!"<<endl;
        return 0;
    }
    triplets<<"NeuronA,NeuronB,NeuronC,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        #pragma omp parallel for
        for (int c = 0; c < neurons; c++) { // Neuron C
            if (a == c) {continue;} // Skip same neurons
            vector<int> time_line_C = spike_trains[c];
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (b == a || b == c) {continue;} // Skip same neurons
                int pos = sgnfcnt_tuplets[c][a] * 4 + 
                        sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b] * 1;
                #pragma omp atomic
                ++motifs_triplets[pos];
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
                    #pragma omp atomic
                    ++ttl_sgnfcnt_triplets;
                    int pos = sgnfcnt_tuplets[c][a] * 4 + 
                            sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b] * 1;
                    #pragma omp atomic
                    ++motifs_sgnfcnts[pos];
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    pos = 0; 
                    while (pos < circ_shifts_num && 
                                            shifted_res_arr[pos] <= trip_sttc) {
                        ++pos;
                    }
                    #pragma omp critical
                    triplets<<a+1<<','<<b+1<<','<<c+1<<','<<trip_sttc<<','
                                        <<pos/double(circ_shifts_num)<<'\n';
                }
            }
        }
    }
    triplets.close();
    cout<<"\nNumber of total significant triplets: "<<ttl_sgnfcnt_triplets<<" ( "
            <<(ttl_sgnfcnt_triplets*100/double(neurons*(neurons-1)*(neurons-2)))
            <<"% )"<<endl;
    
    
// Print Motifs
    print_motifs(motifs_triplets, motifs_sgnfcnts);
    
// Print the data structure and total number of firings in experiment
    // print_all_spikes(spike_trains, neurons);
    
    return 0;
}
