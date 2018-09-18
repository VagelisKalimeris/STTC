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
#include <cstdio>

#include "../INCLUDE/common.hpp"
#include "../INCLUDE/p_p_null_dist.hpp"
#include "../INCLUDE/cond_null_dist.hpp"
#include "../INCLUDE/tuplets_STTC.hpp"
#include "../INCLUDE/triplets_STTC.hpp"
#include "../INCLUDE/print.hpp"

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
// Command line arguments check
    if (argc != 4) {
        cout<<"Error! Wrong parameter count!"<<endl;
        return 0;
    }
// Command Line Arguments. First give random sample size, then tile size. 
    const int circ_shifts_num = atoi(argv[1]), Dt = atoi(argv[2]);
    
// Caclulation variables
    int ttl_sgnfcnt_tuplets = 0, ttl_sgnfcnt_triplets = 0;
    
// Open File
    ifstream data, astros;
    data.open((string("DATASETS/") + argv[3]).c_str(), ifstream::in);
    if (!data.is_open()) {
        cout<<"Error opening dataset file!"<<endl;
        return 0;
    }
    astros.open((string("ASTROCYTES/") + argv[3]).c_str(), ifstream::in);
    if (!astros.is_open()) {
        cout<<"Error opening astrocytes file!"<<endl;
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
    vector<int> spike_trains[neurons + astro_size];
    
// Store each neuron's firing (1's) to the data structure
    int total_time_samples = 0;
    while (getline(data, line)) {
        int astros_count = 0;
        int astrocyte = astrocytes[0];
        int push_count = 0;
        for (int neur = 0; neur < neurons + astro_size; ++neur) {
            int pos;
            if (neur == astrocyte) {
                pos = neurons + astros_count++;
                astrocyte = astrocytes[astros_count];
            }
            else {
                pos = push_count++;
            }
            if (line[neur] == '1') {
                spike_trains[pos].push_back(total_time_samples);
            }
        }
        total_time_samples++;
    }
    
// Make the mapping from virtual to real neuron's number
    int map[neurons];
    int astro = 0;
    int astrocyte = astrocytes[0];
    for (int neur = 0; neur < neurons; ++neur) {
        while ((neur + astro) == astrocyte) {
            astrocyte = astrocytes[(++astro) % astro_size];
        }
        map[neur] = neur + astro;
    }

// Close input files
    data.close();
    astros.close();
    
// Print the data structure and total number of firings in experiment
    char str[33];
    sprintf(str, "%d", circ_shifts_num);
    const string shifts_s = string(str);
    sprintf(str, "%d", Dt);
    const string Dt_s = string(str);
    ofstream info;
    info.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + "-shifts_" + 
                                    Dt_s + "-dt_neurons_info.txt").c_str());
    if (!info.is_open()) {
        cout<<"Error opening results neurons info file!"<<endl;
        return 0;
    }
    print_all_spikes(spike_trains, neurons + astro_size, astrocytes, info, 
                                            string(argv[3]), shifts_s, Dt_s);
    
// Start random sequence
    srand(time(NULL));
    
// Time lines' arrays
    int tl_sizes[neurons];
    int* tl_array[neurons];
    int tl_size_max = 0;
    vector<int> spike_train;
    for (int neur = 0; neur < neurons; ++neur) {
        spike_train = spike_trains[neur];
        int tl_size = spike_train.size();
        tl_sizes[neur] = tl_size;
        tl_array[neur] = (int *)malloc(tl_size * sizeof(int));
        for (int ts = 0; ts < tl_size; ++ts) {
            tl_array[neur][ts] = spike_train[ts];
        }
        if (tl_size_max < tl_size) {
            tl_size_max = tl_size;
        }
    }
    
// All T for tuplets
    double T_Aplus[neurons];
    double T_Bminus[neurons];
    #pragma omp parallel for
    for(int neur = 0; neur < neurons; ++neur) {
        int* tl = tl_array[neur];
        int tl_size = tl_sizes[neur];
        T_Aplus[neur] = T_A_plus(tl, tl_size, total_time_samples, Dt);
        T_Bminus[neur] = T_B_minus(tl, tl_size, total_time_samples, Dt);
    }
    
// Significant tuplets
    bool sgnfcnt_tuplets[neurons][neurons];
    
// Significant limit
    bool sgnfcnt_limit[neurons][neurons];
    
// Reduced spiketrain T for triplets
    double T_Aplus_tripl[neurons][neurons];
    
// Calculate per pair STTC
    ofstream tuplets;
    tuplets.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + "-shifts_" + 
                                            Dt_s + "-dt_tuplets.csv").c_str());
    if (!tuplets.is_open()) {
        cout<<"Error opening results tuplets file!"<<endl;
        return 0;
    }
    tuplets<<"NeuronA,NeuronB,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        int* tl_A = tl_array[a];
        int tl_A_size = tl_sizes[a];
        if (tl_A_size == 0) {continue;}
        int a_real = map[a];
        #pragma omp parallel
        {
            double tAp = T_Aplus[a];
        // Shifted spike trains will be copied here
            int* to_shift = (int *)malloc(tl_size_max * sizeof(int));
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            #pragma omp for
            for (int b = 0; b < neurons; b++) { // Neuron B
            // It will be used to help in categorization of motifs
                sgnfcnt_tuplets[a][b] = false;
                if (a == b) {continue;} // Skip same neurons
                int* tl_B = tl_array[b];
                int tl_B_size = tl_sizes[b];
                if (tl_B_size == 0) {continue;}
                sgnfcnt_limit[a][b] = sign_trpl_limit(tl_A, tl_A_size, tl_B, 
                                                                tl_B_size, Dt);
                if (sgnfcnt_limit[a][b]) {
                    T_Aplus_tripl[a][b] = T_A_plus_tripl(tl_A, tl_A_size, 
                                    tl_B, tl_B_size, total_time_samples, Dt);
                }
                double tBm = T_Bminus[b];
                double tupl_sttc = STTC_A_B(tl_A, tl_A_size, tl_B, tl_B_size, 
                                                                Dt, tBm, tAp);
                if (tupl_sttc == 2.0) {continue;}
                int denominator = circ_shifts_num;
                double mean = 0;
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(to_shift, tl_A, tl_A_size, random, 
                                                        total_time_samples);
                    tAp = T_A_plus(to_shift, tl_A_size, total_time_samples, 
                                                                        Dt);
                    shifted_res_arr[shift] = STTC_A_B(to_shift, tl_A_size, 
                                                tl_B, tl_B_size, Dt, tBm, tAp);
                    if (shifted_res_arr[shift] == 2.0) {
                        --denominator;
                    }
                    else {
                        mean += shifted_res_arr[shift];
                    }
                }
                mean /= denominator;
                double st_dev = std_STTC_dir(shifted_res_arr, circ_shifts_num, 
                                                                        mean);
                double threshold = sign_thresh(mean, st_dev);
                if (tupl_sttc > threshold) {
                    #pragma omp atomic
                    ++ttl_sgnfcnt_tuplets;
                    sgnfcnt_tuplets[a][b] = true;
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    int pos = 0; 
                    while (pos < denominator && 
                                        shifted_res_arr[pos] <= tupl_sttc) {
                        ++pos;
                    }
                    int b_real = map[b];
                    double percentile = pos / double(denominator);
                    #pragma omp critical
                    tuplets<<a_real + 1<<','<<b_real + 1<<','<<tupl_sttc<<','
                                                            <<percentile<<'\n';
                }
            }
            free(to_shift);
        }
    }
    tuplets.close();
    info<<"\nNumber of total significant tuplets: "<<ttl_sgnfcnt_tuplets<<" ( "
                            <<(ttl_sgnfcnt_tuplets * 100 / double(neurons * 
                            (neurons - 1)))<<"% )"<<endl;
    
    
// Motif arrays
    int motifs_triplets[8] = {0};
    int motifs_sgnfcnts[8] = {0};
    
// Calculate conditional STTC
    ofstream triplets;
    triplets.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + 
                            "-shifts_" + Dt_s + "-dt_triplets.csv").c_str());
    if (!triplets.is_open()) {
        cout<<"Error opening results triplets file!"<<endl;
        return 0;
    }
    triplets<<"NeuronA,NeuronB,NeuronC,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        int* tl_A = tl_array[a];
        int tl_A_size = tl_sizes[a];
        if (tl_A_size == 0) {continue;}
        int a_real = map[a];
        #pragma omp parallel
        {
        // Shifted spike trains will be copied here
            int* to_shift = (int *)malloc(tl_size_max * sizeof(int));
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            #pragma omp for
            for (int c = 0; c < neurons; c++) { // Neuron C
                if (a == c) {continue;} // Skip same neurons
                int* tl_C = tl_array[c];
                int tl_C_size = tl_sizes[c];
                if (tl_C_size == 0) {continue;}
                bool sign_trplt_limit = sgnfcnt_limit[a][c];
                double tApt = T_Aplus_tripl[a][c];
                int c_real = map[c];
                for (int b = 0; b < neurons; b++) { // Neuron B
                    if (b == a || b == c) {continue;} // Skip same neurons
                    int pos = sgnfcnt_tuplets[c][a] * 4 + 
                            sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b];
                    #pragma omp atomic
                    ++motifs_triplets[pos];
                    if (!sign_trplt_limit) {
                        continue; // Reduced A spike train has < 5 spikes
                    }
                    int* tl_B = tl_array[b];
                    int tl_B_size = tl_sizes[b];
                    if (tl_B_size == 0) {continue;}
                    double tBm = T_Bminus[b];
                    double trip_sttc = STTC_AB_C(tl_A, tl_A_size, tl_B, 
                                    tl_B_size, tl_C, tl_C_size, Dt, tBm, tApt);
                    if (trip_sttc == 2.0) {
                        #pragma omp atomic
                        --motifs_triplets[pos];
                        continue;
                    }
                    int denominator = circ_shifts_num;
                    double mean = 0;
                    for (int shift = 0; shift < circ_shifts_num; shift++) {
                        unsigned int random = random_gen(total_time_samples);
                        circular_shift(to_shift, tl_C, tl_C_size, random, 
                                                        total_time_samples);
                        tApt = T_A_plus_tripl(tl_A, tl_A_size, to_shift, 
                                            tl_C_size, total_time_samples, Dt);
                        shifted_res_arr[shift] = STTC_AB_C(tl_A, tl_A_size, 
                                                    tl_B, tl_B_size, to_shift, 
                                                    tl_C_size, Dt, tBm, tApt);
                        if (shifted_res_arr[shift] == 2.0) {
                            --denominator;
                        }
                        else {
                            mean += shifted_res_arr[shift];
                        }
                    }
                    mean /= denominator;
                    double st_dev = std_STTC_dir(shifted_res_arr, 
                                                        circ_shifts_num, mean);
                    double threshold = sign_thresh(mean, st_dev);
                    if ( trip_sttc > threshold) {
                        #pragma omp atomic
                        ++ttl_sgnfcnt_triplets;
                        #pragma omp atomic
                        ++motifs_sgnfcnts[pos];
                        sort(shifted_res_arr, (shifted_res_arr + 
                                                            circ_shifts_num));
                        pos = 0; 
                        while (pos < denominator && 
                                        shifted_res_arr[pos] <= trip_sttc) {
                            ++pos;
                        }
                        int b_real = map[b];
                        double percentile = pos / double(denominator);
                        #pragma omp critical
                        triplets<<a_real + 1<<','<<b_real + 1<<','<<c_real + 1
                                    <<','<<trip_sttc<<','<<percentile<<'\n';
                    }
                }
            }
            free(to_shift);
        }
    }
    triplets.close();
    info<<"\nNumber of total significant triplets: "<<ttl_sgnfcnt_triplets
                    <<" ( "<<(ttl_sgnfcnt_triplets * 100 / double(neurons * 
                    (neurons - 1) * (neurons - 2)))<<"% )"<<endl;
    
    
// Free memory
    for (int neur = 0; neur < neurons; ++neur) {
        free(tl_array[neur]);
    }
    
// Print Motifs
    print_motifs(motifs_triplets, motifs_sgnfcnts, info, string(argv[3]), 
                                                            shifts_s, Dt_s);
    
// Close output files
    info.close();
    
    return 0;
}
