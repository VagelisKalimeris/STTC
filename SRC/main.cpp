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
        int astro = 0;
        int astrocyte = astrocytes[0];
        for (int neur = 0; neur < neurons; ++neur) {
            while ((neur + astro) == astrocyte) {
                if (line[neur + astro] == '1') {
                    spike_trains[neurons + astro].push_back(total_time_samples);
                }
                astrocyte = astrocytes[(++astro) % astro_size];
            }
            if (line[neur + astro] == '1') {
                spike_trains[neur].push_back(total_time_samples);
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
    ofstream info;
    info.open(("RESULTS/" + string(argv[3]) + "_neurons_info.txt").c_str());
    if (!info.is_open()) {
        cout<<"Error opening results neurons info file!"<<endl;
        return 0;
    }
    print_all_spikes(spike_trains, neurons + astro_size, astrocytes, info, 
                                                            string(argv[3]));
    
// Start random sequence
    srand(time(NULL));
    
// All T for tuplets
    double T_Aplus[neurons];
    double T_Bminus[neurons];
    for(int neur = 0; neur < neurons; ++neur) {
        vector<int> time_line = spike_trains[neur];
        int time_line_size = time_line.size();
        T_Aplus[neur] = T_A_plus(&time_line[0], time_line_size, total_time_samples, Dt);
        T_Bminus[neur] = T_B_minus(&time_line[0], time_line_size, total_time_samples, Dt);
    }
    
// Significant tuplets
    bool sgnfcnt_tuplets[neurons][neurons];
    
// Significant limit
    bool sgnfcnt_limit[neurons][neurons];
    
// Reduced spiketrain T for triplets
    double T_Aplus_tripl[neurons][neurons];
    
// Calculate per pair STTC
    ofstream tuplets;
    tuplets.open(("RESULTS/" + string(argv[3]) + "_tuplets.csv").c_str());
    if (!tuplets.is_open()) {
        cout<<"Error opening results tuplets file!"<<endl;
        return 0;
    }
    tuplets<<"NeuronA,NeuronB,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        double tAp = T_Aplus[a];
        int a_real = map[a];
        int time_line_A_size = time_line_A.size();
        for (int b = 0; b < neurons; b++) { // Neuron B
        // It will be used to help in categorization of motifs
            sgnfcnt_tuplets[a][b] = false;
            if (a == b) {continue;} // Skip same neurons
            vector<int> time_line_B = spike_trains[b];
            int time_line_B_size = time_line_B.size();
            sgnfcnt_limit[a][b] = sign_trpl_limit(&time_line_A[0], 
                    time_line_A_size, &time_line_B[0], time_line_B_size, Dt);
            if (sgnfcnt_limit[a][b]) {
                T_Aplus_tripl[a][b] = T_A_plus_tripl(&time_line_A[0], 
                                    time_line_A_size, &time_line_B[0], 
                                    time_line_B_size, total_time_samples, Dt);
            }
            double tBm = T_Bminus[b];
            double tupl_sttc = STTC_A_B(&time_line_A[0], time_line_A_size, 
                            &time_line_B[0], time_line_B_size, Dt, tBm, tAp);
            if (tupl_sttc == 2.0) {continue;}
            int denominator = circ_shifts_num;
            double mean = 0;
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            for (int shift = 0; shift < circ_shifts_num; shift++) {
            // Shifted spike trains will be copied here
                vector<int> to_shift = time_line_A;
                unsigned int random = random_gen(total_time_samples);
                circular_shift(&to_shift[0], time_line_A_size, random, 
                                                        total_time_samples);
                tAp = T_A_plus(&to_shift[0], time_line_A_size, 
                                                    total_time_samples, Dt);
                shifted_res_arr[shift] = STTC_A_B(&to_shift[0], 
                                            time_line_A_size, &time_line_B[0], 
                                            time_line_B_size, Dt, tBm, tAp);
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
                tuplets<<a_real + 1<<','<<b_real + 1<<','<<tupl_sttc<<','
                                                            <<percentile<<endl;
            }
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
    triplets.open(("RESULTS/" + string(argv[3]) + "_triplets.csv").c_str());
    if (!triplets.is_open()) {
        cout<<"Error opening results triplets file!"<<endl;
        return 0;
    }
    triplets<<"NeuronA,NeuronB,NeuronC,STTC,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        vector<int> time_line_A = spike_trains[a];
        int a_real = map[a];
        int time_line_A_size = time_line_A.size();
        for (int c = 0; c < neurons; c++) { // Neuron C
            if (a == c) {continue;} // Skip same neurons
            vector<int> time_line_C = spike_trains[c];
            int time_line_C_size = time_line_C.size();
            bool sign_trplt_limit = sgnfcnt_limit[a][c];
            double tApt = T_Aplus_tripl[a][c];
            int c_real = map[c];
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (b == a || b == c) {continue;} // Skip same neurons
                int pos = sgnfcnt_tuplets[c][a] * 4 + 
                        sgnfcnt_tuplets[c][b] * 2 + sgnfcnt_tuplets[a][b] * 1;
                ++motifs_triplets[pos];
                if (!sign_trplt_limit) {
                    continue; // Reduced A spike train has < 5 spikes
                }
                vector<int> time_line_B = spike_trains[b];
                int time_line_B_size = time_line_B.size();
                double tBm = T_Bminus[b];
                double trip_sttc = STTC_AB_C(&time_line_A[0], 
                        time_line_A_size, &time_line_B[0], time_line_B_size, 
                        &time_line_C[0], time_line_C_size, Dt, tBm, tApt);
                if (trip_sttc == 2.0) {
                    --motifs_triplets[pos];
                    continue;
                }
                int denominator = circ_shifts_num;
                double mean = 0;
            // STTC values of shifted spike trains
                double shifted_res_arr[circ_shifts_num];
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                // Shifted spike trains will be copied here
                    vector<int> to_shift = time_line_C;
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(&to_shift[0], time_line_C_size, random, 
                                                        total_time_samples);
                    tApt = T_A_plus_tripl(&time_line_A[0], time_line_A_size, 
                                            &to_shift[0], time_line_C_size, 
                                            total_time_samples, Dt);
                    shifted_res_arr[shift] = STTC_AB_C(&time_line_A[0], 
                        time_line_A_size, &time_line_B[0], time_line_B_size, 
                        &to_shift[0], time_line_C_size, Dt, tBm, tApt);
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
                    ++ttl_sgnfcnt_triplets;
                    ++motifs_sgnfcnts[pos];
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    pos = 0; 
                    while (pos < denominator && 
                                        shifted_res_arr[pos] <= trip_sttc) {
                        ++pos;
                    }
                    int b_real = map[b];
                    double percentile = pos / double(denominator);
                    triplets<<a_real + 1<<','<<b_real + 1<<','<<c_real + 1<<','
                                            <<trip_sttc<<','<<percentile<<endl;
                }
            }
        }
    }
    triplets.close();
    info<<"\nNumber of total significant triplets: "<<ttl_sgnfcnt_triplets
                    <<" ( "<<(ttl_sgnfcnt_triplets * 100 / double(neurons * 
                    (neurons - 1) * (neurons - 2)))<<"% )"<<endl;
    
    
// Print Motifs
    print_motifs(motifs_triplets, motifs_sgnfcnts, info, string(argv[3]));
    
// Close output files
    info.close();
    
    return 0;
}
