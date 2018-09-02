/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.cpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/common.hpp"

/******************************************************************************
* FUNCTION NAME: T_B_minus                                                    *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), the total time       *
*             samples recorded(int) and a time interval(int).                 *
*                                                                             *
* PURPOSE: Calculates the sum of time tiles before a neuron's firing, divided *
*           by the total time.                                                *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double T_B_minus(const vector<int> &time_line_B, int total_time_samples, 
                                                                        int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    
    unsigned int time_line_B_size = time_line_B.size();
    if(time_line_B_size == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_B_size / double(total_time_samples);
    }
    else {
        for(unsigned int b = 0; b < time_line_B_size; ++b) {
            int time_stamp_B = time_line_B[b];
            /* check if last calculated tile is before tile of spike of B */
            if(last < (time_stamp_B - Dt)) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add tB'_curr - tB'_prev */
                s += time_stamp_B - last;
            }
            last = time_stamp_B;
        }

        T = s / double(total_time_samples);
    }
    
    return T;
}


/******************************************************************************
* FUNCTION NAME: sign_thresh_A_B                                              *
*                                                                             *
* ARGUMENTS: The mean(double) and the standard(double) deviations of the      *
*             circ. shifted spike trains.                                     *
*                                                                             *
* RETURNS: The significant threshhold.                                        *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double sign_thresh(double mean, double st_dev)
{
    return mean + (3.0 * st_dev);
}


/******************************************************************************
* FUNCTION NAME: circular_shift                                               *
*                                                                             *
* ARGUMENTS: A vector representing the spikes of a neuron, and a random       *
*             number between 0 and the maximum number of time events.         *
*                                                                             *
* PURPOSE: Shifts forward each firing of a neuron by a random number.         *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void circular_shift(vector<int> &time_line, unsigned int random, 
                                                        int total_time_samples)
{
    unsigned int tl_size = time_line.size();
    
    for (unsigned int i = 0; i < tl_size; i++) {
        time_line[i] = (time_line[i] + random) % total_time_samples;
    }
    sort(time_line.begin(), (time_line.begin() + tl_size));
}



// Helper function. Generates random integers 
// in the range [1, total_time_samples].
unsigned int random_gen(unsigned int max_number)
{
    return 1 + rand() % max_number;
}


/******************************************************************************
* FUNCTION NAME: print_all_spikes                                             *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), the total number of  *
*             neurons(int).                                                   *   
*                                                                             *  
* PURPOSE: Prints all the spikes of each neuron of the dataset, as well as    *
*           the total number of spikes.                                       *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: See PURPOSE.                                                           *
*                                                                             *
******************************************************************************/
void print_all_spikes(const vector<int> spike_trains[], 
                        const int total_neurons, const vector<int> &astrocytes, 
                        ofstream &info, const string output)
{
    int total_firings = 0, max = 0, min = 100000;
    const int astro_size = astrocytes.size();
    const int neur_clean = total_neurons - astro_size;
    int astro = 0, astrocyte = astrocytes[0];
    vector<int> time_lines;
    
    ofstream spikes;
    spikes.open(("RESULTS/" + output + "_neurons_spikes.txt").c_str());
    if (!spikes.is_open()) {
        cout<<"Error opening results neurons spikes file!"<<endl;
        return;
    }
    spikes<<"\nThe data structure: "<<endl;
    for (int neur = 0; neur < total_neurons; neur++) {
        int pos, time_line_size;
        if (neur == astrocyte) {
            pos = neur_clean + astro;
            time_line_size = spike_trains[pos].size();
            spikes<<"No "<<neur + 1<<" astrocyte neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            astrocyte = astrocytes[(++astro) % astro_size];
        }
        else {
            pos = neur - astro;
            time_line_size = spike_trains[pos].size();
            spikes<<"No "<<neur + 1<<" neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            total_firings += time_line_size;
            if (time_line_size > max) {
                max = time_line_size;
            }
            else if (time_line_size < min) {
                min = time_line_size;
            }
            time_lines.push_back(time_line_size);
        }
        for (int fire = 0; fire < time_line_size; fire++) {
            spikes<<spike_trains[pos][fire] + 1<<' ';
        }
        spikes<<endl<<endl;
    }
    spikes.close();
    
    double median;
    sort(time_lines.begin(), (time_lines.begin() + neur_clean));
    if (neur_clean % 2) {
        median = time_lines[neur_clean/2]/1.0;
    }
    else {
        median = (time_lines[neur_clean/2 - 1] + time_lines[neur_clean/2])/2.0;
    }
    
    info<<"\nNeurons' info without astrocytes:"<<endl;
    info<<"Total number of spikes: "<<total_firings<<endl;
    info<<"Max spikes for neurons: "<<max<<endl;
    info<<"Min spikes for neurons: "<<min<<endl;
    info<<"Average spikes in each neuron: "
                                    <<total_firings / double(neur_clean)<<endl;
    info<<"Median spikes for neurons: "<<median<<endl;
}
