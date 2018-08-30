/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.cpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "common.hpp"

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
                        const int total_neurons, const vector<int> &astrocytes)
{
    int total_firings = 0, max_neuron = -1, min_neuron = -1;
    unsigned int max = 0, min = 100000;
    const int astrocytes_size = astrocytes.size();
    const int neur_clean = total_neurons - astrocytes_size;
    int astro = 0, astrocyte = astrocytes[0];

    cout<<"\nThe data structure: "<<endl;
    for (int neur = 0; neur < total_neurons; neur++) {
        unsigned int time_line_size = spike_trains[neur].size();
        int pos;
        if (neur == astrocyte) {
            pos = neur_clean + astro;
            cout<<"No "<<neur + 1<<" astrocyte neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            astrocyte = astrocytes[(++astro) % astrocytes_size];
        }
        else {
            pos = neur - astro;
            cout<<"No "<<neur + 1<<" neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            total_firings += time_line_size;
            if (time_line_size > max) {
                max = time_line_size;
                max_neuron = neur + 1;
            }
            else if (time_line_size < min) {
                min = time_line_size;
                min_neuron = neur + 1;
            }
        }
        for (unsigned int fire = 0; fire < time_line_size; fire++) {
            cout<<spike_trains[pos][fire] + 1<<' ';
        }
        cout<<endl<<endl;
    }
    cout<<"\nNeurons' info without astrocytes:"<<endl;
    cout<<"Total number of spikes: "<<total_firings<<endl;
    cout<<"Neuron "<<max_neuron<<" has max spikes: "<<max<<endl;
    cout<<"Neuron "<<min_neuron<<" has min spikes: "<<min<<endl;
    cout<<"Average spikes in each neuron are: "
                                <<total_firings / double(neur_clean)<<endl;
}

/* comments for later */
void print_sgnfcnt_tuplet_begin(void)
{
    cout<<"\n+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(15)<<'+'
                                                        <<setw(13)<<'+'<<endl;
    cout<<"| Neuron A | Neuron B |"<<setfill(' ')<<setw(15)<<"STTC |"
                                            <<setw(13)<<"Percentile |"<<endl;
    cout<<"+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(15)<<'+'
                                                        <<setw(13)<<'+'<<endl;
}

void print_sgnfcnt_tuplet(const int neuron_A, const int neuron_B, 
                                    const double STTC, const double percentile)
{
    cout<<"| "<<setfill(' ')<<setw(8)<<neuron_A<<" | "<<setw(8)<<neuron_B
        <<" | "<<setw(12)<<STTC<<" | "<<setw(10)<<percentile<<" |"<<endl;
}

void print_sgnfcnt_tuplet_end(void)
{
    cout<<"+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(15)<<'+'
                                                        <<setw(13)<<'+'<<endl;
}

void print_sgnfcnt_triplet_begin(void)
{
    cout<<"\n+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(11)<<'+'
                                        <<setw(15)<<'+'<<setw(13)<<'+'<<endl;
    cout<<"| Neuron A | Neuron B | Neuron C |"<<setfill(' ')
                        <<setw(15)<<"STTC |"<<setw(13)<<"Percentile |"<<endl;
    cout<<"+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(11)<<'+'
                                        <<setw(15)<<'+'<<setw(13)<<'+'<<endl;
}

void print_sgnfcnt_triplet(const int neuron_A, const int neuron_B, 
                const int neuron_C, const double STTC, const double percentile)
{
    cout<<"| "<<setfill(' ')<<setw(8)<<neuron_A<<" | "
                    <<setw(8)<<neuron_B<<" | "<<setw(8)<<neuron_C<<" | "
                    <<setw(12)<<STTC<<" | "<<setw(10)<<percentile<<" |"<<endl;
}

void print_sgnfcnt_triplet_end(void)
{
    cout<<"+"<<setfill('-')<<setw(11)<<'+'<<setw(11)<<'+'<<setw(11)<<'+'
                                        <<setw(15)<<'+'<<setw(13)<<'+'<<endl;
}
