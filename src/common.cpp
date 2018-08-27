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




// Helper function. Generates random integers 
// in the range [1, total_time_samples].



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
void print_all_spikes(const vector<int> spike_trains[], int total_neurons)
{
    int total_firings = 0, max_neuron = -1, min_neuron = -1;
    unsigned int max = 0, min = 100000;

    cout<<"\nThe data structure: "<<endl;
    for (int neur = 0; neur < total_neurons; neur++) {
        unsigned int time_line_size = spike_trains[neur].size();
        cout<<"No "<<neur + 1<<" neuron's spikes ("<<time_line_size<<"):\n";
        for (unsigned int fire = 0; fire < time_line_size; fire++) {
            cout<<spike_trains[neur][fire] + 1<<' ';
            total_firings++;
        }
        if (time_line_size > max) {
            max = time_line_size;
            max_neuron = neur + 1;
        }
        else if (time_line_size < min) {
            min = time_line_size;
            min_neuron = neur + 1;
        }
        cout<<endl<<endl;
    }
    cout<<"\nTotal number of spikes: "<<total_firings<<endl;
    cout<<"Neuron "<<max_neuron<<" has max spikes: "<<max<<endl;
    cout<<"Neuron "<<min_neuron<<" has min spikes: "<<min<<endl;
    cout<<"Average spikes in each neuron are: "<<total_firings/double(total_neurons)<<endl; 
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
