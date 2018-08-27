/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.hpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

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
static double T_B_minus(const vector<int> &time_line_B, int total_time_samples, 
                                        int Dt) __attribute__((always_inline));

static double T_B_minus(const vector<int> &time_line_B, int total_time_samples, 
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
static double sign_thresh(double mean, double st_dev) 
                                                __attribute__((always_inline));

static double sign_thresh(double mean, double st_dev)
{
    return mean + (3.0 * st_dev);
}


/******************************************************************************
* FUNCTION NAME: circular_shift                                               *
*                                                                             *
* ARGUMENTS: A vector representing the firings of a neuron, and a random      *
*             number between 0 and the maximum number of time events.         *
*                                                                             *
* PURPOSE: Shifts forward each firing of a neuron by a random number.         *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
static void circular_shift(vector<int> &time_line, unsigned int random, 
                        int total_time_samples) __attribute__((always_inline));

static void circular_shift(vector<int> &time_line, unsigned int random, 
                                                      int total_time_samples)
{
    vector<int>::iterator front_it = time_line.begin();

    for (unsigned int i = 0; i < time_line.size(); i++) {
        int temp = time_line[i] + random;

        if ((temp) < total_time_samples) {
            time_line[i] = temp;
        }
        else {
            time_line.erase(time_line.begin() + i);
            temp = temp - total_time_samples;
            time_line.insert(front_it, temp);
            ++front_it;
        }
    }
}

// Helper function. Generates random integers 
// in the range 0 - (total_time_samples-1).
static unsigned int random_gen(unsigned int max_number) 
                                                __attribute__((always_inline));

static unsigned int random_gen(unsigned int max_number)
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
void print_all_spikes(const vector<int> spike_trains[], int total_neurons);

/* comments for later */
void print_sgnfcnt_tuplet_begin(void);

void print_sgnfcnt_tuplet(const int neuron_A, const int neuron_B, 
                                    const double STTC, const double percentile);

void print_sgnfcnt_tuplet_end(void);

void print_sgnfcnt_triplet_begin(void);

void print_sgnfcnt_triplet(const int neuron_A, const int neuron_B, 
                const int neuron_C, const double STTC, const double percentile);

void print_sgnfcnt_triplet_end(void);
