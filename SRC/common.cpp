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
double T_B_minus(const int time_line_B[], int time_line_B_size, 
                                                int total_time_samples, int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    
    if(time_line_B_size == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_B_size / double(total_time_samples);
    }
    else {
        for(int b = 0; b < time_line_B_size; ++b) {
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
void circular_shift(int to_shift[], const int time_line[], int tl_size, 
                                    unsigned int random, int total_time_samples)
{
    for (int i = 0; i < tl_size; i++) {
        to_shift[i] = (time_line[i] + random) % total_time_samples;
    }
    sort(to_shift, (to_shift + tl_size));
}



// Helper function. Generates random integers 
// in the range [1, total_time_samples].
unsigned int random_gen(unsigned int max_number)
{
    return 1 + rand() % max_number;
}
