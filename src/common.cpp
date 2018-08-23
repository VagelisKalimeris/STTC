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
* FUNCTION NAME: T_A_plus                                                     *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), and a time           *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the sum of time tiles after a neuron's firing, divided  *
*           by the total time.                                                *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double T_A_plus(const vector<int> &time_line_A, int total_time_samples, 
                                                                        int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    
    if(time_line_A.size() == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_A.size() / double(total_time_samples);
    }
    else {
        for(unsigned int a = 0; a < time_line_A.size(); ++a) {
            /* check if last calculated tile is before tile of spike of A */
            if(last < time_line_A[a]) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add Dt + 1 - (tA'_prev + Dt + 1 - tA'_curr) */
                s += Dt + time_line_A[a] - last;
            }
            last = time_line_A[a] + Dt;
        }
        if((last != -1) && (last >= total_time_samples)) {
            s -= last + 1 - total_time_samples;
        }

        T = s / double(total_time_samples);
    }
    
    return T;
}


/******************************************************************************
* FUNCTION NAME: T_B_minus                                                    *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), and a time           *
*             interval(int).                                                  *
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
    
    if(time_line_B.size() == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_B.size() / double(total_time_samples);
    }
    else {
        for(unsigned int b = 0; b < time_line_B.size(); ++b) {
            /* check if last calculated tile is before tile of spike of B */
            if(last < (time_line_B[b] - Dt)) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add tB'_curr - tB'_prev */
                s += time_line_B[b] - last;
            }
            last = time_line_B[b];
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
void circular_shift(vector<int> &time_line, unsigned int random, 
                                                      int total_time_samples) {
    vector<int>::iterator front_it = time_line.begin();

    for (int i = 0; i < static_cast<int> (time_line.size()); i++) {
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
unsigned int random_gen(unsigned int max_number) {
    return 1 + rand() % max_number;
}

/*{
    auto machine = uniform_int_distribution<unsigned int>(0, max_number);
    return machine(mt19937(random_device()()));
}*/