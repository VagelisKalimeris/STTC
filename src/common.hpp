/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.hpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>
using namespace std;

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
                                                                       int Dt);


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
                                                                       int Dt);


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
void circular_shift(vector<int> &time_line, int random);


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
double sign_thresh(double mean, double st_dev);


/******************************************************************************
* FUNCTION NAME: sign_trpl_limit                                              *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates if the number of firing events of ‘reduced A’ is        *
*           greater than 5 or not.                                            *
*                                                                             *
* RETURNS: True or False.                                                     *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
bool sign_trpl_limit(const vector<int> &time_line_A, 
                                        const vector<int> &time_line_C, int Dt);
