/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.hpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>

#include "common.hpp"
#include "tuplets_STTC.hpp"
#include "triplets_STTC.hpp"

using namespace std;

/******************************************************************************
* FUNCTION NAME: circ_STTC_A_B_C                                              *
*                                                                             *
* ARGUMENTS: A pre-existing array to store the results. Three neuron's        *
*             timelines(references to vectors), and a time interval(int).     *
*                                                                             *
* PURPOSE: Calculates circ_shifts_num random STTC values.                     *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void circ_STTC_A_B(double results_arr[], int circ_shifts_num, 
                const vector<int> &time_line_A, const vector<int> &time_line_B, 
                const vector<int> &time_line_C, int total_time_samples, int Dt);


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
