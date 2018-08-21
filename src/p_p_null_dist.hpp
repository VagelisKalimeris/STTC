/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: p_p_null_dist.hpp                                                *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include <vector>

#include "common.hpp"
#include "tuplets_STTC.hpp"

using namespace std;

/******************************************************************************
* FUNCTION NAME: circ_STTC_A_B                                                *
*                                                                             *
* ARGUMENTS: A pre-existing array to store the results. Two neuron's          *
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
        const vector<int> &time_line_A, const vector<int> &time_line_B, int Dt);


/******************************************************************************
* FUNCTION NAME: mean_STTC_dir                                                *
*                                                                             *
* ARGUMENTS: An array containing the shifted per pair STTC results.           *
*                                                                             *
* PURPOSE: Calculates the mean value.                                         *
*                                                                             *
* RETURNS: The mean value.                                                    *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double mean_STTC_dir(double const arr[], int circ_shifts_num);


/******************************************************************************
* FUNCTION NAME: std_STTC_dir                                                 *
*                                                                             *
* ARGUMENTS: An array containing the shifted per pair STTC results.           *
*                                                                             *
* PURPOSE: Calculates the standard deviation.                                 *
*                                                                             *
* RETURNS: The standard deviation.                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double std_STTC_dir(double const arr[], int circ_shifts_num);


// We also use the functions STTC_A_B, sign_thresh_A_B from common.hpp
