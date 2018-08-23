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
