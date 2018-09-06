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

using namespace std;

/******************************************************************************
* FUNCTION NAME: mean_STTC_dir                                                *
*                                                                             *
* ARGUMENTS: An array containing the shifted per pair STTC results as well as *
*             the total number of circular shifts.                            *
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
* ARGUMENTS: An array containing the shifted per pair STTC results as well as *
*             the total number of circular shifts.                            *
*                                                                             *
* PURPOSE: Calculates the standard deviation.                                 *
*                                                                             *
* RETURNS: The standard deviation.                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double std_STTC_dir(double const arr[], int circ_shifts_num, const double mean);
