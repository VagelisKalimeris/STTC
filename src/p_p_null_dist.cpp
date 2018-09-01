/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: p_p_null_dist.cpp                                                *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/p_p_null_dist.hpp"

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
double mean_STTC_dir(double const arr[], int circ_shifts_num)
{
    double sum = 0.0;

    for (int i = 0; i < circ_shifts_num; i++) {
        sum += arr[i];
    }
    return sum / (double) circ_shifts_num;
}


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
double std_STTC_dir(double const arr[], int circ_shifts_num)
{
    double mean = mean_STTC_dir(arr, circ_shifts_num), st_dev = 0.0;

        for (int i = 0; i < circ_shifts_num; i++) {
            double val = arr[i];
            st_dev += (val - mean) * (val - mean) / circ_shifts_num;
        }
    return sqrt(st_dev);
}
