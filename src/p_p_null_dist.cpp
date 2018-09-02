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
    double denominator = circ_shifts_num;

    for (int i = 0; i < circ_shifts_num; i++) {
        double val = arr[i];
        if (val == 2.0) {
            --denominator;
            continue;
        }
        sum += val;
    }
    return sum / denominator;
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
double std_STTC_dir(double const arr[], int circ_shifts_num, const double mean)
{
    double st_dev = 0.0;
    double denominator = circ_shifts_num;

    for (int i = 0; i < circ_shifts_num; i++) {
        double val = arr[i];
        if (val == 2.0) {
            --denominator;
            continue;
        }
        st_dev += (val - mean) * (val - mean);
    }
    return sqrt(st_dev / denominator);
}
