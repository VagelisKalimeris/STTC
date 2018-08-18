/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: p_p_null_dist.cpp                                                *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>
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
	    const vector<int> &time_line_A, const vector<int> &time_line_B, int Dt)
{
	for (int i = 0; i < circ_shifts_num; i++) {
		vector<int> to_shift = time_line_A;

		circular_shift(to_shift, circ_shifts_num);
		results_arr[i] = STTC_A_B(to_shift, time_line_B, Dt);
	}
}


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
double mean_STTC_dir(double const arr[], int circ_shifts_num) {
	double sum = 0.0;

	for (int i = 0; i < circ_shifts_num; i++) {
		sum += arr[i];
	}
	return sum / (double) circ_shifts_num;
}


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
double std_STTC_dir(double const arr[], int circ_shifts_num) {
	double mean = mean_STTC_dir(arr, circ_shifts_num), st_dev = 0.0;

		for (int i = 0; i < circ_shifts_num; i++) {
			st_dev += pow((arr[i] - mean), 2) / circ_shifts_num;
		}
	return sqrt(st_dev);
}

