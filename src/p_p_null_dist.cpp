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
* PURPOSE: Calculates SHIFTS_NUM random STTC values.                          *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void circ_STTC_A_B(double results_arr[SHIFTS_NUM], 
	    const vector<int> &time_line_A, const vector<int> &time_line_B, int Dt)
{
	for (int i = 0; i < SHIFTS_NUM; i++) {
		vector<int> to_shift = time_line_A;

		circular_shift(to_shift, SHIFTS_NUM);
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
double mean_STTC_dir(double const arr[SHIFTS_NUM]) {
	double sum = 0.0;

	for (int i = 0; i < SHIFTS_NUM; i++) {
		sum += arr[i];
	}
	return sum / (double) SHIFTS_NUM;
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
double std_STTC_dir(double const arr[SHIFTS_NUM]) {
	double mean = mean_STTC_dir(arr), st_dev = 0.0;

		for (int i = 0; i < SHIFTS_NUM; i++) {
			st_dev += pow((arr[i] - mean), 2) / SHIFTS_NUM;
		}
	return sqrt(st_dev);
}

