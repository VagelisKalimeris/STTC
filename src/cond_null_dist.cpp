/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.cpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>
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
void circ_STTC_A_B_C(double results_arr[], int circ_shifts_num, 
                 const vector<int> &time_line_A, const vector<int> &time_line_B, 
                                         const vector<int> &time_line_C, int Dt)
{
	if (sign_trpl_limit(time_line_A, time_line_C)) { // NOT FINISHED!
		for (int i = 0; i < circ_shifts_num; i++) {
			vector<int> to_shift = time_line_C;

			circular_shift(to_shift, circ_shifts_num);
			results_arr[i] = STTC_AB_C(time_line_A, time_line_B, 
                                                                 to_shift, Dt);
		}
	}
}


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
                                        const vector<int> &time_line_C, int Dt)
{
// To Do..
}
