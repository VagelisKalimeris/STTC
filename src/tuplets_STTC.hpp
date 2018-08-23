/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: tuplets.hpp                                                      *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include <iostream>
#include <cmath>
#include <vector>

#include "common.hpp"

using namespace std;

/******************************************************************************
* FUNCTION NAME: P_A_B_minus                                                  *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the fraction of the number of the firing events of A    *
*           which fall within Δt before each firing event of B by the number  *
*            of firing events of A.                                           *
*                                                                             *
* RETURNS: A double > 0 and < 1.                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double P_A_B_minus(const vector<int> &time_line_A,
                                        const vector<int> &time_line_B, int Dt);


/******************************************************************************
* FUNCTION NAME: P_B_A_plus                                                   *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the fraction of the number of the firing events of B    *
*           which fall within Δt after each firing event of A by the number   *
*            of firing events of A.                                           *
*                                                                             *
* RETURNS: A double > 0 and < 1.                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double P_B_A_plus(const vector<int> &time_line_A,
                                        const vector<int> &time_line_B, int Dt);


/******************************************************************************
* FUNCTION NAME: STTC_A_B                                                     *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the the correlation between spike trains for the spikes *
*           of B that follows spikes of A and the spikes of A that proceeds   *
*            spikes of B.                                                     *
*                                                                             *
* RETURNS: The STTC value(double) of the pair A,B.                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double STTC_A_B(const vector<int> &time_line_A, const vector<int> &time_line_B,
                        int total_time_samples, int Dt, double tBm, double tAp);
