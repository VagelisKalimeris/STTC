/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: triplets.hpp                                                     *
*                                                                             *
*******************************************************************************
******************************************************************************/




#include <cmath>
#include<vector>

#include "common.hpp"

using namespace std;

/******************************************************************************
* FUNCTION NAME: T_A_plus_tripl                                               *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total number  *
*             of time samples recorded(int), and a time interval(int).        *
*                                                                             *
* PURPOSE: Calculates the fraction of the total recording time which is       *
*           covered by the tiles +Δt after each spike of A, that fall within  *
*            the tiles Δt after each spike of C.                              *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double T_A_plus_tripl(const vector<int> &time_line_A,
                                       const vector<int> &time_line_C, int Dt);



/******************************************************************************
* FUNCTION NAME: N_BminusA_CA                                                 *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
*                                                                             *
* PURPOSE: The number of firing events of the reduced spike train A that      *
*           falls within the tiles Δt before the firing events of spike       *
*            train B.                                                         *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int N_BminusA_CA(const vector<int> &time_line_A, const vector<int> &time_line_B,
                                        const vector<int> &time_line_C, int Dt);


/******************************************************************************
* FUNCTION NAME: N_AplusB_CA                                                  *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
*                                                                             *
* PURPOSE: Calculates the fraction of the total recording time which is       *
*           covered by the tiles +Δt after each spike of B, that fall within  *
*            the tiles Δt after each spike of B.                              *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int N_AplusB_CA(const vector<int> &time_line_A, const vector<int> &time_line_B,
                                       const vector<int> &time_line_C, int Dt);


/******************************************************************************
* FUNCTION NAME: STTC_AB_C                                                    *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
*                                                                             *
* PURPOSE: Estimates the temporal correlation of two neurons, given that a    *
*           third neuron is firing.                                           *
*                                                                             *
* RETURNS: The STTC value(double) of the pair A,B,C.                          *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double STTC_AB_C(const vector<int> &time_line_A, const vector<int> &time_line_B,
               const vector<int> &time_line_C, int total_time_samples, int Dt);
