/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.hpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>

#include "common.hpp"
#include "tuplets_STTC.hpp"
#include "triplets_STTC.hpp"

using namespace std;

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
                                        const vector<int> &time_line_C, int Dt);
