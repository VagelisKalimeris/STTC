/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.cpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include "cond_null_dist.hpp"

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
    int s = 0, a = 0, c = 0;
    
    /* all spikes of A are before tiles of C */
    if(time_line_A.back() < time_line_C.front()) {
        return false;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C.back() + Dt) < time_line_A.front()) {
        return false;
    }
    
    while((a < time_line_A.size()) && (c < time_line_C.size())) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_line_A[a] >= time_line_C[c]) && 
                                    (time_line_A[a] <= (time_line_C[c] + Dt))) {
            s++;
            a++;
        }
        /* spike of A is before tile of spike of C [tC, tC + Dt] */
        else if(time_line_A[a] < time_line_C[c]) {
            a++;
        }
        /* spike of A is after tile of spike of C [tC, tC + Dt] */
        else if(time_line_A[a] > (time_line_C[c] + Dt)) {
            c++;
        }
    }
    
    return (s > 5);
}