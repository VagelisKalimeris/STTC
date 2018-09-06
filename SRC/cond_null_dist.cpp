/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: cond_null_dist.cpp                                               *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/cond_null_dist.hpp"

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
bool sign_trpl_limit(const int time_line_A[], int time_line_A_size, 
                        const int time_line_C[], int time_line_C_size, int Dt)
{
    int s = 0;
    int a = 0, c = 0;
    
    if(time_line_A_size == 0 || time_line_C_size == 0) {
        return false;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_C = time_line_C[0];
    /* all spikes of A are before tiles of C */
    if(time_line_A[time_line_A_size - 1] < time_stamp_C) {
        return false;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C[time_line_C_size - 1] + Dt) < time_stamp_A) {
        return false;
    }
    
    while((a < time_line_A_size) && (c < time_line_C_size)) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_stamp_A >= time_stamp_C) && 
                                    (time_stamp_A <= (time_stamp_C + Dt))) {
            ++s;
            time_stamp_A = time_line_A[++a];
        }
        /* spike of A is before tile of spike of C [tC, tC + Dt] */
        else if(time_stamp_A < time_stamp_C) {
            time_stamp_A = time_line_A[++a];
        }
        /* spike of A is after tile of spike of C [tC, tC + Dt] */
        else if(time_stamp_A > (time_stamp_C + Dt)) {
            time_stamp_C = time_line_C[++c];
        }
    }
    
    return (s > 5);
}