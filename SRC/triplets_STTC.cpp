/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: triplets.cpp                                                     *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/triplets_STTC.hpp"

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
double T_A_plus_tripl(const int time_line_A[], int time_line_A_size, 
                                const int time_line_C[], int time_line_C_size, 
                                int total_time_samples, int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    int a = 0, c = 0;
    
    if(time_line_A_size == 0 || time_line_C_size == 0) {
        return T;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_C = time_line_C[0];
    /* all spikes of A are before tiles of C */
    if (time_line_A[time_line_A_size - 1] < time_stamp_C) {
        return T;
    }
    /* all spikes of A are after tiles of C */
    if ((time_line_C[time_line_C_size - 1] + Dt) < time_stamp_A) {
        return T;
    }
    
    while((a < time_line_A_size) && (c < time_line_C_size)) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_stamp_A >= time_stamp_C) && 
                                    (time_stamp_A <= (time_stamp_C + Dt))) {
            /* check if last calculated tile is before spike of A */
            if (last < time_stamp_A) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add Dt + 1 - (tA'_prev + Dt + 1 - tA'_curr) */
                s += time_stamp_A + Dt - last;
            }
            last = time_stamp_A + Dt;
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
    if((last != -1) && (last >= total_time_samples)){
        s -= last + 1 - total_time_samples;
    }

    T = s / double(total_time_samples);

    return T;
}


/******************************************************************************
* FUNCTION NAME: N_BminusA_CA                                                 *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
*                                                                             *
* PURPOSE: Finds the number of firing events of the reduced spike train A     *
*           that fall within the tiles Δt before the firing events of         *
*            spike train B.                                                   *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int N_BminusA_CA(const int time_line_A[], int time_line_A_size, 
                        const int time_line_B[], int time_line_B_size, 
                        const int time_line_C[], int time_line_C_size, int Dt)
{
    int N = 0;
    int a = 0, b = 0, c = 0;
    
    if(time_line_A_size == 0 || time_line_B_size == 0 || 
                                                    time_line_C_size == 0) {
        return N;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0], 
                                                time_stamp_C = time_line_C[0];
    /* all spikes of A are before tiles of C */
    if(time_line_A[time_line_A_size - 1] < time_stamp_C) {
        return N;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C[time_line_C_size - 1] + Dt) < time_stamp_A) {
        return N;
    }
    /* all spikes of A are before tiles of B */
    if(time_line_A[time_line_A_size - 1] < (time_stamp_B - Dt)) {
        return N;
    }
    /* all spikes of A are after tiles of B */
    if(time_line_B[time_line_B_size - 1] < time_stamp_A) {
        return N;
    }
    
    while((a < time_line_A_size) && (b < time_line_B_size) && 
                                                    (c < time_line_C_size)) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_stamp_A >= time_stamp_C) && 
                                    (time_stamp_A <= (time_stamp_C + Dt))) {
            /* spike of A is within tile of spike of B [tB, tB + Dt] */
            if((time_stamp_A >= (time_stamp_B - Dt)) && 
                                            (time_stamp_A <= time_stamp_B)) {
                ++N;
                time_stamp_A = time_line_A[++a];
            }
            /* spike of A is before tile of spike of B [tB, tB + Dt] */
            else if(time_stamp_A < (time_stamp_B - Dt)) {
                time_stamp_A = time_line_A[++a];
            }
            /* spike of A is after tile of spike of B [tB, tB + Dt] */
            else if(time_stamp_A > time_stamp_B) {
                time_stamp_B = time_line_B[++b];
            }
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
    
    return N;
}


/******************************************************************************
* FUNCTION NAME: N_AplusB_CA                                                  *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
*                                                                             *
* PURPOSE: Finds the number of firing events of the spike train B that        *
*           fall within the tiles Δt after the firing events of reduced       *
*            spike train A.                                                   *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int N_AplusB_CA(const int time_line_A[], int time_line_A_size, 
                        const int time_line_B[], int time_line_B_size, 
                        const int time_line_C[], int time_line_C_size, int Dt)
{
    int N = 0;
    int a = 0, b = 0, c = 0;
    
    if(time_line_A_size == 0 || time_line_B_size == 0 || 
                                                    time_line_C_size == 0) {
        return N;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0], 
                                                time_stamp_C = time_line_C[0];
    /* all spikes of A are before tiles of C */
    if(time_line_A[time_line_A_size - 1] < time_stamp_C) {
        return N;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C[time_line_C_size - 1] + Dt) < time_stamp_A) {
        return N;
    }
    /* all spikes of B are before tiles of A */
    if(time_line_B[time_line_B_size - 1] < time_stamp_A) {
        return N;
    }
    /* all spikes of B are after tiles of A */
    if((time_line_A[time_line_A_size - 1] + Dt) < time_stamp_B) {
        return N;
    }
    
    while((a < time_line_A_size) && (b < time_line_B_size) && 
                                                    (c < time_line_C_size)) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_stamp_A >= time_stamp_C) && 
                                    (time_stamp_A <= (time_stamp_C + Dt))) {
            /* spike of B is within tile of spike of A [tA, tA + Dt] */
            if((time_stamp_B >= time_stamp_A) && 
                                    (time_stamp_B <= (time_stamp_A + Dt))) {
                ++N;
                time_stamp_B = time_line_B[++b];
            }
            /* spike of B is before tile of spike of A [tA, tA + Dt] */
            else if(time_stamp_B < time_stamp_A) {
                time_stamp_B = time_line_B[++b];
            }
            /* spike of B is after tile of spike of A [tA, tA + Dt] */
            else if(time_stamp_B > (time_stamp_A + Dt)) {
                time_stamp_A = time_line_A[++a];
            }
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
    
    return N;
}


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
double STTC_AB_C(const int time_line_A[], int time_line_A_size, 
                const int time_line_B[], int time_line_B_size, 
                const int time_line_C[], int time_line_C_size, 
                int Dt, double tBm, double tApt)
{
    int nBmACA =  N_BminusA_CA(time_line_A, time_line_A_size, time_line_B, 
                        time_line_B_size, time_line_C, time_line_C_size, Dt);
    int nApBCA = N_AplusB_CA(time_line_A, time_line_A_size, time_line_B, 
                        time_line_B_size, time_line_C, time_line_C_size, Dt);
    double nA = double(time_line_A_size), nB = double(time_line_B_size);
    
    if ((nBmACA == nA && tBm == 1.0) || (nApBCA == nB && tApt == 1.0)) {
        return 2.0;
    }
    
    return 0.5 * ((((nBmACA / nA) - tBm) / (1.0 - ((nBmACA / nA) * tBm))) + 
                    (((nApBCA / nB) - tApt) / (1.0 - ((nApBCA / nB) * tApt))));
}
