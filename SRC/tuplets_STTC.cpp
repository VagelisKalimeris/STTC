/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: tuplets.cpp                                                      *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/tuplets_STTC.hpp"

/******************************************************************************
* FUNCTION NAME: T_A_plus                                                     *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), the total time       *
*             samples recorded(int) and a time interval(int).                 *
*                                                                             *
* PURPOSE: Calculates the sum of time tiles after a neuron's firing, divided  *
*           by the total time.                                                *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double T_A_plus(const int time_line_A[], int time_line_A_size, 
                                                int total_time_samples, int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    
    if(time_line_A_size == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_A_size / double(total_time_samples);
    }
    else {
        for(int a = 0; a < time_line_A_size; ++a) {
            int time_stamp_A = time_line_A[a];
            /* check if last calculated tile is before tile of spike of A */
            if(last < time_stamp_A) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add Dt + 1 - (tA'_prev + Dt + 1 - tA'_curr) */
                s += Dt + time_stamp_A - last;
            }
            last = time_stamp_A + Dt;
        }
        if((last != -1) && (last >= total_time_samples)) {
            s -= last + 1 - total_time_samples;
        }

        T = s / double(total_time_samples);
    }
    
    return T;
}


/******************************************************************************
* FUNCTION NAME: N_A_B_minus                                                  *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Finds the number of firing events of the spike train A that fall   *
*           within the tiles Δt before the firing events of spike train B.    *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double N_A_B_minus(const int time_line_A[], int time_line_A_size, 
                        const int time_line_B[], int time_line_B_size, int Dt)
{
    int N = 0;
    int a = 0, b = 0;
    
    if(time_line_A_size == 0 || time_line_B_size == 0) {
        return N;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0];
    /* all spikes of A are before tiles of B */
    if(time_line_A[time_line_A_size - 1] < (time_stamp_B - Dt)) {
        return N;
    }
    /* all spikes of A are after tiles of B */
    if(time_line_B[time_line_B_size - 1] < time_stamp_A) {
        return N;
    }
    
    while((a < time_line_A_size) && (b < time_line_B_size)) {
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
    
    return N;
}


/******************************************************************************
* FUNCTION NAME: P_B_A_plus                                                   *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Finds the number of firing events of the spike train B that fall   *
*           within the tiles Δt after the firing events of spike train A.     *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double N_B_A_plus(const int time_line_A[], int time_line_A_size, 
                        const int time_line_B[], int time_line_B_size, int Dt)
{
    int N = 0;
    int a = 0, b = 0;
    
    if(time_line_A_size == 0 || time_line_B_size == 0) {
        return N;
    }
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0];
    /* all spikes of B are before tiles of A */
    if(time_line_B[time_line_B_size - 1] < time_stamp_A) {
        return N;
    }
    /* all spikes of B are after tiles of A */
    if((time_line_A[time_line_A_size - 1] + Dt) < time_stamp_B) {
        return N;
    }
    
    while((a < time_line_A_size) && (b < time_line_B_size)) {
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
    
    return N;
}


/******************************************************************************
* FUNCTION NAME: STTC_A_B                                                     *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), the total time    *
*             samples recorded(int)and a time interval(int).                  *
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
double STTC_A_B(const int time_line_A[], int time_line_A_size, 
                                const int time_line_B[], int time_line_B_size, 
                                int Dt, double tBm, double tAp)
{
    double nABm = N_A_B_minus(time_line_A, time_line_A_size, time_line_B, 
                                                        time_line_B_size, Dt);
    double nBAp = N_B_A_plus(time_line_A, time_line_A_size, time_line_B, 
                                                        time_line_B_size, Dt);
    double nA = double(time_line_A_size), nB = double(time_line_B_size);
    
    if ((nABm == nA && tBm == 1.0) || (nBAp == nB && tAp == 1.0)) {
        return 2.0;
    }
    
    return 0.5 * ((((nABm / nA) - tBm) / (1.0 - ((nABm / nA) * tBm))) + 
                        (((nBAp / nB) - tAp) / (1.0 - ((nBAp / nB) * tAp))));
}
