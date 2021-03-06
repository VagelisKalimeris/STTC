/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: triplets.cpp                                                     *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "triplets_STTC.hpp"

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
                const vector<int> &time_line_C, int total_time_samples, int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    unsigned int a = 0, c = 0;
    
    if(time_line_A.size() == 0 || time_line_C.size() == 0) {
        return T;
    }
    /* all spikes of A are before tiles of C */
    if (time_line_A.back() < time_line_C.front()) {
        return T;
    }
    /* all spikes of A are after tiles of C */
    if ((time_line_C.back() + Dt) < time_line_A.front()) {
        return T;
    }

    while ((a < time_line_A.size()) && (c < time_line_C.size())) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if ((time_line_A[a] >= time_line_C[c]) && 
                                    (time_line_A[a] <= (time_line_C[c] + Dt))) {
            /* check if last calculated tile is before spike of A */
            if (last < time_line_A[a]) {
                /* add Dt + 1 */
                s += Dt + 1;
            }
            else {
                /* add Dt + 1 - (tA'_prev + Dt + 1 - tA'_curr) */
                s += time_line_A[a] + Dt - last;
            }
            last = time_line_A[a] + Dt;
            a++;
        }
        /* spike of A is before tile of spike of C [tC, tC + Dt] */
        else if (time_line_A[a] < time_line_C[c]) {
            a++;
        }
        /* spike of A is after tile of spike of C [tC, tC + Dt] */
        else if (time_line_A[a] > (time_line_C[c] + Dt)) {
            c++;
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
* PURPOSE: The number of firing events of the reduced spike train A that      *
*           falls within the tiles Δt before the firing events of spike       *
*            train B.                                                         *
*                                                                             *
* RETURNS: A number(int) of events.                                           *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
int N_BminusA_CA(const vector<int> &time_line_A, 
        const vector<int> &time_line_B, const vector<int> &time_line_C, int Dt)
{
    int N = 0;
    unsigned int a = 0, b = 0, c = 0;
    
    if(time_line_A.size() == 0 || time_line_B.size() == 0 || 
                                                    time_line_C.size() == 0) {
        return N;
    }
    /* all spikes of A are before tiles of C */
    if(time_line_A.back() < time_line_C.front()) {
        return N;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C.back() + Dt) < time_line_A.front()) {
        return N;
    }
    /* all spikes of A are before tiles of B */
    if(time_line_A.back() < (time_line_B.front() - Dt)) {
        return N;
    }
    /* all spikes of A are after tiles of B */
    if(time_line_B.back() < time_line_A.front()) {
        return N;
    }
    
    while((a < time_line_A.size()) && (b < time_line_B.size()) && 
                                                    (c < time_line_C.size())) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_line_A[a] >= time_line_C[c]) && 
                                    (time_line_A[a] <= (time_line_C[c] + Dt))) {
            /* spike of A is within tile of spike of B [tB, tB + Dt] */
            if((time_line_A[a] >= (time_line_B[b] - Dt)) && 
                                        (time_line_A[a] <= time_line_B[b])) {
                N++;
                a++;
            }
            /* spike of A is before tile of spike of B [tB, tB + Dt] */
            else if(time_line_A[a] < (time_line_B[b] - Dt)) {
                a++;
            }
            /* spike of A is after tile of spike of B [tB, tB + Dt] */
            else if(time_line_A[a] > time_line_B[b]) {
                b++;
            }
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
    
    return N;
}


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
int N_AplusB_CA(const vector<int> &time_line_A, 
        const vector<int> &time_line_B, const vector<int> &time_line_C, int Dt)
{
    int N = 0;
    unsigned int a = 0, b = 0, c = 0;
    
    if(time_line_A.size() == 0 || time_line_B.size() == 0 || 
                                                    time_line_C.size() == 0) {
        return N;
    }
    /* all spikes of A are before tiles of C */
    if(time_line_A.back() < time_line_C.front()) {
        return N;
    }
    /* all spikes of A are after tiles of C */
    if((time_line_C.back() + Dt) < time_line_A.front()) {
        return N;
    }
    /* all spikes of B are before tiles of A */
    if(time_line_B.back() < time_line_A.front()) {
        return N;
    }
    /* all spikes of B are after tiles of A */
    if((time_line_A.back() + Dt) < time_line_B.front()) {
        return N;
    }
    
    while((a < time_line_A.size()) && (b < time_line_B.size()) && 
                                                    (c < time_line_C.size())) {
        /* spike of A is within tile of spike of C [tC, tC + Dt] */
        if((time_line_A[a] >= time_line_C[c]) && 
                                    (time_line_A[a] <= (time_line_C[c] + Dt))) {
            /* spike of B is within tile of spike of A [tA, tA + Dt] */
            if((time_line_B[b] >= time_line_A[a]) && 
                                    (time_line_B[b] <= (time_line_A[a] + Dt))) {
                N++;
                b++;
            }
            /* spike of B is before tile of spike of A [tA, tA + Dt] */
            else if(time_line_B[b] < time_line_A[a]) {
                b++;
            }
            /* spike of B is after tile of spike of A [tA, tA + Dt] */
            else if(time_line_B[b] > (time_line_A[a] + Dt)) {
                a++;
            }
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
double STTC_AB_C(const vector<int> &time_line_A, const vector<int> &time_line_B,
                const vector<int> &time_line_C, int total_time_samples, int Dt) 
{
    int nBmACA =  N_BminusA_CA(time_line_A, time_line_B, time_line_C, Dt);
    int nApBCA = N_AplusB_CA(time_line_A, time_line_B, time_line_C, Dt);
    double tBm = T_B_minus(time_line_B, total_time_samples, Dt);
    double tApt = T_A_plus_tripl(time_line_A, time_line_C, 
                                                        total_time_samples, Dt);
    double nA = double(time_line_A.size()), nB = double(time_line_B.size());
    
    return 0.5 * ((((nBmACA / nA) - tBm) / (1.0 - ((nBmACA / nA) * tBm))) + 
                    (((nApBCA / nB) - tApt) / (1.0 - ((nApBCA / nB) * tApt))));
}