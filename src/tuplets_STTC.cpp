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
double T_A_plus(const vector<int> &time_line_A, int total_time_samples, 
                                                                        int Dt)
{
    double T = 0.0;
    int s = 0, last = -1;
    
    unsigned int time_line_A_size = time_line_A.size();
    if(time_line_A_size == 0) {
        return T;
    }
    if(Dt == 0) {
        T = time_line_A_size / double(total_time_samples);
    }
    else {
        for(unsigned int a = 0; a < time_line_A_size; ++a) {
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
* FUNCTION NAME: P_A_B_minus                                                  *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the fraction of the number of the firing events of A    *
*           which fall within Δt before each firing event of B by the number  *
*            of firing events of A.                                           *
*                                                                             *
* RETURNS: A double >= 0 and <= 1.                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double N_A_B_minus(const vector<int> &time_line_A, 
                                        const vector<int> &time_line_B, int Dt)
{
    int N = 0;
    unsigned int a = 0, b = 0;
    
    unsigned int time_line_A_size = time_line_A.size();
    unsigned int time_line_B_size = time_line_B.size();
    if(time_line_A_size == 0 || time_line_B_size == 0) {
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
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0];
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
* PURPOSE: Calculates the fraction of the number of the firing events of B    *
*           which fall within Δt after each firing event of A by the number   *
*            of firing events of B.                                           *
*                                                                             *
* RETURNS: A double >= 0 and <= 1.                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double N_B_A_plus(const vector<int> &time_line_A,
                                        const vector<int> &time_line_B, int Dt)
{
    int N = 0;
    unsigned int a = 0, b = 0;
    
    unsigned int time_line_A_size = time_line_A.size();
    unsigned int time_line_B_size = time_line_B.size();
    if(time_line_A_size == 0 || time_line_B_size == 0) {
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
    
    int time_stamp_A = time_line_A[0], time_stamp_B = time_line_B[0];
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
double STTC_A_B(const vector<int> &time_line_A, const vector<int> &time_line_B,
                                                int Dt, double tBm, double tAp)
{
    double nABm = N_A_B_minus(time_line_A, time_line_B, Dt);
    double nBAp = N_B_A_plus(time_line_A, time_line_B, Dt);
    double nA = double(time_line_A.size()), nB = double(time_line_B.size());
    
    if (nA == 0 || nB == 0 || (nABm == nA && tBm == 1) || 
                                                    (nBAp == nB && tAp == 1)) {
        return 2.0;
    }
    return 0.5 * ((((nABm / nA) - tBm) / (1.0 - ((nABm / nA) * tBm))) + 
                        (((nBAp / nB) - tAp) / (1.0 - ((nBAp / nB) * tAp))));
}
