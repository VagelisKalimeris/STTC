/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: tuplets.cpp                                                      *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "tuplets_STTC.hpp"

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
double P_A_B_minus(const vector<int> &time_line_A, 
                                        const vector<int> &time_line_B, int Dt)
{
    double P = 0.0;
    int N = 0;
    unsigned int a = 0, b = 0;
    
    unsigned int time_line_A_size = time_line_A.size();
    unsigned int time_line_B_size = time_line_B.size();
    if(time_line_A_size == 0 || time_line_B_size == 0) {
        return P;
    }
    /* all spikes of A are before tiles of B */
    if(time_line_A.back() < (time_line_B.front() - Dt)) {
        return P;
    }
    /* all spikes of A are after tiles of B */
    if(time_line_B.back() < time_line_A.front()) {
        return P;
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
    
    P = N / double(time_line_A_size);
    
    return P;
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
double P_B_A_plus(const vector<int> &time_line_A,
                                        const vector<int> &time_line_B, int Dt)
{
    double P = 0.0;
    int N = 0;
    unsigned int a = 0, b = 0;
    
    unsigned int time_line_A_size = time_line_A.size();
    unsigned int time_line_B_size = time_line_B.size();
    if(time_line_A_size == 0 || time_line_B_size == 0) {
        return P;
    }
    /* all spikes of B are before tiles of A */
    if(time_line_B.back() < time_line_A.front()) {
        return P;
    }
    /* all spikes of B are after tiles of A */
    if((time_line_A.back() + Dt) < time_line_B.front()) {
        return P;
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
    
    P = N / double(time_line_B_size);
    
    return P;
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
                        int total_time_samples, int Dt, double tBm, double tAp)
{
    double pABm = P_A_B_minus(time_line_A, time_line_B, Dt);
    double pBAp = P_B_A_plus(time_line_A, time_line_B, Dt);
    
    return 0.5 * (((pABm - tBm) / (1.0 - (pABm * tBm))) + 
                                        ((pBAp - tAp) / (1.0 - (pBAp * tAp))));
}
