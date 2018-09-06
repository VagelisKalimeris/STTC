/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: triplets.hpp                                                     *
*                                                                             *
*******************************************************************************
******************************************************************************/


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
double T_A_plus_tripl(const int time_line_A[], int time_line_A_size, 
                                const int time_line_C[], int time_line_C_size, 
                                int total_time_samples, int Dt);


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
                        const int time_line_C[], int time_line_C_size, int Dt);


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
                        const int time_line_C[], int time_line_C_size, int Dt);


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
                int Dt, double tBm, double tApt);
