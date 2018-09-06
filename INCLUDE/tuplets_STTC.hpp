/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: tuplets.hpp                                                      *
*                                                                             *
*******************************************************************************
******************************************************************************/


using namespace std;

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
                                                int total_time_samples, int Dt);


/******************************************************************************
* FUNCTION NAME: P_A_B_minus                                                  *
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
                        const int time_line_B[], int time_line_B_size, int Dt);


/******************************************************************************
* FUNCTION NAME: N_A_B_minus                                                  *
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
                        const int time_line_B[], int time_line_B_size, int Dt);


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
                                int Dt, double tBm, double tAp);
