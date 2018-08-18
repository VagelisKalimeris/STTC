/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: tuplets.cpp                                                      *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>
using namespace std;

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
	int i = 0, j = 0, s = 0;
	
	/* all the spikes of A are before or after the tiles of B */
	if(((time_line_A.back() < (time_line_B.front() - Dt)) ||
                               ((time_line_B.back() < time_line_A.front())))) {
	    return P;
	}
	while((i < time_line_A.size()) && (j < time_line_B.size())){
		/* the spike of A is in the tile of spike of B, 
			where tile of B is [tB - Dt, tB] */
		if((time_line_A[i] >= (time_line_B[j] - Dt)) &&
                                          (time_line_A[i] <= time_line_B[j])) {
			s++;
			i++;
		}
		/* the spike of A is before the tile of spike of B */
		else if(time_line_A[i] < (time_line_B[j] - Dt)){
			i++;
		}
		/* the spike of A is after the tile of spike of B */
		else if(time_line_A[i] > time_line_B[j]){
			j++;
		}
	}
	P = s / double(time_line_A.size());
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
	int i = 0, j = 0, s = 0;
	
	/* all spikes of B are before or after tiles of A */
	if((time_line_B.back() < time_line_A.front()) || 
                           ((time_line_A.back() + Dt) < time_line_B.front())) {
	    return P;
	}
	while((i < time_line_B.size()) && (j < time_line_A.size())){
		/* spike of B is within tile of spike of A [tA, tA + Dt] */
		if((time_line_B[i] >= time_line_A[j]) &&
                                   (time_line_B[i] <= (time_line_B[j] + Dt))) {
			s++;
			i++;
		}
		/* spike of B is before tile of spike of A [tA, tA + Dt] */
		else if(time_line_B[i] < time_line_A[j]) {
			i++;
		}
		/* spike of B is after tile of spike of A [tA, tA + Dt] */
		else if(time_line_B[i] > (time_line_A[j] + Dt)) {
			j++;
		}
	}
	
	P = s / double(time_line_B.size());
	
	return P;
}


/******************************************************************************
* FUNCTION NAME: STTC_A_B                                                     *
*                                                                             *
* ARGUMENTS: Two neuron's timelines(references to vectors), and a time        *
*             interval(int).                                                  *
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
																		int Dt)
{
  double pABm = P_A_B_minus(time_line_A, time_line_B, Dt);
  double tBm = T_B_minus(time_line_B, Dt);
  double pBAp = P_B_A_plus(time_line_B, time_line_A, Dt);
  double tAp = T_A_plus(time_line_A, Dt);

  return (1/2) * (((pABm - tBm) / (1 - (pABm * tBm))) +
                                                ((pBAp - tAp) / (pBAp * tAp)));
}
