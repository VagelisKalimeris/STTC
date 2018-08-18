/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.hpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <cmath>
#include<vector>
using namespace std;

/******************************************************************************
* FUNCTION NAME: T_A_plus                                                     *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), and a time           *
*             interval(int).                                                  *
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
	const int Dt_1 = Dt + 1;
	if (Dt == 0) {
		// if Dt is zero then return mean
		T = time_line_A.size() / double(total_time_samples);
	}
	else {
		int s = 0, last_spike = 0;
	    for (auto &spike : time_line_A){ // for each spike
	       	for (int j = 0; j < Dt_1; ++j){ // check all the next
                if((spike + j <= total_time_samples) && (spike+j > last_spike)){
					++s;
	          }
	       	}
		   	last_spike = spike + Dt; //  keep the last spike
	    }
	    T = s / double(total_time_samples);
	}
	return T;
}


/******************************************************************************
* FUNCTION NAME: T_B_minus                                                    *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), and a time           *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Calculates the sum of time tiles before a neuron's firing, divided *
*           by the total time.                                                *
*                                                                             *
* RETURNS: The total time(double).                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double T_B_minus(const vector<int> &time_line_B, int total_time_samples, 
                                                                        int Dt)
{
	double T = 0.0;
	const int Dt_1 = Dt + 1;
	if (Dt == 0) {
		// if Dt is zero then return mean
		T = time_line_B.size() / double(total_time_samples);
	}
	else {
        int s = 0, last_spike = -1; // -1 counts the case: first spike-j = zero
	    for (auto &spike : time_line_B){ // for each spike 
	       for (int j = 0; j < Dt_1; ++j){ // check all the previous spikes
	          if((spike - j) > last_spike){
				  ++s;
	          }
	       }
		   last_spike = spike; // keep the first spike
	    }
	    T = s / double(total_time_samples);
	}
	return T;
}


/******************************************************************************
* FUNCTION NAME: sign_thresh_A_B                                              *
*                                                                             *
* ARGUMENTS: The mean(double) and the standard(double) deviations of the      *
*             circ. shifted spike trains.                                     *
*                                                                             *
* RETURNS: The significant threshhold.                                        *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double sign_thresh_A_B(double mean, double st_dev)
{
	return mean + (3 * st_dev);
}


/******************************************************************************
* FUNCTION NAME: circular_shift                                               *
*                                                                             *
* ARGUMENTS: A vector representing the firings of a neuron, and a random      *
*             number between 0 and the maximum number of time events.         *
*                                                                             *
* PURPOSE: Shifts forward each firing of a neuron by a random number.         *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void circular_shift(vector<int> &time_line, int random) {
	int max = time_line.back();
	vector<int>::iterator front_it = time_line.begin();

	for (int i = 0; i < time_line.size(); i++) {
		int temp = time_line[i] + random;

		if ((temp) < max) {
			time_line[i] = temp;
		}
		else {
			time_line.erase(time_line.begin() + i);
			temp = temp - max - 1;
			time_line.insert(front_it, temp);
			++front_it;
		}
	}
}
