/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: common.cpp                                                       *
*                                                                             *
*******************************************************************************
******************************************************************************/



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
double T_A_plus(const vector<int> &time_line_A, int Dt){
	double T = 0.0;
	const int Dt_1 = Dt + 1;
	if (Dt == 0) {
		// if Dt is zero then return mean
		T = time_line_A.size() / double(TIME_STAMPS);
	}
	else {
	    int s = 0, last_spike = 0;
	    for (auto &spike : time_line_A){ // for each spike
	       	for (int j = 0; j < Dt_1; ++j){ // check all the next
			if((spike + j <= TIME_STAMPS) && (spike+j > last_spike))
					++s;
	       	}
		last_spike = spike + Dt; //  keep the last spike
	    }
	    T = s / double(TIME_STAMPS);
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
double T_B_minus(const vector<int> &time_line_B, int Dt){
	double T = 0.0;
	const int Dt_1 = Dt + 1;
	if (Dt == 0) {
		// if Dt is zero then return mean
		T = time_line_B.size() / double(TIME_STAMPS);
	}
	else {
	    int s = 0, last_spike = -1; // -1 counts the case: 1st spike-j=0
	    for (auto &spike : time_line_B){ // for each spike 
	       for (int j = 0; j < Dt_1; ++j){ // check previous spikes
	          if((spike - j) > last_spike){
				  ++s;
	          }
	       }
	       last_spike = spike; // keep the first spike
	    }
	    T = s / double(TIME_STAMPS);
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
