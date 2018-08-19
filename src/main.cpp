/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: main.cpp                                                         *
*                                                                             *
*******************************************************************************
******************************************************************************/



#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#define NEURONS 183
// total_time_samples 11970
// circ_shifts_num 50

using namespace std;


/******************************************************************************
* FUNCTION NAME: main                                                         *
*                                                                             *
* ARGUMENTS: The total numbers of: Neurons, Time Samples, Circular Shifts.    *
*             And the size of the tile Δt.                                    *
*                                                                             *
* PURPOSE: This main function is for testing purposes only.                   *
*                                                                             *
* RETURNS: 0 on program comletion.                                            *
*                                                                             *
* I/O: Opens/Reads a file containing the spike trains.                        *
*                                                                             *
******************************************************************************/
int main()
{
	vector<int> spike_trains[NEURONS];
	ifstream data;
	data.open("../psm_avalanche", ifstream::in);
	string line;

	int count = 0;
	while (getline(data, line)) {
		for (int n = 0; n < NEURONS; n++) {
			if (line[n] == '1') {
				spike_trains[n].push_back(count);
			}
			count++;
		}
	}

 	int l = 0;
	for (int neur = 0; neur < NEURONS; neur++) {
		for (int fire = 0; fire < spike_trains[neur].size(); fire++) {
			cout<<spike_trains[neur][fire]<<' '<<endl;
      	l++;
		}
		cout<<endl;
	}
  
  	cout<<endl<<l<<endl;
	data.close();
	return 0;
}
