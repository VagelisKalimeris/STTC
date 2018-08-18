#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#define NEURONS 183
// total_time_samples 11970
// circ_shifts_num 50

using namespace std;

// This main function is for testing purposes only.
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
