/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: motif.hpp                                                        *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/******************************************************************************
* FUNCTION NAME: print_all_spikes                                             *
*                                                                             *
* ARGUMENTS: A neuron's timeline(reference to a vector), the total number of  *
*             neurons(int).                                                   *   
*                                                                             *  
* PURPOSE: Prints all the spikes of each neuron of the dataset, as well as    *
*           the total number of spikes.                                       *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: See PURPOSE.                                                           *
*                                                                             *
******************************************************************************/
void print_all_spikes(const vector<int> spike_trains[], 
                        const int total_neurons, const vector<int> &astrocytes, 
                        ofstream &info, const string output, 
                        const string shifts, const string Dt);


/******************************************************************************
* FUNCTION NAME: print_motifs                                                 *
*                                                                             *
* ARGUMENTS: The totals of each motif category of triplets (reference to      *
*             an array), the totals of each motif category of significant     *
*              triplets (reference to an array).                              *
*                                                                             *
* PURPOSE: Prints each motif category of triplets.                            *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void print_motifs(const int *triplets, const int *significants, 
                                    ofstream &info, const string output, 
                                    const string shifts, const string Dt);
