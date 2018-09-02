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

using namespace std;

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
                                        ofstream &info, const string output);
