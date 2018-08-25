/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: motif.hpp                                                        *
*                                                                             *
*******************************************************************************
******************************************************************************/


using namespace std;

/******************************************************************************
* FUNCTION NAME: categorization                                               *
*                                                                             *
* ARGUMENTS: .                                                                *
*                                                                             *
* PURPOSE: Categorizes the triplet of neurons to one motif category.          *
*                                                                             *
* RETURNS: None.                                                              *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
void categorization(const bool CtoA, const bool CtoB, const bool AtoB, 
                                                        unsigned int *motifs);
