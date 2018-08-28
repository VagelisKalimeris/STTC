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
static void categorization(const bool CtoA, const bool CtoB, const bool AtoB, 
                                    int *motifs) __attribute__((always_inline));

static void categorization(const bool CtoA, const bool CtoB, const bool AtoB, 
                                                                int *motifs)
{
    if (CtoA) {
        if (CtoB) {
            if (AtoB) {
                ++motifs[7];
            }
            else {
                ++motifs[6];
            }
        }
        else {
            if (AtoB) {
                ++motifs[5];
            }
            else {
                ++motifs[4];
            }
        }
    }
    else {
        if (CtoB) {
            if (AtoB) {
                ++motifs[3];
            }
            else {
                ++motifs[2];
            }
        }
        else {
            if (AtoB) {
                ++motifs[1];
            }
            else {
                ++motifs[0];
            }
        }
    }
}


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
void print_motifs(const int *triplets, const int *significants);
