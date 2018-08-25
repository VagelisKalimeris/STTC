/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: motif.cpp                                                        *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "motif.hpp"

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
                                                        unsigned int *motifs)
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
