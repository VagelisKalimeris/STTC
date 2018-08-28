/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: motif.cpp                                                        *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/motif.hpp"

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
void print_motifs(const int *triplets, const int *significants)
{
    int ttl_triplets = 0, ttl_sgnfcnt_triplets = 0;
    
    cout<<"\n+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    cout<<"| Motif |"<<setfill(' ')<<setw(16)<<"Triplets |"
                                            <<setw(16)<<"Significant |"<<endl;
    cout<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    for (int m = 0; m < 8; m++) {
        int triplet = triplets[m];
        int sgnfcnt = significants[m];
        ttl_triplets += triplet;
        ttl_sgnfcnt_triplets += sgnfcnt;
        cout<<"|   "<<m<<"   | "<<setfill(' ')<<setw(13)<<triplet<<" | "
                                                <<setw(13)<<sgnfcnt<<" |"<<endl;
    }
    cout<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    cout<<"| Total | "<<setfill(' ')<<setw(13)<<ttl_triplets<<" | "
                                <<setw(13)<<ttl_sgnfcnt_triplets<<" |"<<endl;
    cout<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
}
