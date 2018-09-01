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
