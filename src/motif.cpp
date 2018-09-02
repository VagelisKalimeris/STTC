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
void print_motifs(const int *triplets, const int *significants, 
                                        ofstream &info, const string output)
{
    int ttl_triplets = 0, ttl_sgnfcnt_triplets = 0;
    
    ofstream motifs;
    motifs.open(("RESULTS/" + output + "_motifs.csv").c_str());
    if (!motifs.is_open()) {
        cout<<"Error opening results motifs file!"<<endl;
        return;
    }
    motifs<<"Motif,Triplets,Significants\n";
    info<<"\n+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    info<<"| Motif |"<<setfill(' ')<<setw(16)<<"Triplets |"
                                            <<setw(16)<<"Significant |"<<endl;
    info<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    for (int m = 0; m < 8; m++) {
        int triplet = triplets[m];
        int sgnfcnt = significants[m];
        motifs<<m<<','<<triplet<<','<<sgnfcnt<<endl;
        ttl_triplets += triplet;
        ttl_sgnfcnt_triplets += sgnfcnt;
        info<<"|   "<<m<<"   | "<<setfill(' ')<<setw(13)<<triplet<<" | "
                                                <<setw(13)<<sgnfcnt<<" |"<<endl;
    }
    motifs.close();
    info<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
    info<<"| Total | "<<setfill(' ')<<setw(13)<<ttl_triplets<<" | "
                                <<setw(13)<<ttl_sgnfcnt_triplets<<" |"<<endl;
    info<<"+"<<setfill('-')<<setw(8)<<'+'<<setw(16)<<'+'<<setw(16)<<'+'<<endl;
}
