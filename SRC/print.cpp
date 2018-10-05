/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: motif.cpp                                                        *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include "../INCLUDE/print.hpp"

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
                        const string shifts, const string Dt)
{
    int ttl_firings = 0, max = 0, min = 100000;
    const int astro_size = astrocytes.size();
    const int neur_clean = total_neurons - astro_size;
    int astro = 0, astrocyte = 0;
    vector<int> time_lines;
    if (astro_size) {
        astrocyte = astrocytes[0];
    }
    
    ofstream spikes;
    spikes.open(("RESULTS/" + output + "_" + shifts + "-shifts_" + Dt + 
                                            "-dt_neurons_spikes.txt").c_str());
    if (!spikes.is_open()) {
        cout<<"Error opening results neurons spikes file!"<<endl;
        return;
    }
    spikes<<"\nThe data structure: "<<endl;
    for (int neur = 0; neur < total_neurons; neur++) {
        int pos, time_line_size;
        if (astro_size && neur == astrocyte) {
            pos = neur_clean + astro;
            time_line_size = spike_trains[pos].size();
            spikes<<"No "<<neur + 1<<" astrocyte neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            astrocyte = astrocytes[(++astro) % astro_size];
        }
        else {
            pos = neur - astro;
            time_line_size = spike_trains[pos].size();
            spikes<<"No "<<neur + 1<<" neuron's spikes ("
                                                    <<time_line_size<<"):\n";
            ttl_firings += time_line_size;
            if (time_line_size > max) {
                max = time_line_size;
            }
            else if (time_line_size < min) {
                min = time_line_size;
            }
            time_lines.push_back(time_line_size);
        }
        for (int fire = 0; fire < time_line_size; fire++) {
            spikes<<spike_trains[pos][fire] + 1<<' ';
        }
        spikes<<endl<<endl;
    }
    spikes.close();
    
    double median;
    sort(time_lines.begin(), (time_lines.begin() + neur_clean));
    if (neur_clean % 2) {
        median = time_lines[neur_clean/2]/1.0;
    }
    else {
        median = (time_lines[neur_clean/2 - 1] + time_lines[neur_clean/2])/2.0;
    }
    
    info<<"\nNeurons' analytics:"<<endl;
    info<<"Total number of neurons: "<<total_neurons<<endl;
    info<<"Total number of astrocytes: "<<astro_size<<endl;
    info<<"\nSpikes info excluding astrocytes:"<<endl;
    info<<"Total number of spikes in all neurons: "<<ttl_firings<<endl;
    info<<"Maximum number of spikes in a neuron: "<<max<<endl;
    info<<"Minimum number of spikes in a neuron: "<<min<<endl;
    info<<"Average number of spikes: "<<ttl_firings / double(neur_clean)<<endl;
    info<<"Median number of spikes: "<<median<<endl;
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
void print_motifs(const int *triplets, const int *significants, 
                                    ofstream &info, const string output, 
                                    const string shifts, const string Dt)
{
    int ttl_triplets = 0, ttl_sgnfcnt_triplets = 0;
    
    ofstream motifs;
    motifs.open(("RESULTS/" + output + "_" + shifts + "-shifts_" + Dt + 
                                                    "-dt_motifs.txt").c_str());
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
