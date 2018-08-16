/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: triplets.cpp                                                     *
*                                                                             *
*******************************************************************************
******************************************************************************/



/******************************************************************************
* FUNCTION NAME: STTC_AB_C                                                     *
*                                                                             *
* ARGUMENTS: Three neuron's timelines(references to vectors), and a time      *
*             interval(int).                                                  *
*                                                                             *
* PURPOSE: Estimates the temporal correlation of two neurons, given that a    *
*           third neuron is firing.                                           *
*                                                                             *
* RETURNS: The STTC value(double) of the pair A,B,C.                          *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double STTC_AB_C(const vector<int> &time_line_A, const vector<int> &time_line_B,
                                             const vector<int> &time_line_C, Dt) 
{
  int nApBCA = N_AplusB_CA(time_line_A, time_line_B, time_line_C, Dt);
  int nBmACA =  N_BminusA_CA(time_line_A, time_line_B, time_line_C, Dt);
  double tApt = T_A_plus_tripl(time_line_A, time_line_C, Dt);
  double tBm = T_B_minus(time_line_B, Dt);
  int nA = time_line_A.size(), nB = time_line_A.size();

  return (1/2) * ((((nBmACA / nA) - tBm) / (1 - ((nBmACA / nA) * tBm))) +
                             (((nApBCA / nB) - tApt) / ((nApBCA / nB) * tApt)));
}
