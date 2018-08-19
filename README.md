# STTC ANALYSES

# GENERAL DESCRIPTION

1. Directional STTC

To incorporate the temporal order of the firing events of two neurons, we extended the STTC as follows: the fraction of the number of the firing events of A which fall within Δt before each firing event of B by the number of firing events of A is computed (𝑃𝐴𝐵−). The firing events that have this property increase the positive correlation (Fig. 1). Similarly, the fraction of the number of firing events of B that fall within Δt after each firing event of A by the number of firing events of B is estimated(𝑃𝐵𝐴+). In this way we estimate the correlation between spike trains for the spikes of B that follows spikes of A and the spikes of A that proceeds spikes of B.


2. Conditional STTC of two neurons, given the firing of a third one

To estimate the temporal correlation of two neurons, given that a third neuron is firing, we defined the
conditional STTC as follows: we identify the firing events of A which follow within an interval of Δt of a
firing event of C (Fig. 2). This new sequence of firing events of A forms the “reduced” spike train A. The
number of firing events of B that falls within the tiles Δt after the firing events of the reduced spike train
𝑪𝑨
A (𝑵𝑨+𝑩) and the number of firing events of the reduced spike train A that falls within the tiles Δt after
𝑪𝑨
the firing events of spike train B is estimated (𝑵𝑩−𝑨). The spikes that have this property increase the
positive correlation between A and B given the spike train C.


3. Per pair threshold null distribution test for directional STTC

For a given pair (A,B) we circular shift the spike train of the neuron A 50 times and estimate the directional STTC between the circular shifted spike trains A and the spike train B . Based on these 50 values we estimate the mean value (mean_sttc_direct) and the standard deviation (std_sttc_direct). The significant threshold corresponds to:
Significant threshold = mean_sttc_direct + 3 std_sttc_direct
If the directional STTC value of the given pair (A, B) is greater than the threshold, then we consider the edge between this pair as significant. In this approach we set one significant threshold value per pair, based on the 50 produced values and compare each pair (A,B) with the corresponding significant threshold value.


4. Null distribution for conditional STTC

For a given triplet A->B|C we circular shift the spike train of the neuron C 50 times and estimate the conditional STTC between A->B and the circular shifted spike trains C. Based on these 50 values we estimate the mean value (mean_sttc_cond) and the standard deviation (std_sttc_cond). The significant threshold corresponds to:
Significant threshold = mean_sttc_cond + 3 std_sttc_cond
A second limitation in order to consider a triplet as significant is that the number of firing events of A that follows ΔΤ after each spike of C (‘reduced A’).
If the conditional STTC of the given triplet is greater than the significant threshold and the number of firing events of ‘reduced A’ is greater than 5, then we consider this triplet as significant.
In this approach we set one significant threshold value per triplet, based on the 50 produced values and compare each triplet (A-> B)|C with the corresponding significant threshold value.



# Ερωτήσεις για Κα Παπαδοπούλη:

1. Τα double_STTC και triple_STTC θα υπολογίζονται πάντα ταυτόχρονα ή και ξεχωριστά? 
  Στην πρώτη περίπτωση θα μπορούσαμε να γλυτώνουμε τον υπολογισμό του double_STTC δεύτερη φορα.

2. Ποιές από τις εξής μεταβλητές θα είναι #define μεσα στο αρχείο main, και ποιές θα δίνονται από το χρήστη 
  ως command line arguments? : number of total neurons, length of spike trains(total time samples), number of circular 
  shifts, size of tile Dt. 

3. Υπάρχει περίπτωση το Dt να είναι 0?

4. Standard deviation διαιρούμε με Ν η Ν-1?

5. C++11?

