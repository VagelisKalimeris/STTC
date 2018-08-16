/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: p_p_null_dist.cpp                                                *
*                                                                             *
*******************************************************************************
******************************************************************************/




/******************************************************************************
* FUNCTION NAME: mean_STTC_dir                                                *
*                                                                             *
* ARGUMENTS: An array containing the shifted per pair STTC results.           *
*                                                                             *
* PURPOSE: Calculates the mean value.                                         *
*                                                                             *
* RETURNS: The mean value.                                                    *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double mean_STTC_dir(double const arr[SHIFTS_NUM]) {
	double sum = 0.0;

	for (int i = 0; i < SHIFTS_NUM; i++) {
		sum += arr[i];
	}
	return sum / (double) SHIFTS_NUM;
}


/******************************************************************************
* FUNCTION NAME: std_STTC_dir                                                 *
*                                                                             *
* ARGUMENTS: An array containing the shifted per pair STTC results.           *
*                                                                             *
* PURPOSE: Calculates the standard deviation.                                 *
*                                                                             *
* RETURNS: The standard deviation.                                            *
*                                                                             *
* I/O: None.                                                                  *
*                                                                             *
******************************************************************************/
double std_STTC_dir(double const arr[SHIFTS_NUM]) {
	double mean = mean_STTC_dir(arr), st_dev = 0.0;

		for (int i = 0; i < SHIFTS_NUM; i++) {
			st_dev += pow((arr[i] - mean), 2) / SHIFTS_NUM;
		}
	return st_dev;
}
