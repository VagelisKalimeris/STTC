#pragma once

#include <functional>
#include "Times.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <string>
#include <array>

#define Function std::bind

using namespace std;

static unsigned int get_unit_time_measure(vector<Times> &v) {
	unsigned int unit_time = 2;
	if (any_of(v.begin(), v.end(), [&v](Times &x) { return x.tim > 1.0; })) {
		unit_time = 0;
	}
	else if (any_of(v.begin(), v.end(), [&v](Times &x) { return x.tim*1e+03 > 1.0; })) {
		unit_time = 1;
	}
	return unit_time;
}

static string operator*(string s,int x) {
	string res;
	for (int i = 0; i < x; ++i)
		res += s;
	return res;
}

namespace benchmark{

	using std::chrono::steady_clock;
	using std::vector;

	template<class F>
	vector<Times> benchmark(F funct, Times time){
		vector<Times> tim;
		steady_clock sc;
		(void)sc; //  an to afaireso bgazei warning
		for (unsigned int i = 0; i < time.tim; ++i){
			auto start = sc.now();
			funct();
			auto end = sc.now();
			auto t = static_cast<chrono::duration<double>>(end - start).count();
			tim.push_back(t);
		}
		return tim;
	}
	
	void print(vector<Times> v, unsigned int number_of_precission = 4) {
		auto mnmx = minmax_element(v.begin(), v.end(), [&v](Times &x, Times &y) { return x.tim < y.tim; });
		auto mean_t = std::accumulate(v.begin(), v.end(), 0.0, [&v](double &s, Times &x) { return s + x.tim; }) / v.size();
		auto mn = mnmx.first->tim;
		auto mx = mnmx.second->tim;

		std::cout << std::setprecision(number_of_precission) << std::fixed;

		string empty_space, unit_time,header_space;

		switch (get_unit_time_measure(v)) {
		case 0: {  //seconds
			unit_time = "seconds:	";
			header_space = "        	";
			break;
		}
		case 1: {  //milliseconds
			unit_time = "milliseconds:	";
			header_space = "             	";
			mn *= 1000;
			mean_t *= 1000;
			mx *= 1000;
			break;
		}
		default: {  //microseconds
			unit_time = "microseconds:	";
			header_space = "             	";
			mn *= 1e+06;
			mean_t *= 1e+06;
			mx*= 1e+06;
			break;
		}
		}

		auto length_mn = static_cast<unsigned int>(to_string(static_cast<unsigned int>(mn)).size());
		auto length_mx = static_cast<unsigned int>(to_string(static_cast<unsigned int>(mx)).size());
		auto length_mean_t = static_cast<unsigned int>(to_string(static_cast<unsigned int>(mean_t)).size());

		std::array<unsigned int, 3> all_the_values = { length_mn, length_mx, length_mean_t };

		auto num_int_digits = *max_element(all_the_values.begin(), all_the_values.end());
		auto max_length = num_int_digits + number_of_precission + 1; // +1 for the dot

		for (unsigned int i = 0; i < max_length; ++i)
			empty_space += " ";		
		
		std::cout << header_space << "min" << empty_space + string(" ")*(num_int_digits + 2) << "mean"
				  << empty_space + string(" ")*(num_int_digits + 1) << "max\n";
		std::cout << unit_time << mn << empty_space << mean_t << empty_space << mx << "\n\n";
	}

}