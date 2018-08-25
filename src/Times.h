#pragma once

#include <iostream>

template<class T>
struct times{
	T tim;
	times(T tim=0){
		this->tim = tim;
	}
	static void print(){
		std::cout << "              " << "time\n";
		std::cout << "seconds:      " << times::tim << std::endl;
		std::cout << "milliseconds: " << times::tim * 1000 << std::endl;
		std::cout << "microseconds: " << times::tim * 1000000 << std::endl;
	}
};

using Times = times<double>;
