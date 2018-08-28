###############################################################################
#                                                                             #
# PROJECT NAME: STTC Analyses                                                 #
#                                                                             #
# FILE NAME: Makefile                                                         #
#                                                                             #
###############################################################################



COMPILE_FLAGS =  -Wall -Wextra -g -O3 -fopenmp


sttc : main.o common.o tuplets_STTC.o triplets_STTC.o p_p_null_dist.o  \
                                                        cond_null_dist.o motif.o
	g++ ${COMPILE_FLAGS} -o sttc main.o common.o tuplets_STTC.o triplets_STTC.o \
                                    p_p_null_dist.o cond_null_dist.o motif.o

main.o : main.cpp common.hpp tuplets_STTC.hpp triplets_STTC.hpp p_p_null_dist.hpp \
                                                    cond_null_dist.hpp motif.hpp
	g++ -c ${COMPILE_FLAGS} main.cpp

common.o : common.cpp common.hpp
	g++ -c ${COMPILE_FLAGS} common.cpp

tuplets.o : tuplets_STTC.cpp tuplets_STTC.hpp
	g++ -c ${COMPILE_FLAGS} tuplets_STTC.cpp

triplets.o : triplets_STTC.cpp triplets_STTC.hpp
	g++ -c ${COMPILE_FLAGS} triplets_STTC.cpp

p_p_null_dist.o : p_p_null_dist.cpp p_p_null_dist.hpp
	g++ -c ${COMPILE_FLAGS} p_p_null_dist.cpp

cond_null_dist.o : cond_null_dist.cpp cond_null_dist.hpp
	g++ -c ${COMPILE_FLAGS} cond_null_dist.cpp

motif.o : motif.cpp motif.hpp
	g++ -c ${COMPILE_FLAGS} motif.cpp


clean :
	-rm *.o sttc output *.csv *.txt
