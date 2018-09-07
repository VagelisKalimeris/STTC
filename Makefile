###############################################################################
#                                                                             #
# PROJECT NAME: STTC Analyses                                                 #
#                                                                             #
# FILE NAME: Makefile                                                         #
#                                                                             #
###############################################################################



CXXFLAGS =  -Wall -Wextra -g -O3 -fopenmp

INCLUDES = INCLUDE/
OBJS = OBJ/
SRCS = SRC/

sttc : ${OBJS}main.o ${OBJS}common.o ${OBJS}cond_null_dist.o ${OBJS}print.o \
			${OBJS}p_p_null_dist.o ${OBJS}triplets_STTC.o ${OBJS}tuplets_STTC.o
	g++ ${CXXFLAGS} -o sttc $^

${OBJS}%.o : ${SRCS}%.cpp ${INCLUDES}%.hpp
	g++ $(CXXFLAGS) -c $< -o $@

${OBJS}%.o : ${SRCS}%.cpp 
	g++ $(CXXFLAGS) -c $< -o $@

clean :
	-rm ${OBJS}*.o sttc 

