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

sttc : ${OBJS}*.o 
	g++ ${CXXFLAGS} -o sttc $< 

${OBJS}%.o : ${SRCS}%.cpp ${INCLUDES}%.hpp
	g++ $(CXXFLAGS) -c $< -o $@

${OBJS}%.o : ${SRCS}%.cpp 
	g++ $(CXXFLAGS) -c $< -o $@

clean :
	-rm ${OBJS}*.o sttc 

