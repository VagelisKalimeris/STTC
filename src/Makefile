COMPILE_FLAGS =  -Wall -Wextra -g


sttc : main.o common.o tuplets.o triplets.o pp_null_dist.o cond_null_dist.o
	gcc ${COMPILE_FLAGS} main.o common.o tuplets.o triplets.o pp_null_dist.o -o sttc


main.o : main.c
	gcc -c ${COMPILE_FLAGS} main.c

common.o : common.c common.h
	gcc -c ${COMPILE_FLAGS} common.c

tuplets.o : tuplets_STTC.c tuplets_STTC.h common.h
	gcc -c ${COMPILE_FLAGS} tuplets_STTC.c

triplets.o : triplets_STTC.c triplets_STTC.h common.h
	gcc -c ${COMPILE_FLAGS} triplets_STTC.c

pp_null_dist.o : pp_null_dist.c pp_null_dist.h tuplets_STTC.h common.h
	gcc -c ${COMPILE_FLAGS} pp_null_dist.c

cond_null_dist.o : cond_null_dist.c triplets_STTC.h common.h
	gcc -c ${COMPILE_FLAGS} cond_null_dist.c



clean :
	-rm *.o sttc
