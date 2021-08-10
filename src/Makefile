#mac
CC = gcc-9 -Os -fopenmp 
#linux
#CC = gcc -Os -fopenmp
CFLAGS = -L/usr/local/sac/lib -lsac -lsacio -lm
BIN = ../bin

FDTCC: FDTCC.o sacio.o 
	${CC} -o ${BIN}/$@ $^ ${CFLAGS}
clean:
	rm -f *.o
