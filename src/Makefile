CC = gcc-9 -Os -fopenmp 
#CFLAGS = -L/usr/local/sac/lib -lsac -lsacio -lm
CFLAGS = -L./ -lsac -lsacio -lm
CC_SO = gcc-9 -Os -fopenmp -shared -fPIC
BIN = ../bin

FDTCC: FDTCC.o sacio.o 
	${CC} -o ${BIN}/$@ $^ ${CFLAGS}
FDTCC.so: FCC.o sacio.o
	${CC_SO} -o ${BIN}/$@ $^ ${CFLAGS}
clean:
	rm -f *.o
