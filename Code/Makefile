.PHONY: all clean

CC = g++ # compiler to use
CFLAGS = -c -O3 -I. -IRoutines #-Wall -Werror -Wextra -pedantic# flags to use at the compliation
LINKERFLAG = -lm

all: prog

prog: obj/bose_shooter.o obj/shooter_lib.o
	${CC} -o shooter $(LIB_BINS) obj/bose_shooter.o ${LINKERFLAG}
	${CC} -fPIC -shared -o shooter_lib.so $(CORE_BINS) obj/shooter_lib.o ${LINKERFLAG}

obj:
	@echo "Making obj directory..."
	mkdir obj

obj/bose_shooter.o: obj bose_shooter.cpp
	@echo "Compiling bose_shooter.cpp..."
	${CC} ${CFLAGS} -o obj/bose_shooter.o bose_shooter.cpp

obj/shooter_lib.o: obj shooter_lib.cpp
	@echo "Compiling shooter_lib.cpp..."
	${CC} ${CFLAGS} -fPIC -o obj/shooter_lib.o shooter_lib.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -rvf Data
	rm -vf shooter
	rm -vf shooter_lib.so