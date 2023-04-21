.PHONY: all clean

CC = g++ # compiler to use
CFLAGS = -c -O3 -I. -IRoutines -I.. -I../Routines #-Wall -Werror -Wextra -pedantic# flags to use at the compliation
LINKERFLAG = -lm -lfftw3

all: prog

prog: obj/mode_solver.o obj/mode_lib.o
	${CC} -o solver $(LIB_BINS) obj/mode_solver.o ${LINKERFLAG}
	${CC} -fPIC -shared -o ../mode_lib.so $(CORE_BINS) obj/mode_lib.o ${LINKERFLAG}

obj:
	@echo "Making obj directory..."
	mkdir obj

obj/mode_solver.o: obj mode_solver.cpp
	@echo "Compiling mode_solver.cpp..."
	${CC} ${CFLAGS} -o obj/mode_solver.o mode_solver.cpp

# eigen_binary_tmp.o: eigen_binary_tmp.cpp
# 	@echo "Compiling eigen_binary_tmp.cpp..."
# 	${CC} ${CFLAGS} -o eigen_binary_tmp.o eigen_binary_tmp.cpp

obj/mode_lib.o: obj mode_lib.cpp
	@echo "Compiling mode_lib.cpp..."
	${CC} ${CFLAGS} -fPIC -o obj/mode_lib.o mode_lib.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -vf solver
	rm -vf ../mode_lib.so