CPP = g++
FLAGS = -Wall -O2
OBJ = main.o matrix_func.o

all: ${OBJ}
	${CPP} ${FLAGS} -o main ${OBJ}

main.o: main.cpp matrix_func.h
	${CPP} ${FLAGS} -c main.cpp

matrix_ops.o: matrix_func.cpp matrix_func.h
	${CPP} ${FLAGS} -c matrix_func.cpp
