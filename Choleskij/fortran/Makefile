FORT = ifort -mkl=sequential #gfortran
FLAGS = -O3 -heap-arrays 2048 #-fbounds-check -ffree-line-length-none -lblas -llapack

all: header.o subroutines.o main.o
	$(FORT) $(FLAGS) header.o functions.o main.o -o main
	
header.o: header.f90
	$(FORT) $(FLAGS) -c header.f90 -o header.o

subroutines.o: functions.f90
	$(FORT) $(FLAGS) -c functions.f90 -o functions.o

main.o: main.f90
	$(FORT) $(FLAGS) -c main.f90 -o main.o
