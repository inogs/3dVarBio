all: run

run: a.out
	mpirun -np 2 ./a.out

a.out: mytest.f90
	mpifort mytest.f90
