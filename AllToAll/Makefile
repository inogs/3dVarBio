all: run

run: a.out
	mpirun -np 3 ./a.out

a.out: mytest.f90
	mpifort mytest.f90
