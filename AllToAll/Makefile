all: def  # run

def: allv
	mpirun -np 3 ./allv

run: mytest
	mpirun -np 3 ./mytest

mytest: mytest.f90
	mpifort mytest.f90 -o mytest

allv: allv.f90
	mpifort allv.f90 -o allv

clean:
	rm -rf allv mytest
