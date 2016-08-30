CC=mpiifort
all: def  # run

def: allv
	mpirun -np 4 ./allv

run: mytest
	mpirun -np 3 ./mytest

mytest: mytest.f90
	$(CC) mytest.f90 -o mytest

allv: allv.f90
	$(CC) allv.f90 -o allv

clean:
	rm -rf allv mytest
