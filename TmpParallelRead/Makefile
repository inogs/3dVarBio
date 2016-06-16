CC=mpif90
FFLAGS=-O2
LFLAGS=-L$(PNETCDF_LIB) -lpnetcdf
INCLUDE=-I$(PNETCDF_INC)
EXE = MyRead
all: $(EXE)

MyRead: MyRead.f90
	$(CC) $(INCLUDE) MyRead.f90 $(FFLAGS) $(LFLAGS) -o $@

run: all
	@echo ""
	mpirun -np 4 ./MyRead
	@echo ""

clean:
	rm -rf MyRead