CC=mpiifort #f90
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

test: MyReadTestfile WriteBuff
	@echo ""
	mpirun -np 4 ./WriteBuff
	mpirun -np 4 ./MyReadTestfile
	@echo ""

grid: ProcGrid WriteBuff
	@echo ""
	mpirun -np 4 ./WriteBuff
	mpirun -np 4 ./ProcGrid
	@echo ""

MyReadTestfile: MyReadTestfile.f90
	$(CC) $(INCLUDE) $^ $(FFLAGS) $(LFLAGS) -o $@

ProcGrid: ProcGrid.f90
	$(CC) $(INCLUDE) $^ $(FFLAGS) $(LFLAGS) -o $@

WriteBuff: pnetcdf-write-bufferedf.f90
	$(CC) $(INCLUDE) $^ $(FFLAGS) $(LFLAGS) -o $@

clean:
	rm -rf MyRead WriteBuff MyReadTestfile ProcGrid