include ../compiler.inc

OBJS	= c_exit.o f_exit.o test_f_exit.o test_f_exit_message.o
FFLAGS	= -O2 -I$(NETCDF_INC)
CFLAGS	= -O2 -I$(NETCDF_INC) -DFORTRAN_UNDERSCORE=_
LD	= $(FC)
MPLD	= $(MPFC)
LDFLAGS	+= -L$(NETCDF_LIB) -lnetcdff
AR      = ar

EXE	= test_f_exit test_f_exit_message
LIB     = libf_exit.a

target: $(EXE) $(LIB)

$(LIB):	$(OBJS)
	$(AR) cru $(LIB) $(OBJS)

test_f_exit_message:	test_f_exit_message.o $(LIB)
			$(LD) $(FFLAGS) $< $(LIB) -o $@

test_f_exit:		test_f_exit.o $(LIB)
			$(LD) $(FFLAGS) $< $(LIB) -o $@

#f_exit.o:	f_exit.f90
#		$(F90) $(FFLAGS) $< -c
#
c_exit.o:	c_exit.c c_exit.h
		$(CC) $(CFLAGS) $< -c

clean:
	@rm -f $(OBJS)

erase:	clean
	@rm -f $(LIB)

install:	$(LIB) $(EXE)
	cp -p $(LIB) $(LIBDIR)

uninstall:	$(LIB)
	rm -f $(LIBDIR)/$(LIB)

.SUFFIXES: .f90 .c .o
.f90.o:
	$(FC) $(FFLAGS) $< -c

.c.o:
	$(CC) $(CFLAGS) $< -c
