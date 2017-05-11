
include ../compiler.inc

OBJS    = nc-med-level-lib.o
EXE     =
LIB     = libnc-medlevel.a

target: $(LIB)

libnc-medlevel.a: 	 $(OBJS)
			 $(AR) cru $(LIB) $(OBJS)

nc-med-level-lib.o:     nc-med-level-lib.f90
			$(FC) $(FFLAGS) nc-med-level-lib.f90 -c



clean:
	@rm -f $(OBJS)

erase:  clean
	@rm -f $(LIB)

install:        $(LIB)
	cp -p $(LIB) $(LIBDIR)

uninstall:      $(LIB)
	rm -f $(LIBDIR)/$(LIB)

.SUFFIXES: .f90 .c .o
.f90.o:
	$(FC) $(FFLAGS) $< -c
