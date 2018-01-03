  #---------------------------------------------------------------------------
  #                                                                          #
  #    Copyright 2018 Anna Teruzzi, OGS, Trieste                         #
  #                                                                          #
  #    This file is part of 3DVarBio.
  #                                             #
  #                                                                          #
  #    3DVarBio is  free software: you can redistribute it and/or modify.     #
  #    it under the terms of the GNU General Public License as published by  #
  #    the Free Software Foundation, either version 3 of the License, or     #
  #    (at your option) any later version.                                   #
  #                                                                          #
  #    3DVarBio is  distributed in the hope that it will be useful,           #
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
  #    GNU General Public License for more details.                          #
  #                                                                          #
  #    You should have received a copy of the GNU General Public License     #
  #    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       #
  #                                                                          #
  #---------------------------------------------------------------------------

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
