#ifndef __c_exit__h__
#define __c_exit__h__

#define XCAT(x,y) x##y

#define CAT(x,y) XCAT(x,y)

#define C_EXIT CAT(c_exit,FORTRAN_UNDERSCORE)

void C_EXIT(int * exit_code);

#endif
