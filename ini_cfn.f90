subroutine ini_cfn
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
  !                                                                          !
  !    This file is part of OceanVar.                                          !
  !                                                                          !
  !    OceanVar is free software: you can redistribute it and/or modify.     !
  !    it under the terms of the GNU General Public License as published by  !
  !    the Free Software Foundation, either version 3 of the License, or     !
  !    (at your option) any later version.                                   !
  !                                                                          !
  !    OceanVar is distributed in the hope that it will be useful,           !
  !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
  !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
  !    GNU General Public License for more details.                          !
  !                                                                          !
  !    You should have received a copy of the GNU General Public License     !
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
  !                                                                          !
  !---------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Initialise the minimisation                                          !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str

#ifdef _USE_MPI
  use mpi
  use mpi_str
#endif
  
  implicit none
  
  INTEGER(i4)  :: i

#ifdef _USE_MPI
  INTEGER(i4)  :: ierr
#endif
  
  ctl%task = 'START'
  
  ctl%iprint = 0 ! if 1 prints file iterate.dat

  ctl%factr=1.0d+7
  
  !  ctl%pgtol=1.0d-2
  
  
  if( drv%ktr.eq.1 .or. drv%ratio(drv%ktr).ne.1.0 ) then
     
     ! ---
     ! Allocate memory for optimization arrays
     
     ctl%n = nSurfaceWaterPoints * ros%neof
#ifdef _USE_MPI
     call MPI_Allreduce(ctl%n, ctl%n_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)

     if (MyRank .eq. 0) write(drv%dia,*) 'Size of the control vector: ',ctl%n_global

#else
     write(drv%dia,*) 'Size of the control vector: ',ctl%n
#endif

     ALLOCATE( ctl%nbd(ctl%n)) ; ctl%nbd = huge(ctl%nbd(1))
     ALLOCATE(ctl%iwa(3*ctl%n)); ctl%iwa = huge(ctl%iwa(1))
     ALLOCATE( ctl%x_c(ctl%n)) ; ctl%x_c = huge(ctl%x_c(1))
     ALLOCATE( ctl%g_c(ctl%n)) ; ctl%g_c = huge(ctl%g_c(1))
     ALLOCATE( ctl%l_c(ctl%n)) ; ctl%l_c = huge(ctl%l_c(1))
     ALLOCATE( ctl%u_c(ctl%n)) ; ctl%u_c = huge(ctl%u_c(1))
     ALLOCATE( ctl%wa(8*ctl%m)); ctl%wa  = huge(ctl%wa(1))
     ALLOCATE( ctl%sg(ctl%m))  ; ctl%sg  = huge(ctl%sg(1))
     ALLOCATE( ctl%sgo(ctl%m)) ; ctl%sgo = huge(ctl%sgo(1))
     ALLOCATE( ctl%yg(ctl%m))  ; ctl%yg  = huge(ctl%yg(1))
     ALLOCATE( ctl%ygo(ctl%m)) ; ctl%ygo = huge(ctl%ygo(1))
     
     
     
     ALLOCATE( ctl%ws(ctl%n,ctl%m)) ; ctl%ws = huge(ctl%ws(1,1)) ;
     ALLOCATE( ctl%wy(ctl%n,ctl%m)) ; ctl%wy = huge(ctl%wy(1,1))
     ALLOCATE( ctl%sy(ctl%m,ctl%m)) ; ctl%sy = huge(ctl%sy(1,1))
     ALLOCATE( ctl%ss(ctl%m,ctl%m)) ; ctl%ss = huge(ctl%ss(1,1))
     ALLOCATE( ctl%yy(ctl%m,ctl%m)) ; ctl%yy = huge(ctl%yy(1,1))
     
     ALLOCATE( ctl%wt(ctl%m,ctl%m))      ; ctl%wt  = huge(ctl%wt(1,1))
     ALLOCATE( ctl%wn(2*ctl%m,2*ctl%m))  ; ctl%wn  = huge(ctl%wn(1,1))
     ALLOCATE( ctl%snd(2*ctl%m,2*ctl%m)) ; ctl%snd = huge(ctl%snd(1,1))
     
     ALLOCATE( ctl%z_c(ctl%n)) ; ctl%z_c = huge(ctl%z_c(1))
     ALLOCATE( ctl%r_c(ctl%n)) ; ctl%r_c = huge(ctl%r_c(1))
     ALLOCATE( ctl%d_c(ctl%n)) ; ctl%d_c = huge(ctl%d_c(1))
     ALLOCATE( ctl%t_c(ctl%n)) ; ctl%t_c = huge(ctl%t_c(1))
     
     
     do i=1,ctl%n
        ctl%nbd(i)=0
        ctl%l_c(i)=-1.0d3
        ctl%u_c(i)= 1.0d3
        ctl%x_c(i)= 0.0d0
     enddo
     
  endif
  ctl%f_c = 0.0
  ctl%csave = 'gnignigni'
  ctl%lsave = .true.
  ctl%isave = 0
  ctl%dsave = 0.0
  
end subroutine ini_cfn
