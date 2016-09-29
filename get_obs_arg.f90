subroutine get_obs_arg
  
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
  ! Load ARGO observations                                                !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  
  implicit none
  
  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, kk, i
  
  arg%no = 0
  arg%nc = 0

  
!  open(511,file='arg_mis.dat',form='formatted')
  open(511,file='arg_datnew.dat') !,form='formatted')
  
  ! ---
  ! Allocate memory for observations   
  read(511,'(I5)') arg%no
  write(drv%dia,*)'Number of ARGO observations: ',arg%no
  
  if(arg%no.eq.0)then
     close(511)
     return
  endif
  ALLOCATE ( arg%ino(arg%no), arg%flg(arg%no), arg%flc(arg%no), arg%par(arg%no))
  ALLOCATE ( arg%lon(arg%no), arg%lat(arg%no), arg%dpt(arg%no), arg%tim(arg%no))
  ALLOCATE ( arg%val(arg%no), arg%bac(arg%no), arg%inc(arg%no))
  ALLOCATE ( arg%bia(arg%no), arg%err(arg%no))
  ALLOCATE ( arg%res(arg%no), arg%b_a(arg%no))
  ALLOCATE ( arg%ib(arg%no), arg%jb(arg%no), arg%kb(arg%no))
  ALLOCATE ( arg%pb(arg%no), arg%qb(arg%no), arg%rb(arg%no))
  ALLOCATE ( arg%pq1(arg%no), arg%pq2(arg%no), arg%pq3(arg%no), arg%pq4(arg%no))
  ALLOCATE ( arg%pq5(arg%no), arg%pq6(arg%no), arg%pq7(arg%no), arg%pq8(arg%no))
  
  
  
  arg%bia(:) = 0.0
  do k=1,arg%no
     read (511,*) & !'(I5,I5,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,I8)') &
          arg%flc(k), arg%par(k), &
          arg%lon(k), arg%lat(k), &
          arg%dpt(k), arg%tim(k), &
          arg%res(k), arg%err(k), arg%ino(k)
  end do
  close (511)


  ! DECOMMENT FOLLOWING TWO LINES TO MAKE FILTER TEST
  !arg%res(:) = 1
  !arg%err(:) = 1.d-2

  
  ! ---
  ! Initialise quality flag
  arg%flg(:) = 1
  
  ! ---
! Vertical interpolation parameters
  do k = 1,arg%no
     if(arg%flg(k).eq.1)then
        arg%kb(k) = grd%km-1
        do kk = 1,grd%km-1
           if( arg%dpt(k).ge.grd%dep(kk) .and. arg%dpt(k).lt.grd%dep(kk+1) ) then
              arg%kb(k) = kk
              arg%rb(k) = (arg%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
           endif
        enddo
     endif
  enddo
  
  
  ! residual check
  ! do k=1,arg%no
  !    if(arg%par(k).eq.1 .and. abs(arg%res(k)).gt.5.0) arg%flg(k) = 0
  ! enddo
    
  ! ---
  ! Count good observations
  arg%nc = 0
  do k=1,arg%no
     if(arg%flg(k).eq.1)then
        arg%nc = arg%nc + 1
     else
        arg%bia(k) = 0.
        arg%res(k) = 0.
        arg%inc(k) = 0.
        arg%b_a(k) = 0.
        arg%pq1(k) = 0.
        arg%pq2(k) = 0.
        arg%pq3(k) = 0.
        arg%pq4(k) = 0.
        arg%pq5(k) = 0.
        arg%pq6(k) = 0.
        arg%pq7(k) = 0.
        arg%pq8(k) = 0.
     endif
  enddo
  arg%flc(:) = arg%flg(:)
  
end subroutine get_obs_arg



subroutine int_par_arg
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Get interpolation parameters for a grid                              !
  !                                                                      !
  ! Version 1: S.Dobricic 2006                                           !
  !-----------------------------------------------------------------------
  
  use set_knd
  use drv_str
  use grd_str
  use obs_str
  
  implicit none
  
  INTEGER(i4)   ::  i, j, k
  INTEGER(i4)   ::  i1, j1, k1, idep
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y
  LOGICAL       ::  ins
  
  ins(i,i1) = i.ge.1 .and. i.lt.i1
  
  if(arg%no.gt.0) then
     
     arg%flc(:) = arg%flg(:)

     ! ---
     ! Horizontal interpolation parameters
     do k = 1,arg%no
        do j=1,grd%jm-1
           do i=1,grd%im-1
              if( grd%lat(i,j).le.arg%lat(k) .and. grd%lat(i,j+1).gt.arg%lat(k) .and.   &
                   grd%lon(i,j).le.arg%lon(k) .and. grd%lon(i+1,j).gt.arg%lon(k) ) then
                 j1 = j
                 i1 = i
                 q1 = j1 + (arg%lat(k) - grd%lat(i,j)) / (grd%lat(i,j+1) - grd%lat(i,j))
                 p1 = i1 + (arg%lon(k) - grd%lon(i,j)) / (grd%lon(i+1,j) - grd%lon(i,j))
              endif
           enddo
        enddo
        !     q1 = (arg%lat(k) - grd%lat(1,1)) / grd%dlt + 1.0
        !     j1 = int(q1)
        !     p1 = (arg%lon(k) - grd%lon(1,1)) / grd%dln + 1.0
        !     i1 = int(p1)
        if(ins(j1,grd%jm) .and. ins(i1,grd%im)) then
           arg%ib(k) = i1
           arg%jb(k) = j1
           arg%pb(k) = (p1-i1)
           arg%qb(k) = (q1-j1)
        else
           arg%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Undefine masked for multigrid
     do k = 1,arg%no
        if(arg%flc(k).eq.1)then
           i1 = arg%ib(k)
           j1 = arg%jb(k)
           idep = arg%kb(k)+1
           msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
           if(msk4.lt.1.) arg%flc(k) = 0
        endif
     enddo
     
     ! ---
     ! Horizontal interpolation parameters for each masked grid
     do k = 1,arg%no
        if(arg%flc(k) .eq. 1) then
           
           i1=arg%ib(k)
           p1=arg%pb(k)
           j1=arg%jb(k)
           q1=arg%qb(k)
           
           
           k1=arg%kb(k)
           div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
           arg%pq1(k) = grd%msk(i1,j1,k1)                                       &
                * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           arg%pq2(k) = grd%msk(i1+1,j1,k1)                                     &
                * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
           arg%pq3(k) = grd%msk(i1,j1+1,k1)                                     &
                * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           arg%pq4(k) = grd%msk(i1+1,j1+1,k1)                                   &
                * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           k1=arg%kb(k) + 1
           div_y =  (1.-q1) * max(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
                +    q1  * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
           div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
           arg%pq5(k) = grd%msk(i1,j1,k1)                                       &
                * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                * (1.-p1) * (1.-q1)                                     &
                /( div_x * div_y + 1.e-16 )
           arg%pq6(k) = grd%msk(i1+1,j1,k1)                                     &
                * max(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))             &
                *     p1  * (1.-q1)                                      &
                /( div_x * div_y + 1.e-16 )
           div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
           arg%pq7(k) = grd%msk(i1,j1+1,k1)                                     &
                * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                * (1.-p1) *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           arg%pq8(k) = grd%msk(i1+1,j1+1,k1)                                   &
                * max(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))         &
                *     p1  *     q1                                       &
                /( div_x * div_y + 1.e-16 )
           
           r1=arg%rb(k)
           arg%pq1(k) = (1.-r1) * arg%pq1(k)
           arg%pq2(k) = (1.-r1) * arg%pq2(k)
           arg%pq3(k) = (1.-r1) * arg%pq3(k)
           arg%pq4(k) = (1.-r1) * arg%pq4(k)
           arg%pq5(k) =     r1  * arg%pq5(k)
           arg%pq6(k) =     r1  * arg%pq6(k)
           arg%pq7(k) =     r1  * arg%pq7(k)
           arg%pq8(k) =     r1  * arg%pq8(k)
           
        endif
     enddo
     
     
     ! ---
     ! Count good observations
     arg%nc = 0
     do k=1,arg%no
        if(arg%flc(k).eq.1)then
           arg%nc = arg%nc + 1
        endif
     enddo
     write(drv%dia,*)'Real number of ARGO observations: ',arg%nc
     
  endif
  
end subroutine int_par_arg
