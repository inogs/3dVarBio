      SUBROUTINE mod_trj_ad( jpi,jpj,umod,vmod,e1u,e2v,flg,        &
                             jpt,jpn,pimod,pjmod,ptime,uadj,vadj,xadj,yadj )

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Vincent Taillandier, Locean, Paris                     !
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
! Trajectory model (adjoint)                                           !
!                                                                      !
! Version 1: V. Taillandier 2007                                       !
!                                                                      !
!-----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Arguments
!! =========
      implicit none
      INTEGER*8 :: jpt,jpn
      INTEGER*8 :: jpi,jpj
      INTEGER*8 :: flg(jpn)
      REAL*8, DIMENSION(jpt+1,jpn) :: pimod,pjmod
      REAL*8, DIMENSION(jpn) :: ptime
      REAL*8, DIMENSION(jpi,jpj) :: umod,vmod,e1u,e2v
      REAL*8, DIMENSION(jpn) :: xadj,yadj
!!----------------------------------------------------------------------
!! Local declarations
!! ==================
      INTEGER :: jn,jt
      REAL*8, DIMENSION(jpi,jpj) :: uadj,vadj
      REAL*8 :: zpin,zpjn,zuadj,zvadj,zrdt,zpiadj,zpjadj
!!----------------------------------------------------------------------
!*... Trajectories computation
      uadj(:,:) = 0.e0
      vadj(:,:) = 0.e0
      DO jn = 1, jpn
       if(flg(jn).eq.1)then
        zpiadj = xadj(jn)
        zpjadj = yadj(jn)
        zrdt = ptime(jn)*3600./(jpt)
        DO jt = jpt, 1, -1
          zpin = pimod(jt,jn)
          zpjn = pjmod(jt,jn)
          zuadj  = zpiadj
          zvadj  = zpjadj
          zuadj  = zuadj * zrdt
          zvadj  = zvadj * zrdt
          CALL floitpadj( jpi,jpj,umod,vmod,uadj,vadj,e1u,e2v, &
                          zpin,zpjn,zpiadj,zpjadj,zuadj,zvadj )
        END DO
       endif
      END DO

      RETURN
      END SUBROUTINE mod_trj_ad

!!======================================================================

      SUBROUTINE floitpadj( jpi,jpj,umod,vmod,uadj,vadj,e1u,e2v, &
                             pifl,pjfl,piflad,pjflad,puflad,pvflad )
!!----------------------------------------------------------------------
!! Arguments
!! =========
      implicit none
      INTEGER*8 :: jpi,jpj
      REAL*8 :: pifl,pjfl,piflad,pjflad,puflad,pvflad
      REAL*8, DIMENSION(jpi,jpj) :: umod,vmod,uadj,vadj,e1u,e2v
!!----------------------------------------------------------------------
!! Local declarations
!! ==================
      INTEGER :: iil,ijl,jind1,jind2
      INTEGER, DIMENSION(2) :: iid,ijd
      REAL*8, DIMENSION(2) :: zlagx,zlagy,zlagxad,zlagyad
      REAL*8, DIMENSION(2,2) :: zuv,zuvad
!!----------------------------------------------------------------------
 
!! 1. Interpolation of the zonal velocity
!! ======================================
 
!*... Neighbooring points (background)
      iil = INT(pifl-.5)
      ijl = INT(pjfl   )
      DO jind1 = 1, 2
        iid(jind1) = iil + jind1 - 1
        ijd(jind1) = ijl + jind1 - 1
      END DO
 
!*... Lagrange coefficients (background)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        IF( jind1.NE.jind2 ) THEN
          zlagx(jind1) = ( pifl - ((iid(jind2))+.5) ) / ( iid(jind1)-iid(jind2) )
          zlagy(jind1) = ( pjfl -  (ijd(jind2))     ) / ( ijd(jind1)-ijd(jind2) )
        ENDIF
      END DO
      END DO
 
!*... Value of the zonal velocity (background)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        zuv(jind1,jind2) = umod(iid(jind1),ijd(jind2)) / e1u(iid(jind1),ijd(jind2))
      END DO
      END DO
 
!*... Interpolation of the zonal velocity
      zlagxad(:) = 0.e0
      zlagyad(:) = 0.e0
      zuvad(:,:) = 0.e0
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        zlagxad(jind1)     = zlagxad(jind1) + puflad * zuv(jind1,jind2) * zlagy(jind2)
        zlagyad(jind2)     = zlagyad(jind2) + puflad * zuv(jind1,jind2) * zlagx(jind1)
        zuvad(jind1,jind2) = zuvad(jind1,jind2) + puflad * zlagx(jind1) * zlagy(jind2)
      END DO
      END DO
      puflad = 0.e0
 
!*... Value of the zonal velocity (adjoint)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        uadj(iid(jind1),ijd(jind2)) = uadj(iid(jind1),ijd(jind2)) &
                                    + zuvad(jind1,jind2) / e1u(iid(jind1),ijd(jind2))
      END DO
      END DO
 
!*... Lagrange coefficients (adjoint)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        IF( jind1.NE.jind2 ) THEN
          piflad = piflad + zlagxad(jind1) / ( iid(jind1)-iid(jind2) )
!          pjflad = pjflad + zlagyad(jind1) / ( iid(jind1)-iid(jind2) )
          pjflad = pjflad + zlagyad(jind1) / ( ijd(jind1)-ijd(jind2) )
        ENDIF
      END DO
      END DO
 
!! 2. Interpolation of the meridian velocity
!! =========================================
 
!*... Neighbooring points (background)
      iil = INT(pifl   )
      ijl = INT(pjfl-.5)
      DO jind1 = 1, 2
        iid(jind1) = iil + jind1 - 1
        ijd(jind1) = ijl + jind1 - 1
      END DO
 
!*... Lagrange coefficients (background)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        IF( jind1.NE.jind2 ) THEN
          zlagx(jind1) = ( pifl -  (iid(jind2))     ) / ( iid(jind1)-iid(jind2) )
          zlagy(jind1) = ( pjfl - ((ijd(jind2))+.5) ) / ( ijd(jind1)-ijd(jind2) )
        ENDIF
      END DO
      END DO
 
!*... Value of the meridian velocity (background)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        zuv(jind1,jind2) = vmod(iid(jind1),ijd(jind2)) / e2v(iid(jind1),ijd(jind2))
      END DO
      END DO
 
!*... Interpolation of the meridian velocity
      zlagxad(:) = 0.e0
      zlagyad(:) = 0.e0
      zuvad(:,:) = 0.e0
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        zlagxad(jind1)     = zlagxad(jind1) + pvflad * zuv(jind1,jind2) * zlagy(jind2)
        zlagyad(jind2)     = zlagyad(jind2) + pvflad * zuv(jind1,jind2) * zlagx(jind1)
        zuvad(jind1,jind2) = zuvad(jind1,jind2) + pvflad * zlagx(jind1) * zlagy(jind2)
      END DO
      END DO
      pvflad = 0.e0
 
!*... Value of the meridian velocity (adjoint)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        vadj(iid(jind1),ijd(jind2)) = vadj(iid(jind1),ijd(jind2)) &
                                    + zuvad(jind1,jind2) / e2v(iid(jind1),ijd(jind2))
      END DO
      END DO
 
!*... Lagrange coefficients (adjoint)
      DO jind1 = 1, 2
      DO jind2 = 1, 2
        IF( jind1.NE.jind2 ) THEN
          piflad = piflad + zlagxad(jind1) / ( iid(jind1)-iid(jind2) )
!          pjflad = pjflad + zlagyad(jind1) / ( iid(jind1)-iid(jind2) )
          pjflad = pjflad + zlagyad(jind1) / ( ijd(jind1)-ijd(jind2) )
        ENDIF
      END DO
      END DO

      RETURN
      END SUBROUTINE floitpadj
