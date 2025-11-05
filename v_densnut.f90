subroutine v_densnut(NutArray)
  
  
  !---------------------------------------------------------------------------
  !                                                                          !
  !    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
  !                                                                          !
  !    This file is part of OceanVar.                                        !
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
  !    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
  !                                                                          !
  !---------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  !                                                                      !
  ! Apply density nutrient increments and unterpolare aon vertical grid  !
  !                                                                      !
  ! Version 1: A. Teruzzi 2025                                           !
  !-----------------------------------------------------------------------
  
  
  use set_knd
  use grd_str
  ! use eof_str
  ! use obs_str
  use mpi_str
  ! use filenames
  use drv_str
  use bio_str
  use dnc_str

  implicit none
  
  INTEGER(i4)   ::  i, j, k, kk, my_km
  REAL(r8), DIMENSION(grd%jm,grd%km)  :: slicevar
  REAL(r8), DIMENSION(grd%im+1,grd%jm,grd%km)  :: N3nExt
  REAL(r8) :: NutArray(grd%im,grd%jm,grd%km)
  REAL(8) :: nit_dinc

  my_km = grd%km
  
  ! call EXTEND_2D( grd%n3n_ad, grd%km, N3nExtended_3d )

  dnc%k = 0
  N3nExt(:,:,:) = 0.0

  do kk = 1,dnc%no

    i=dnc%ib(kk)
    j=dnc%jb(kk)
    k=dnc%kb(kk)

    if(dnc%flc(kk).eq.1) then

          dnc%k = dnc%k + 1
          nit_dinc = dnc%inc(dnc%k) * dnc%corr(dnc%k) * dnc%err(dnc%k) / dnc%std(dnc%k)

          N3nExt(i  ,j  ,k  ) = N3nExt(i  ,j  ,k  ) + dnc%pq1(kk) * nit_dinc
          N3nExt(i+1,j  ,k  ) = N3nExt(i+1,j  ,k  ) + dnc%pq2(kk) * nit_dinc
          N3nExt(i  ,j+1,k  ) = N3nExt(i  ,j+1,k  ) + dnc%pq3(kk) * nit_dinc
          N3nExt(i+1,j+1,k  ) = N3nExt(i+1,j+1,k  ) + dnc%pq4(kk) * nit_dinc
          N3nExt(i  ,j  ,k+1) = N3nExt(i  ,j  ,k+1) + dnc%pq5(kk) * nit_dinc
          N3nExt(i+1,j  ,k+1) = N3nExt(i+1,j  ,k+1) + dnc%pq6(kk) * nit_dinc
          N3nExt(i  ,j+1,k+1) = N3nExt(i  ,j+1,k+1) + dnc%pq7(kk) * nit_dinc
          N3nExt(i+1,j+1,k+1) = N3nExt(i+1,j+1,k+1) + dnc%pq8(kk) * nit_dinc


    endif

  enddo

!  we apply contribution in grd%variable

  slicevar(:,1:grd%km) = NutArray(1,:,:)
  call ADD_PREVCORE_CONTRIB(N3nExt, grd%km, NutArray, slicevar)
  ! call ADD_PREVCORE_CONTRIB(N3nExtended_3d,  grd%km, grd%N3n_ad, grd%n3n_ad(1,:,:))



end subroutine v_densnut
