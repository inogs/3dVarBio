      REAL(r8) FUNCTION mean_rad(k,rad_xy)
      use rcfl
      use grd_str
      
      IMPLICIT NONE
      INTEGER(i4):: j,i,count
      INTEGER(i4),INTENT(IN):: k
      REAL(r8),INTENT(IN):: rad_xy(grd%im,grd%jm,1)
      REAL(r8):: meanxy
     
      print*, k
      count=0
      meanxy=0
      do j=1,grd%jm
        do i=1,grd%im
          if(rad_xy(i,j,1)>0.0001) then
            count=count+1
            meanxy=meanxy+rad_xy(i,j,1)
          end if
        enddo
      enddo
      mean_rad=meanxy/count
      return     
      END FUNCTION mean_rad
