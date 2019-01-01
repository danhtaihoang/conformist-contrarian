   PROGRAM value_delta_R
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=27
   INTEGER (KIND=8):: i
   REAL (KIND=8):: Rmin,Rmax,dev_R
   REAL    (KIND=8),DIMENSION(n):: R,pa,delta_R

!!!=======================================================================================
!!!======================================================================================= 

   OPEN(unit=11,file='R_L24.txt')
   OPEN(unit=12,file='delta_R_Q005_L24.txt')

   DO i=1,n
      READ(11,*)pa(i),R(i),Rmin,Rmax,dev_R
   END DO

   delta_R(:)=0.
   DO i=1,n-1
      delta_R(i)=(R(i+1)-R(i))/(pa(i+1)-pa(i))
   END DO

   DO i=1,n
      WRITE(12,*)pa(i),delta_R(i)
   END DO

   CLOSE(11)
   CLOSE(12)
    
   END PROGRAM value_delta_R

