   PROGRAM value_delta_R
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=27
   INTEGER (KIND=8):: i
   REAL    (KIND=8),DIMENSION(n):: ra,rb,R,v,pa,dev,delta_R

!!!=======================================================================================
!!!======================================================================================= 

   OPEN(unit=11,file='average_Q005_L12.txt')
   OPEN(unit=12,file='delta_R_Q005_L12.txt')

   DO i=1,n
      READ(11,*)pa(i),ra(i),rb(i),R(i),v(i),dev(i)
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

