      PROGRAM divide_by
      IMPLICIT NONE

      INTEGER (KIND=8),PARAMETER::n=2601

      INTEGER (KIND=8):: i
      REAL    (KIND=8):: Jbb,Jab,cost

!!!=======================================================================================
!!!======================================================================================= 

      OPEN(unit=11,file='a.dat')
      OPEN(unit=12,file='a_new.dat')

      DO i=1,n
            READ(11,*)Jbb,Jab,cost
            WRITE(12,*)Jbb/100,Jab/100,cost
      END DO
 
      CLOSE(11)
      CLOSE(12)
    
      END PROGRAM divide_by

