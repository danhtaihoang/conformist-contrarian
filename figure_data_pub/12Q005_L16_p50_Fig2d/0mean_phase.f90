!!!! Chuong trinh tinh standard deviation

   PROGRAM standard_deviation
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::N=32768
   INTEGER (KIND=8):: i,j,ia,ib
   REAL    (KIND=8):: Sa_av,dev_a,Sb_av,dev_b,pi

   REAL (KIND=8),DIMENSION(N) :: Ks,S,Sa,Sb

!!!=======================================================================================
!!!======================================================================================= 

   pi=acos(-1.)

   OPEN(unit=11,file='phi_equi.dat')
   OPEN(unit=12,file='phi_av.txt')


   S(:)=0. ; Sa(:)=0. ; Sb(:)=0.

!!! Doc gia tri
   ia=0 ; ib=0 
   DO i=1,N
      READ(11,*)j,Ks(j),S(j)

      S(i)=S(i)-int(S(i)/2./pi)*2.*pi
   
      IF (S(i)<0.) THEN
         S(i)=S(i)+2.*pi
      END IF

      IF (Ks(j)>0.) THEN
         ia=ia+1
         Sa(ia)=S(j)
      ELSE
         ib=ib+1
         Sb(ib)=S(j)
      END IF
   END DO

!!! Tinh gia tri trung binh
   Sa_av=0.   
   DO i=1,ia
      Sa_av=Sa_av+Sa(i)
   END DO
   Sa_av=Sa_av/ia

   Sb_av=0.   
   DO i=1,ib
      Sb_av=Sb_av+Sb(i)
   END DO
   Sb_av=Sb_av/ib

!!! Tinh standard deviation
   dev_a=0.   
   DO i=1,ia
      dev_a=dev_a+(Sa(i)-Sa_av)**2.
   END DO
   dev_a=dev_a/ia

   dev_b=0.   
   DO i=1,ib
      dev_b=dev_b+(Sb(i)-Sb_av)**2.
   END DO
   dev_b=dev_b/ib

   WRITE(12,*)N,ia,ib
   WRITE(12,*)Sa_av,dev_a,Sb_av,dev_b

   CLOSE(11)
   CLOSE(12)
    
   END PROGRAM standard_deviation

