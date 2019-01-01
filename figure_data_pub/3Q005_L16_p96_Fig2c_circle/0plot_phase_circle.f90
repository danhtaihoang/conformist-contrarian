!!!! Chuong trinh tinh ve phase tren 1 duong tron

   PROGRAM phase_circle
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER:: N=4096
   INTEGER (KIND=8),PARAMETER:: Na1=50,Nb1=50
   REAL (KIND=8),PARAMETER::  r0=10.,delr=1.

   INTEGER (KIND=8):: i,j,ia,ib,na,nb,na2,nb2
   REAL    (KIND=8):: pi,r,rx,ry,rdn_a,rdn_b,rdn_ra,rdn_rb

   INTEGER (KIND=8),DIMENSION(N) :: site_a,site_b
   REAL (KIND=8),DIMENSION(N) :: Ks,S,Sa,Sb

!!!=======================================================================================
!!!======================================================================================= 

   pi=acos(-1.) ; site_a(:)=0 ; site_b(:)=0

   OPEN(unit=11,file='phi_equi.txt')
   OPEN(unit=12,file='phi_a.txt')
   OPEN(unit=13,file='phi_b.txt')

   S(:)=0. ; Sa(:)=0. ; Sb(:)=0.

!!! Doc gia tri
   ia=0 ; ib=0 
   DO i=1,N
      READ(11,*)j,Ks(j),S(j)

      S(i)=S(i)-int(S(i)/2./pi)*2.*pi
   
      IF (S(i)<0.) THEN
         S(i)=S(i)+2.*pi
      END IF

      IF (Ks(i)>0.) THEN
         ia=ia+1
         Sa(ia)=S(i)
      ELSE
         ib=ib+1
         Sb(ib)=S(i)
      END IF
   END DO

   na=ia ; nb=ib
   WRITE(*,*)na,nb

!!! Select na1 at random
   na2=0
   DO WHILE (na2<na1)
      CALL random_number(rdn_a)
      ia=int(rdn_a*(na-1)+1)

      IF (site_a(ia)==0) THEN

         CALL random_number(rdn_ra)
         r=r0+delr*(rdn_ra-0.5)
         rx=r*cos(Sa(ia))
         ry=r*sin(Sa(ia))

         WRITE(12,*)ia,Sa(ia),rx,ry
         na2=na2+1
         site_a(ia)=1
      END IF

   END DO

!!! Select nb1 at random
   nb2=0
   DO WHILE (nb2<nb1)
      CALL random_number(rdn_b)
      ib=int(rdn_b*(nb-1)+1)

      IF (site_b(ib)==0) THEN
         CALL random_number(rdn_rb)
         r=r0+delr*(rdn_rb-1.)
         rx=r*cos(Sb(ib))
         ry=r*sin(Sb(ib))

         WRITE(13,*)ib,Sb(ib),rx,ry
         nb2=nb2+1
         site_b(ib)=1
      END IF

   END DO

   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
    
   END PROGRAM phase_circle

