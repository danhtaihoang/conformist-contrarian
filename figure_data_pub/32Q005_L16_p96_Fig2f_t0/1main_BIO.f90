!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!  Hogil Kim Memorial Building #501 POSTECH,
!  San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!  E-mail: hoangdanhtai@gmail.com    Personal site: hoangdanhtai.com 
!-----------------------------------------------------------------------------------------!
!!! 2014.03.15: Interact with nearest neighbor, using 1 index
!!! 2014.04.14: Bo sung structure 3D
!!! 2014.04.30: Bo sung GS1 config (CHUA hoan chinh)
!!! 2014.06.10: Sua lai cach tinh Zx, Zy de chay nhanh hon
!!! 2014.06.12: Plot Wave Speed
!!! 2014.08.22: Bo sung average phase, delta average phase,pa1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
   PROGRAM main_BIO
   IMPLICIT NONE

   CHARACTER (LEN=150):: INI_CONFIG,LOAD_PHI
   CHARACTER (LEN=3)  :: SAt
   REAL    (KIND=8),PARAMETER :: nul=0.,aaa=4.0
   CHARACTER (256)    :: Ligne20,Ligne21,Ligne22,Ligne23,Ligne40

   INTEGER (KIND=8) :: i,j,i1,i2,nx,ny,nz,nx_p,ny_p,nz_p,na,nb
   INTEGER (KIND=8):: i_loop,i_loop1,i_av,n_equi1,n_equi2,n_average,i_times,n_times
   INTEGER (KIND=8) :: time_2,time_3,time_5,time_6,time_7,natom
   
   REAL (KIND=8) :: rdn_config,n_atom,number_a,number_b,rdn_s,pa1
   REAL (KIND=8) :: pi,delt,delt1,pa,omega,A,Ka,Kb,Zax,Zay,Za,Zbx,Zby,Zb,Za_av,Zb_av
   REAL (KIND=8) :: Zx,Zy,Z_total,Z_total_av,run_time,S_tmp,Ks_load,S_load,WS,WS_av
   REAL (KIND=8) :: phi_a,phi_b,delta_ab,phi_total

   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: S,delS,Ks,Ks1,S1,S2
   INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn1,nn2,nn3,nn4,nn5,nn6,x,y,z

   !!! Khai bao phan Histogram
   INTEGER (KIND=8),PARAMETER :: n_histo=101
   INTEGER (KIND=8) :: i_histo
   REAL    (KIND=8) :: del_histo

   REAL (KIND=8),DIMENSION(0:n_histo) :: S_histo,H1_histo,H2_histo
      
!!!=======================================================================================   
!!!=======================================================================================
   CALL system('rm config_ini_3D.pdb')
   CALL system('rm *.dat*')

   CALL ini_rdm_number()
   CALL read_input_parameter()

   nx_p=nx+1 ; ny_p=ny+1 ; nz_p=nz+1
   natom=nx*ny*nz
   n_atom=real(natom)
   number_a=real(nint(n_atom*pa)) ; number_b=n_atom-number_a
   pa1=number_a/n_atom

   !WRITE(*,*)number_a,number_b

   ALLOCATE(S(natom),delS(natom),Ks(natom),Ks1(natom),S1(natom),S2(natom))
   ALLOCATE(nn1(natom),nn2(natom),nn3(natom),nn4(natom),nn5(natom),nn6(natom))
   ALLOCATE(x(natom),y(natom),z(natom))

   delt1=delt+(1./2.)*delt**2. ! +(1./6.)*delt**3.+(1./24.)*delt**4.
   
   pi=acos(-1.)

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   OPEN(unit=20,file='phi_equi.dat')
   OPEN(unit=21,file='r_time.dat')
   OPEN(unit=22,file='average.dat')
   OPEN(unit=23,file='phi_time.dat')

   IF (LOAD_PHI=='YES') THEN
      WRITE(*,*)'Load phi equi'
      CALL load_phi_equi()
   ELSE
      CALL load_config_ini()
      CALL generate_ini_s()
   END IF

   Ks1(:)=Ks(:)/6.
 
   CALL index_coordinate()
   CALL nearest_neighbor()
   !!CALL write_config_ini_3D()

   !CALL value_thermal()
   
   CALL equi_lattice1()
   WRITE(*,*)'Finish Equilibre'

   CALL average_thermal()

   CALL save_phi_equi()

   CALL computation_time()

   CONTAINS
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ini_rdm_number()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time
   INTEGER (KIND=8),DIMENSION(50) :: seed

   CALL DATE_AND_TIME(values=time)     ! Get the current time
   seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
   CALL RANDOM_SEED(PUT=seed)

   time_2=time(2) ; time_3=time(3) ; time_5=time(5) ; time_6=time(6) ; time_7=time(7)

   END SUBROUTINE ini_rdm_number
   
!!!=======================================================================================  
   SUBROUTINE computation_time()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time1
   INTEGER (KIND=8) :: run_date,run_hour,run_minute,run_second

   OPEN (90,file='time_run.dat')

   CALL DATE_AND_TIME(values=time1)     ! Get the current time
      
   run_time = (time1(2)-time_2)*1296000.+(time1(3)-time_3)*86400.&
             +(time1(5)-time_5)*3600.+(time1(6)-time_6)*60.+(time1(7)-time_7) !! second (s)
   run_date=int(run_time/86400.)
   run_hour=int(run_time/3600.-run_date*24.)
   run_minute=int(run_time/60.-run_date*1440.-run_hour*60.)
   run_second=int(run_time-run_date*86400.-run_hour*3600.-run_minute*60.)

   WRITE(90,*)'run_date  :',run_date
   WRITE(90,*)'run_hour  :',run_hour
   WRITE(90,*)'run_minute:',run_minute
   WRITE(90,*)'run_second:',run_second

   CLOSE(90)

   END SUBROUTINE computation_time   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE read_input_parameter()
   IMPLICIT NONE

   CHARACTER (LEN=150) :: tamp
   OPEN(11,file='1parameter.in')   
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(A10))')   tamp, LOAD_PHI
   READ(11, '(A30,(A10))')   tamp, INI_CONFIG
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(I5))')    tamp, nx
   READ(11, '(A30,(I5))')    tamp, ny
   READ(11, '(A30,(I5))')    tamp, nz
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, pa
   READ(11, '(A30,(F12.6))') tamp, omega
   READ(11, '(A30,(F12.6))') tamp, A
   READ(11, '(A30,(F12.6))') tamp, Ka
   READ(11, '(A30,(F12.6))') tamp, Kb
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, delt
   READ(11, '(A30,(I12))')   tamp,n_equi1
   READ(11, '(A30,(I12))')   tamp,n_equi2
   READ(11, '(A30,(I12))')   tamp,n_average
   READ(11, '(A30,(I12))')   tamp,n_times
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp

   CLOSE(11) 

   END SUBROUTINE read_input_parameter

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice1()
   IMPLICIT NONE

   DO i_loop=1,n_equi1
      CALL equi_lattice()
   END DO

   END SUBROUTINE equi_lattice1
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE load_config_ini()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE load_config_ini()
   IMPLICIT NONE 

   Ks(:)=0.
   
   !!! GS random : -----------------------------------------------
   IF (INI_CONFIG == 'NO') THEN
      na=0 ; nb=0
      
      DO WHILE (na<number_a)
            i=0

            DO WHILE(i==0)
                  CALL random_number(rdn_config)
                  i=int(rdn_config*n_atom)
            ENDDO

            IF (Ks(i)==0.) THEN                      
                  Ks(i)=Ka
                  na=na+1
            END IF

      END DO

      DO i=1,natom

         IF (Ks(i)==0.) THEN   
            Ks(i)=Kb
            nb=nb+1
         END IF

      ENDDO

   WRITE(*,*)'GS random, na=',na,'nb=',nb

   END IF

   !!! GS1: 8 contrarians randomly --------------------
   IF (INI_CONFIG == 'GS1') THEN
      na=0 ; nb=8
      
      Ks(1)=Kb  ;  Ks(64)=Kb 
      Ks(253)=Kb ; Ks(288)=Kb 
      Ks(512)=Kb ; Ks(557)=Kb 
      Ks(867)=Kb ; Ks(852)=Kb 

      DO i=1,natom

         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF

      ENDDO

   WRITE(*,*)'GS1, na=',na,'nb=',nb

   END IF

   !!! GS2: 1 clusters + 4 contrarians randomly --------------------
   IF (INI_CONFIG == 'GS2') THEN
      na=0 ; nb=8
      
      Ks(5)=Kb  ;  Ks(6)=Kb ; Ks(7)=Kb ; Ks(16)=Kb
      Ks(512)=Kb ; Ks(557)=Kb 
      Ks(867)=Kb ; Ks(852)=Kb 

      DO i=1,natom

         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF

      ENDDO

   WRITE(*,*)'GS2, na=',na,'nb=',nb

   END IF

   !!! GS3: 2 clusters --------------------
   IF (INI_CONFIG == 'GS3') THEN
      na=0 ; nb=8
      
      Ks(5)=Kb  ;  Ks(6)=Kb ; Ks(7)=Kb ; Ks(16)=Kb
      Ks(512)=Kb ; Ks(513)=Kb ;  Ks(514)=Kb ;  Ks(523)=Kb

      DO i=1,natom

         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF

      ENDDO

   WRITE(*,*)'GS3, na=',na,'nb=',nb

   END IF

   !!! GS4: 1 clusters --------------------
   IF (INI_CONFIG == 'GS4') THEN
      na=0 ; nb=4
      
      Ks(5)=Kb  ;  Ks(6)=Kb ; Ks(7)=Kb ; Ks(16)=Kb

      DO i=1,natom

         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF

      ENDDO

   WRITE(*,*)'GS4, na=',na,'nb=',nb

   END IF

   !!! GS23: 16 clusters for L=20 --------------------
   IF (INI_CONFIG == 'GS23') THEN
      na=0 ; nb=64
      
      Ks(4)=Kb  ;  Ks(5)=Kb ; Ks(6)=Kb ; Ks(25)=Kb
      Ks(75)=Kb ; Ks(76)=Kb ;  Ks(77)=Kb ;  Ks(96)=Kb
      Ks(204)=Kb  ;  Ks(205)=Kb ; Ks(206)=Kb ; Ks(225)=Kb
      Ks(275)=Kb ; Ks(276)=Kb ;  Ks(277)=Kb ;  Ks(296)=Kb

      Ks(2004)=Kb  ;  Ks(2005)=Kb ; Ks(2006)=Kb ; Ks(2025)=Kb
      Ks(2075)=Kb ; Ks(2076)=Kb ;  Ks(2077)=Kb ;  Ks(2096)=Kb
      Ks(2204)=Kb  ;  Ks(2205)=Kb ; Ks(2206)=Kb ; Ks(2225)=Kb
      Ks(2275)=Kb ; Ks(2276)=Kb ;  Ks(2277)=Kb ;  Ks(2296)=Kb

      Ks(4004)=Kb  ;  Ks(4005)=Kb ; Ks(4006)=Kb ; Ks(4025)=Kb
      Ks(4075)=Kb ; Ks(4076)=Kb ;  Ks(4077)=Kb ;  Ks(4096)=Kb
      Ks(4204)=Kb  ;  Ks(4205)=Kb ; Ks(4206)=Kb ; Ks(4225)=Kb
      Ks(4275)=Kb ; Ks(4276)=Kb ;  Ks(4277)=Kb ;  Ks(4296)=Kb

      Ks(6004)=Kb  ;  Ks(6005)=Kb ; Ks(6006)=Kb ; Ks(6025)=Kb
      Ks(6075)=Kb ; Ks(6076)=Kb ;  Ks(6077)=Kb ;  Ks(6096)=Kb
      Ks(6204)=Kb  ;  Ks(6205)=Kb ; Ks(6206)=Kb ; Ks(6225)=Kb
      Ks(6275)=Kb ; Ks(6276)=Kb ;  Ks(6277)=Kb ;  Ks(6296)=Kb

      DO i=1,natom
         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF
      ENDDO

   WRITE(*,*)'GS23, na=',na,'nb=',nb

   END IF

   !!! GS22: 8 clusters + 32 randoms for L=20 --------------------
   IF (INI_CONFIG == 'GS22') THEN
      na=0 ; nb=64
      
      Ks(4)=Kb ; Ks(53)=Kb ; Ks(122)=Kb ; Ks(130)=Kb ; Ks(185)=Kb ; Ks(195)=Kb ; Ks(221)=Kb ; Ks(291)=Kb
      Ks(1204)=Kb ; Ks(1253)=Kb ; Ks(1322)=Kb ; Ks(1330)=Kb ; Ks(1385)=Kb ; Ks(1395)=Kb ; Ks(1421)=Kb ; Ks(1491)=Kb
      Ks(4004)=Kb ; Ks(4053)=Kb ; Ks(4122)=Kb ; Ks(4130)=Kb ; Ks(4185)=Kb ; Ks(4195)=Kb ; Ks(4221)=Kb ; Ks(4291)=Kb
      Ks(5204)=Kb ; Ks(5253)=Kb ; Ks(5322)=Kb ; Ks(5330)=Kb ; Ks(5385)=Kb ; Ks(5395)=Kb ; Ks(5421)=Kb ; Ks(5491)=Kb

      Ks(2004)=Kb  ;  Ks(2005)=Kb ; Ks(2006)=Kb ; Ks(2025)=Kb
      Ks(2075)=Kb ; Ks(2076)=Kb ;  Ks(2077)=Kb ;  Ks(2096)=Kb
      Ks(2204)=Kb  ;  Ks(2205)=Kb ; Ks(2206)=Kb ; Ks(2225)=Kb
      Ks(2275)=Kb ; Ks(2276)=Kb ;  Ks(2277)=Kb ;  Ks(2296)=Kb

      Ks(6004)=Kb  ;  Ks(6005)=Kb ; Ks(6006)=Kb ; Ks(6025)=Kb
      Ks(6075)=Kb ; Ks(6076)=Kb ;  Ks(6077)=Kb ;  Ks(6096)=Kb
      Ks(6204)=Kb  ;  Ks(6205)=Kb ; Ks(6206)=Kb ; Ks(6225)=Kb
      Ks(6275)=Kb ; Ks(6276)=Kb ;  Ks(6277)=Kb ;  Ks(6296)=Kb

      DO i=1,natom
         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF
      ENDDO

   WRITE(*,*)'GS22, na=',na,'nb=',nb

   END IF

   !!! GS24: 8 clusters for L=20 --------------------
   IF (INI_CONFIG == 'GS24') THEN
      na=0 ; nb=32
      
      Ks(2004)=Kb  ;  Ks(2005)=Kb ; Ks(2006)=Kb ; Ks(2025)=Kb
      Ks(2075)=Kb ; Ks(2076)=Kb ;  Ks(2077)=Kb ;  Ks(2096)=Kb
      Ks(2204)=Kb  ;  Ks(2205)=Kb ; Ks(2206)=Kb ; Ks(2225)=Kb
      Ks(2275)=Kb ; Ks(2276)=Kb ;  Ks(2277)=Kb ;  Ks(2296)=Kb

      Ks(6004)=Kb  ;  Ks(6005)=Kb ; Ks(6006)=Kb ; Ks(6025)=Kb
      Ks(6075)=Kb ; Ks(6076)=Kb ;  Ks(6077)=Kb ;  Ks(6096)=Kb
      Ks(6204)=Kb  ;  Ks(6205)=Kb ; Ks(6206)=Kb ; Ks(6225)=Kb
      Ks(6275)=Kb ; Ks(6276)=Kb ;  Ks(6277)=Kb ;  Ks(6296)=Kb

      DO i=1,natom
         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF
      ENDDO

   WRITE(*,*)'GS24, na=',na,'nb=',nb

   END IF

   !!! GS25: 2 clusters for L=20 --------------------
   IF (INI_CONFIG == 'GS25') THEN
      na=0 ; nb=8
      
      Ks(2004)=Kb  ;  Ks(2005)=Kb ; Ks(2006)=Kb ; Ks(2025)=Kb
      Ks(6275)=Kb ; Ks(6276)=Kb ;  Ks(6277)=Kb ;  Ks(6296)=Kb

      DO i=1,natom
         IF (Ks(i)==0.) THEN   
            Ks(i)=Ka
            na=na+1
         END IF
      ENDDO

   WRITE(*,*)'GS25, na=',na,'nb=',nb

   END IF

   END SUBROUTINE load_config_ini
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE generate_ini_s() !!! Generate initial value for S:
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE generate_ini_s()
   IMPLICIT NONE 
   
   DO i=1,natom

      CALL random_number(rdn_s)
      S(i)=2.*pi*rdn_s
      !S(i)=0.

   ENDDO

   END SUBROUTINE generate_ini_s

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE save_phi_equi()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE save_phi_equi()
   IMPLICIT NONE 
   
   DO i=1,natom
      WRITE(Ligne20,*)i,Ks(i),S(i)
      WRITE(20,'(a)') trim(Ligne20)
   ENDDO

   CLOSE(20)

   END SUBROUTINE save_phi_equi

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE load_phi_equi()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE load_phi_equi()
   IMPLICIT NONE 

   OPEN(unit=30,file='phi_equi.txt')

   DO i_loop=1,natom
      READ(30,*)i,Ks_load,S_load
      Ks(i)=Ks_load
      S(i)=S_load
 
   END DO

   END SUBROUTINE load_phi_equi

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE write_config_ini_3D()
   IMPLICIT NONE
 
   OPEN(unit=12,file='config_ini_3D.pdb')

   DO i=1,natom
            
      IF (Ks(i)==Ka) THEN 
         SAt='S'  !! yellow
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
                              
      IF (Ks(i)==Kb) THEN 
         SAt='C'  !! white
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul

      ELSE
           SAt='Na'  !! blue
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul

      END IF
      END IF

   END DO

   CLOSE(12)

   END SUBROUTINE write_config_ini_3D

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice()
   IMPLICIT NONE
     
   DO i=1,natom
      S_tmp=S(i)
      delS(i)=Ks1(i)*(sin(S(nn1(i))-S_tmp)+sin(S(nn2(i))-S_tmp)+sin(S(nn3(i))-S_tmp)&
       +sin(S(nn4(i))-S_tmp)+sin(S(nn5(i))-S_tmp)+sin(S(nn6(i))-S_tmp)) +omega-A*sin(S_tmp) 
   END DO

   DO i=1,natom
      S(i)=S(i)+delS(i)*delt1
   END DO

   END SUBROUTINE equi_lattice
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE value_thermal()
   IMPLICIT NONE

!!!===========================================================
   DO i=1,natom
      S_tmp=S(i)
      delS(i)=Ks1(i)*(sin(S(nn1(i))-S_tmp)+sin(S(nn2(i))-S_tmp)+sin(S(nn3(i))-S_tmp)&
       +sin(S(nn4(i))-S_tmp)+sin(S(nn5(i))-S_tmp)+sin(S(nn6(i))-S_tmp)) +omega-A*sin(S_tmp) 
   END DO

   DO i=1,natom
      S(i)=S(i)+delS(i)*delt1
   END DO

!!!===========================================================
   Zax=0. ; Zay=0. ; Zbx=0. ; Zby=0. ; WS=0.
   DO i=1,natom
         S_tmp=S(i)
         !!! Order parameter:     
         IF (Ks(i)>0.9) THEN
            Zax=Zax+cos(S_tmp) 
            Zay=Zay+sin(S_tmp)
         ELSE
            Zbx=Zbx+cos(S_tmp) 
            Zby=Zby+sin(S_tmp)
         END IF
         !!!----------------------
         !!! Wave speed:
         WS=WS+delS(i)

   END DO
   Zx=(Zax+Zbx)/n_atom ; Zy=(Zay+Zby)/n_atom
   Zax=Zax/number_a ; Zay=Zay/number_a ; Zbx=Zbx/number_b ; Zby=Zby/number_b

   Za=sqrt(Zax**2.+Zay**2.)
   Zb=sqrt(Zbx**2.+Zby**2.)
   Z_total=sqrt(Zx**2.+Zy**2.)
   WS=WS/n_atom
     
   END SUBROUTINE value_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE histogram()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE histogram()
   IMPLICIT NONE

   DO i=1,natom

      S1(i)=S(i)-int(S(i)/2./pi)*2.*pi
   
      IF (S1(i)<0.) THEN
         S1(i)=S1(i)+2.*pi
      END IF

      i_histo=int(S1(i)/del_histo)+1

      IF (Ks(i)>0.9) THEN
         H1_histo(i_histo)=H1_histo(i_histo)+1.
      ELSE
         H2_histo(i_histo)=H2_histo(i_histo)+1.
      END IF

   END DO

   END SUBROUTINE histogram

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE average_thermal()
   IMPLICIT NONE   

   Za_av=0. ; Zb_av=0. ; Z_total_av=0. ; WS_av=0.

   DO i_times=1,n_times
      DO i_loop1=1,n_equi2
         CALL equi_lattice()
      END DO

      DO i_av=1,n_average

         CALL value_thermal()
         CALL average_phase()

         !IF (mod(i_av,100)==0) THEN
         WRITE(Ligne21,*)i_av,Za,Zb,Z_total,abs(WS)
         WRITE(21,'(a)') trim(Ligne21)

         WRITE(Ligne23,*)i_av,phi_a,phi_b,delta_ab
         WRITE(23,'(a)') trim(Ligne23)

         !END IF     
                  
         Za_av=Za_av+Za
         Zb_av=Zb_av+Zb
         Z_total_av=Z_total_av+Z_total
         WS_av=WS_av+WS                     

      END DO
      
   END DO

   Za_av=Za_av/real(n_times*n_average)
   Zb_av=Zb_av/real(n_times*n_average)
   Z_total_av=Z_total_av/real(n_times*n_average)
   WS_av=WS_av/real(n_times*n_average)

   WRITE(Ligne22,*)pa1,Za_av,Zb_av,Z_total_av,abs(WS_av)
   WRITE(22,'(a)') trim(Ligne22)

   CLOSE(21)
   CLOSE(22)
   CLOSE(23)

   !!! Histogram
   del_histo = 2.*pi/(n_histo-1)
   H1_histo(:)=0. ; H2_histo(:)=0.

   WRITE(*,*)'Running histogram'
   CALL histogram()
   OPEN(unit=40,file='histogram.dat')
   DO i_histo=1,n_histo
      H1_histo(i_histo)=H1_histo(i_histo)/number_a
      H2_histo(i_histo)=H2_histo(i_histo)/number_b     
      S_histo(i_histo)=(i_histo-0.5)*del_histo
      WRITE(Ligne40,*)S_histo(i_histo),H1_histo(i_histo),H2_histo(i_histo)
      WRITE(40,'(a)') trim(Ligne40)
   ENDDO
   CLOSE(40)
      
   END SUBROUTINE average_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE index_coordinate() : transfer index i to coordinates x, y, z 
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE index_coordinate()
   IMPLICIT NONE

   DO i=1,natom
      z(i)=(i-1)/(nx*ny)+1
      i1=i-(z(i)-1)*nx*ny
      y(i)=(i1-1)/nx+1
      i2=i1-(y(i)-1)*nx
      x(i)=i2
   END DO

   END SUBROUTINE index_coordinate

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE nearest_neighbor() : find nearest neighbor j of every i
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE nearest_neighbor()
   IMPLICIT NONE

   DO i=1,natom
   DO j=1,natom
   
      IF (((x(j)==(x(i)-1)).or.(x(j)==(x(i)-1+nx))).and.(y(j)==y(i)).and.(z(j)==z(i))) THEN
         nn1(i)=j
      END IF

      IF (((x(j)==(x(i)+1)).or.(x(j)==(x(i)+1-nx))).and.(y(j)==y(i)).and.(z(j)==z(i))) THEN
         nn2(i)=j
      END IF

      IF (((y(j)==(y(i)-1)).or.(y(j)==(y(i)-1+ny))).and.(x(j)==x(i)).and.(z(j)==z(i))) THEN
         nn3(i)=j
      END IF

      IF (((y(j)==(y(i)+1)).or.(y(j)==(y(i)+1-ny))).and.(x(j)==x(i)).and.(z(j)==z(i))) THEN
         nn4(i)=j
      END IF

      IF (((z(j)==(z(i)-1)).or.(z(j)==(z(i)-1+nz))).and.(x(j)==x(i)).and.(y(j)==y(i)))  THEN
         nn5(i)=j
      END IF

      IF (((z(j)==(z(i)+1)).or.(z(j)==(z(i)+1-nz))).and.(x(j)==x(i)).and.(y(j)==y(i))) THEN
         nn6(i)=j
      END IF

   END DO
   END DO

   !!OPEN(unit=41,file='neighbor.dat')
   !!DO i=1,natom
   !   WRITE(41,*)x(i),y(i),z(i),x(nn6(i)),y(nn6(i)),z(nn6(i))
   !   
   !   IF ((x(i)==10).and.(y(i)==10).and.(z(i)==10)) THEN
   !      WRITE(*,*)x(nn4(i)),y(nn4(i)),z(nn4(i))
   !   END IF

   !!ENDDO
   !!CLOSE(41)


   END SUBROUTINE nearest_neighbor

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE average phase phi_a, phi_b, phi_total: OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE average_phase()
   IMPLICIT NONE
     
   !!! Tinh phi_a
   IF (Zay >=0.) THEN      
      phi_a=acos(Zax/Za)
   ELSE
      phi_a=2.*pi-acos(Zax/Za)
   END IF

   !!! Tinh phi_b
   IF (Zby >=0.) THEN      
      phi_b=acos(Zbx/Zb)
   ELSE
      phi_b=2.*pi-acos(Zbx/Zb)
   END IF

   !!! Tinh phi_total
   IF (Zy >=0.) THEN
      phi_total=acos(Zx/Z_total)
   ELSE
      phi_total=2.*pi-acos(Zx/Z_total)
   END IF

   delta_ab=phi_a-phi_b
   IF (delta_ab < 0.) THEN
      delta_ab=delta_ab+2.*pi
   END IF

   END SUBROUTINE average_phase

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   END PROGRAM

