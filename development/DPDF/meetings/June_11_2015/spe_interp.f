c       this interpolates group of 3 adjacent points in ARCS spe file
c       for other stuff has to be adapted
	INTEGER*4 n_phi,n_e
	INTEGER*4 new_phi
	INTEGER*4 indata_len
	INTEGER*4 n_step(10),n_gap(10)
	REAL*4 s(500,1000),e(500,1000),s1(500,1000),e1(500,1000)
	REAL*4 phi_grid(500),e_grid(1000)
	REAL*4 phi_grid1(500)
	real*4 d_low(10),d_high(10),d_step(10)
        DIMENSION A(20),B(20),BC(20),T(200)
        CHARACTER*80 filenam, det_file
        write(6,'(a,$)')' Enter SpE filename :'
        read(5,'(a)') filenam
	OPEN (UNIT=20, FILE=filenam, STATUS='OLD')
c	IF (err .GT. 0 ) THEN
c	   WRITE(6,'(1X,A)') 'ERROR - opening the data file'
c	   STOP
c	END IF 
	WRITE(6,'(1X,A)') 'Reading data from the file :: '//filenam
	READ(20, '(2I5)') n_phi,n_e
	write(6,*)n_phi,n_e
	READ(20, '(a)')line
	write(6,'(a)') line
	write(6,*)n_phi,n_e
	READ(20, '(8g12.5)',err=101) (phi_grid(i),i=1,n_phi+1)
	READ(20, '(A)') line
	write(6,'(a)') line
	READ(20, '(8g12.5)',err=102) (e_grid(i),i=1,n_e+1)
c the file format adds zero at the end of phi_grid and
c there is one more point in energy grid so that the correct energy 
c is averege between n and n+1
	DO i=1,n_phi
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8g12.4)') (s(i,j),j=1,n_e)
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8g12.4)') (e(i,j),j=1,n_e)
	   DO j = 1,n_e
	      e(i,j) = e(i,j)**2
	   END DO
	END DO                 
	CLOSE (20)

c input detector angles be interpolated, they have to be 
c ascending and and 3 in a row are bad                         
	write(6,'(a,$)')' mdet angles file::'
	read(5,*) det_file
	open (unit=3,file=det_file,status='old')
	read(3,*)n_det
	do i=1,n_det
	   read(3,*)d_low(i)
	end do
	close(3)
	

c find boundary in phi grid

	do i=1,n_det
	   do 55 j=1,n_phi
	      if(abs((d_low(i)-phi_grid(j))).lt.0.005) then
		 n_gap(i)=j
		 goto 56
	      end if
 55	   continue
 56	   continue
	end do

	      
	do i=1,n_det
	   write (6,*) n_gap(i), phi_grid(n_gap(i))
	end do
c	stop 'for now'




c     loop at given energy over phi 
c im=degree of the polynomial
	IM=2
	NG=3
	NJ=7
	do 58 m=1,n_e
	   do i=1,n_det
c 
	      IX=n_gap(i)

	      do J=1,NG
		 A(j)=phi_grid(IX-NG+J-1)
		 B(j)=s(IX-NG+J-1,m)
	      end do

	      do J=NG+1,NJ
		 A(j)=phi_grid(IX-NG+J+2)
		 B(j)=s(IX-NG+J+2,m)
	      end do
	      
	      QX=phi_grid(IX)
	      s(IX,m)=AKNINT(QX,NJ,IM,A,B,T)
	      
	      do J=1,NG
		 A(j)=phi_grid(IX-NG+J-1)
		 B(j)=e(IX-NG+J-1,m)
	      end do

	      do J=NG+1,NJ
		 A(j)=phi_grid(IX-NG+J+2)
		 B(j)=e(IX-NG+J+2,m)
	      end do
	      
	      QX=phi_grid(IX)
	      e(IX,m)=AKNINT(QX,NJ,IM,A,B,T)

	   end do
 58	continue
c       fills end            

        do i=1,n_phi
           do j=1,n_e
              if((s(i,j).lt.0.0).or.(e(i,j).lt.0.0)) then
                 write (6,'(a)')' Phi and energy and s and err'
                 write(6,*) phi_grid(i)
     +  ,e_grid(j),s(i,j),e(i,j)
		if(s(i,j).lt.0.0) then
		 s(i,j)=0.0 
		 e(i,j)=0.0
		end if
		if(e(i,j).lt.0.0) then
		 e(i,j)=0.0
		end if
              end if
           end do
        end do





	WRITE(6,'(1X,A)') 'Data from the file :: '//filenam
        write(6,'(a,$)')' Enter output SpE filename :'

        read(5,'(a)') filenam
	OPEN (UNIT=21, FILE=filenam, STATUS='new')
c	IF (err .GT. 0 ) THEN
c	   WRITE(6,'(1X,A)') 'ERROR - creating the data file'
c	   STOP
c	END IF 
	WRITE(6,'(1X,A)') 'Writing data to the file :: '//filenam
	write(21, '(2I5)') n_phi,n_e
	write(21, '(a)')'### Phi Grid'
	write(21, '(8g12.5)',err=101) (phi_grid(i),i=1,n_phi+1)
	write(21, '(A)')'### Energy Grid'
	write(6,'(a)') line
	write(21, '(8g12.5)',err=102) (e_grid(i),i=1,n_e+1)
	DO i=1,n_phi
	   write(21,'(a)') '### S(Phi,w)' 
	   write(21, '(8g12.4)') (s(i,j),j=1,n_e)
	   write(21, '(a)') '### Errors' 
	   write(21, '(8g12.4)') (sqrt(e(i,j)),j=1,n_e)
	END DO                 
	CLOSE (20)
	stop
 101	write(6,'(a,2x,i5,f8.4)') 'Error in phi grid  ',i-1,phi_grid(i-1)
        stop
 102    write(6,'(a,2x,i5,f8.4)') 'Error in e grid  ',i-1,e_grid(i-1)
        stop                 
        end
                 






      FUNCTION AKNINT (XBAR,IN,IM,X,Y,T)
C
C
C AITKEN REPEATED INTERPOLATION
C
C   XBAR = ABSCISSA AT WHICH INTERPOLATION IS DESIRED
C   IABS(IN) = NO. OF VALUES IN TABLE
C              IF IN.GT.0, CHK ORDERING OF X(I).
C              IF IN.LT.0, SKIP PRECEEDING TEST.
C   IM   = DEGREE OF APPROXIMATING POLYNOMIAL
C   X    = VECTOR OF IABS(IN) VALUES OF ABSCISSA
C   Y    = VECTOR OF IABS(IN) VALUES OF ORDINATE
C   T    = TEMPORARY STORAGE VECTOR OF 4*(M+1) LOCATIONS)
C
C
C
      DOUBLE PRECISION T, DXBAR
      DIMENSION T(200), X(*), Y(*)
      DXBAR=XBAR
      N=IABS(IN)
      M=IM
      IF (M.GE.N) GO TO 120
   10 K=N-1
      IF (N.LT.2) GO TO 110
      S=X(2)-X(1)
      IF (IN.LT.0) GO TO 30
C CHK IF ORDER MONOTONIC
      IF (N.EQ.2) GO TO 30
      DO 20 I=3,N
      Z=(X(I)-X(I-1))*S
   20 IF (Z.LE.0.) GO TO 100
   30 IF (S.LT.0.) GO TO 50
C INCREASING ORDER
      DO 40 J=1,N
   40 IF (XBAR.LE.X(J)) GO TO 70
      J=N
      GO TO 70
C DECREASING ORDER
   50 DO 60 J=1,N
   60 IF (XBAR.GE.X(J)) GO TO 70
      J=N
   70 K=M
      M=M+1
      J=J-M/2
      J=MAX0(J,1)
      J=MIN0(J,N-K)
      MEND=J+K
      DO 80 I=J,MEND
      KK=I-J+1
      T(KK)=Y(I)
   80 T(KK+M)=X(I)-DXBAR
      DO 90 I=1,K
      KK=I+1
      DO 90 JJ=KK,M
      T(JJ)=(T(I)*T(JJ+M)-T(JJ)*T(I+M))/(X(JJ+J-1)-X(I+J-1))
   90 CONTINUE
      AKNINT=T(M)
      RETURN
  100 PRINT 4
    4 FORMAT (35H AKNINT X(I) NOT SEQUENCED PROPERLY)
  110 PRINT 3
    3 FORMAT (36H AKNINT N.LT.2 YBAR RETURNED AS Y(1))
      AKNINT=Y(1)
      RETURN
  120 PRINT 1
    1 FORMAT (48H AKNINT WARNING ORDER OF INTERPOLATION TOO LARGE)
      M=N-1
      GO TO 10
      END
C
C
      SUBROUTINE SORT (N,A,B)
      DIMENSION A(*),B(*)
      M=N-1
      DO 200 I=1,M
      I1=I+1
      DO 100 J=I1,N
      IF(A(J).GT.A(I))GO TO 100
      X=A(J)
      Y=A(I)
      A(I)=X
      A(J)=Y
      X=B(J)
      Y=B(I)
      B(I)=X
      B(J)=Y
  100 CONTINUE
  200 CONTINUE
      RETURN
      END






