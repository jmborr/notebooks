c       to read MARI spe files
	INTEGER*4 n_phi,n_e
	INTEGER*4 indata_len
	REAL*4 s(500,2000),e(500,2000)
	REAL*4 phi_grid(500),e_grid(2000)
        CHARACTER*80 filenam
        write(6,'(a,$)')' Enter SpE filename :'
        read(5,'(a)') filenam
	write(6,'(a,$)')' Enter norm factor::'
	OPEN (UNIT=20, FILE=filenam, STATUS='OLD')
c	IF (err .GT. 0 ) THEN
c	   WRITE(6,'(1X,A)') 'ERROR - opening the data file'
c	   STOP
c	END IF 
	WRITE(6,'(1X,A)') 'Reading data from the file :: '//filenam
	READ(20, *) n_phi,n_e
	write(6,*)n_phi,n_e
	READ(20, *) 
	READ(20, '(8f10.0)',err=101) (phi_grid(i),i=1,n_phi+1)
	READ(20, *)
	write(6,*) phi_grid(1),phi_grid(9),phi_grid(n_phi),phi_grid(n_phi+1)
	READ(20, '(8f10.0)',err=102) (e_grid(i),i=1,n_e+1)
c the file format adds zero at the end of phi_grid and
c there is one more point in energy grid so that the correct enrgy 
c is averege between n and n+1
	DO i=1,n_phi
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8(g10.0))') (s(i,j),j=1,n_e)
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8(g10.0))') (e(i,j),j=1,n_e)
	   DO j = 1,n_e
	      e(i,j) = e(i,j)**2
	   END DO
	END DO                 
	CLOSE (20)                          
c     loop at given phi over energies
        do i=1,n_phi
           do j=1,n_e
              if(s(i,j).lt.0.0) then
                 write (6,'(a)')' Phi and energy and s and err'
                 write(6,*) phi_grid(i)
     +  ,e_grid(j),s(i,j),sqrt(e(i,j))
c changed to sqrt of (err) since it was confusing 
		 s(i,j)=0.0
		 e(i,j)=0.0
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
                 
