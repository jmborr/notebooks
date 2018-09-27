	INTEGER*4 n_q,n_e
	INTEGER*4 indata_len
	REAL*4 s(1000,1000),e(1000,1000)
	REAL*4 q_grid(1000),e_grid(1000), SQ(1000)
	character*80  arg(10)
        character*80 f_in, f_out, fname, exten, instr
	logical range

1	write(6,*)' qcut f_in f_out'
	n_g=iargc()
	if(n_g.lt.2) then 
         write(6,*) 'Missing arguments: f_in f_out'
         stop 
	endif	
	do i=1,n_g
	   call getarg(i,arg(i))
	enddo
	read(arg(1),'(a)') f_in
	write(6,'(a)') f_in
	read(arg(2),'(a)') f_out
	write(6,'(a)') f_out
	

	OPEN (UNIT=20, FILE=f_in, STATUS='OLD')
c	IF (err .GT. 0 ) THEN
c	   WRITE(6,'(1X,A)') 'ERROR - opening the data file'
c	   STOP
c	END IF 
	WRITE(6,'(1X,A)') 'Reading data from the file :: '//f_in
	READ(20, '(2I5)') n_q,n_e
	write(6,*)n_q,n_e
	READ(20, '(a)')line
	write(6,'(a)') line
	write(6,*)n_q,n_e
	READ(20, '(8g12.5)',err=101) (q_grid(i),i=1,n_q+1)
	READ(20, '(A)') line
	write(6,'(a)') line
	READ(20, '(8g12.5)',err=102) (e_grid(i),i=1,n_e+1)
	
c
	DO i=1,n_q
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8g12.4)') (s(i,j),j=1,n_e)
	   READ(20, *) 
c	   READ(20, '(A)') value
	   READ(20, '(8g12.4)') (e(i,j),j=1,n_e)
c	   DO j = 1,n_e
c	      e(i,j) = e(i,j)**2
c	   END DO
	END DO 
	CLOSE (20)                          

c     loop at given phi over energies
        do i=1,n_q
	   sqsum=0.0
           do j=1,n_e
              if(s(i,j).gt.-999.0) then
		 sqsum=sqsum+s(i,j)
              end if
           end do
	   SQ(i)=sqsum
        end do
	t_sum=0.0
	d_sum=0.0
	do i=1,n_q-1
	   de=q_grid(i+1)-q_grid(i)
	   d_sum=d_sum+de
	   t_sum=t_sum+0.5*(SQ(i+1)+SQ(i))*de
	end do
	t_sum=t_sum/d_sum
	write(6,*)' Total sum=',t_sum
	write(6,'(a,2(2xf8.3))')'Emin and Emax:',e_grid(1),e_grid(n_e)
	emin=e_grid(1)
	emax=e_grid(n_e)
	
	write(fname,'(a,a4)') trim(adjustl(f_out)),'.soq'
	open(unit=11,file=fname)

	WRITE(6,'(1X,A)') 'Writing S(Q) to the file :: '//fname
	write(11, '(a20,2(2x,f8.3))') 'Integrated S(Q)for '
     +	         ,e_grid(1),e_grid(n_e)
	write(11, '(a,a6,g12.5)') trim(adjustl(f_in))
     +	         ,' Tsum=',t_sum
	write(11,'(i5,3(2x,f8.3),2(1x,i2))') n_q,0.5,q_grid(1),
     +          q_grid(n_q),1,0
	do i=1,n_q+1
	write(11, '(2g12.5)',err=101) q_grid(i),SQ(i)/t_sum
	end do
	close(11)

	do j=1,n_e
	  do i=1,n_q
              if(s(i,j).ne.-999.0) then
		SQ(i)=s(i,j)
		else
		SQ(i)=0
	      end if
          end do
	write(exten,'(f6.2)')e_grid(j)
c	write(6,*) trim(adjustl(exten))
c	write(fname,'(a,f6.2,a5)')f_out(:lench(f_out)),e_grid(j),'.qcut'
c	write(fname,'(a3,a,a5)')'Ec_',trim(adjustl(exten)),'.qcut'
	write(instr,'(a,a4)') trim(adjustl(f_out)),'_Ec_'
	write(fname,'(2a,a5)')trim(adjustl(instr)),trim(adjustl(exten)),'.qcut'
	OPEN (UNIT=21, FILE=fname)
c	write(6,*)fname

c	IF (err .GT. 0 ) THEN
c	   WRITE(6,'(1X,A)') 'ERROR - creating the data file'
c	   STOP
c	END IF 
	WRITE(6,'(1X,A)') 'Writing data to the file :: '//fname
	write(21, '(a20,f6.1)') 'Integrated S(Q)for '
     +	         ,e_grid(j)
	write(21, '(a,a11,f6.2,a6,g12.5)') trim(adjustl(f_in)),' S(Q)for E='
     +	         ,e_grid(j),' Tsum=',t_sum
	write(21,'(i5,3(2x,f8.3),2(1x,i2))') n_q,0.5,q_grid(1),
     +          q_grid(n_q),1,0
	do i=1,n_q+1
	write(21, '(2g12.5)',err=101) q_grid(i),SQ(i)/t_sum
c	write(21, '(2g12.5)',err=101) q_grid(i),SQ(i)
	END DO                 
	CLOSE (21)

	end do

	stop
 101	write(6,'(a,2x,i5,f8.4)') 'Error in q grid  ',i-1,q_grid(i-1)
        stop
 102    write(6,'(a,2x,i5,f8.4)') 'Error in e grid  ',i-1,e_grid(i-1)
        stop                 
        end
                 
