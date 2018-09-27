c fourie transform (S-DW factor) for elastic or total part
c	implicit none 
	integer N_soq, NRPTS,NR
	real rmax
	parameter(N_soq=2000,NRPTS=1000,RMAX=10.0,NR=1000)
        REAL*4 D(N_soq),SS(N_soq),DSS(N_soq)
	REAL*4 RR(NR),rho(NR),drho(NR), ac(4)
        CHARACTER ftitle*80,fhistry*80,f_in*80,f_out*80,fname*80
	CHARACTER*80 arg(10)
	n_fit=-1
	PI = 4.0*ATAN(1.0D0)
	DO 400 I=1,NRPTS
	  RR(I)=(RMAX/NRPTS)*I
400	CONTINUE

c	write(6,*)' rdf-b4 f_in f_out Qmin Qmax Qdump'
	n_g=iargc()
	if(n_g.lt.5) then 
         write(6,*) 'Arguments: f_in f_out Qmin Qmax Qdump'
         stop 
	endif	
	do i=1,n_g
	   call getarg(i,arg(i))
	enddo
	read(arg(1),'(a)') f_in
	write(6,'(a)') f_in
	read(arg(2),'(a)') f_out
	write(6,'(a)') f_out
	read(arg(3),*) Qnmin
	write(6,*) Qnmin
	read(arg(4),*) Qnmax
	write(6,*) Qnmax
	read(arg(5),*) Qdump
	write(6,*) Qdump
	

	
	OPEN(UNIT=1,file=f_in,status='OLD',ERR=1998)
	READ(1,'(a)') FTITLE 
	write(6,*) ftitle
	read(1,'(a)') fhistry
	write(6,*) fhistry 
	read(ftitle,*) n_fit,(ac(i),i=1,4)
	if(n_fit.lt.0) stop ' wrong format'
	read(1,*)
	write(6,*) (ac(i),i=1,4)
	  do jj=1,N_soq
           READ(1,*,end=601) d(jj),ss(JJ)
	  end do
601	nh=jj-1
	nl=1
	close(1)

	
        m1=0
        m2=0
        call findQ(d,nl,nh,qnmin,m1)
        call findQ(d,nl,nh,qnmax,m2)
         nl=m1
         nh=m2
	write(6,*)' Qmin and Qmax for transformation:',d(nl),d(nh)         


	write(fname,'(a,a4)') trim(adjustl(f_out)),'.soq'
	open(unit=4,file=fname)
	write(4,'(a)') ftitle
	write(4,'(a)') fhistry
	write(4,'(1x,i5,3(1x,f5.2),2(2xi2))') nh-nl+1,.2,.2,.2,2,0
	sum=0.0
	do i=nl,nh
	   sub=fgas(n_fit,ac,d(i))
	   ss(i)=(ss(i)-sub)
	   sum=sum+ss(i)
	   write(4,*) d(i),ss(i)+1,ss(i)
	end do
	sum=sum/(nh-nl+1)
	write(6,'(a,f10.6)')' Average= ',sum
	close(4)
c	stop 'temp end'

c make damping
	do i=nl,nh
	   smod=1.0
	   if(d(i).le.Qdump) goto 399
	   A=d(nh)-Qdump
	   Q0=d(i)-Qdump
	   SMOD=SIN(PI*Q0/A)/(PI*Q0/A)
 399	   SS(i)=smod*SS(i)
	end do
C


	AFACT=2.0/PI
	DO I=1,NRPTS
	  FS = 0.0
	  RP = RR(I)
	  DO N=nl+1,nh
	    DELD = D(N) - D(N-1)
	    dX=0.5*(D(N)+D(N-1))
	    SINUS = SIN(dX*RP) * DELD*dX
	    FS = FS + SINUS * 0.5*(SS(N)+SS(N-1))
	  ENDDO
	  afor=afact/(4.*pi*RP)
	  rho(I) = FS*afor
	ENDDO
	write(fname,'(a,a4)') trim(adjustl(f_out)),'.rdf'
     
        open(unit=2,file=fname)
	write(2,'(a)') ftitle
	write(2,'(a)') fhistry
	write(2,*) nrpts,' 0.5 0.5 0.5 1 0'
	do i=1,nrpts
	 write(2,*)rr(i),rho(i)
	end do
	close(2)
	stop
1998	write(6,'(2a)')' Error open file',filename
	END

	function fgas(n_fit,ac,x)
	real*4 ac(*)
	x2=x*x
	if(n_fit.eq.2.or.n_fit.eq.3) p=ac(4)*x2
	if(n_fit.eq.3)fgas=ac(1)+(ac(2)*x2+ac(3)*x2*x2)*exp(-p)
	if(n_fit.eq.2)fgas=ac(1)+(ac(2)+ac(3)*x2)*exp(-p)
	if(n_fit.eq.1)fgas=ac(1)+ac(2)+ac(3)*x2*exp(-ac(3)*x2)
	if(n_fit.eq.4) then
	   x1=(x-ac(3))/ac(4)
	   x2=0.5*x1*x1
	   fgas=ac(1)+ac(2)*exp(-x2)
	   return
	endif

c	if(n_fit.eq.4) then
c	   p=(x-ac(3))/ac(4)
c	   fgas=ac(2)+(ac(1)-ac(2))/(1.+exp(p))
c	   return
c	endif

	if(n_fit.eq.0) fgas=ac(1)+ac(2)*x+ac(3)*x2
	return
	end
	
	subroutine findQ(Q,nmin,nmax,qf,m1)
	real*4 Q(*)
        eps=0.0005
 28     eps=2.*eps
        m1=0
        do i=nmin,nmax
           if(abs(q(i)-qf).lt.eps) then
              m1=i
              goto 222
           end if
         end do
         write(6,'(a)')' Decreasing eps twice'
         goto 28
 222     continue
         return
         end
