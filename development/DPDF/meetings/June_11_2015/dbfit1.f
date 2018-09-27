      PARAMETER (MA=3,MFIT=3,NCA=3,Nsz=3000)
      DIMENSION X(Nsz),Y(Nsz),SIG(Nsz),A(MA),LISTA(MA),dyda(ma)
      CHARACTER ftitle*76,fhistry*80 ,FILENAME*60, line1*76, fname*80
      DIMENSION COVAR(NCA,NCA),ALPHA(NCA,NCA)
      DIMENSION X1(Nsz),Y1(Nsz),SIG1(Nsz)      
      character*80 arg(10)
      character*80 fplot
      EXTERNAL FJC
      fplot='~/bin/fplot'
      n_fit=1
 1    write(6,*)' dbfit f_name Qmin Qmax p1 p2 p3'
      n_g=iargc()
      if(n_g.lt.6) then 
         stop 'Missing argument'
      endif
      do i=1,n_g
         call getarg(i,arg(i))
      end do
      read(arg(1),'(a)') filename
      write(6,'(a)') filename
      read(arg(2),*) Qmin
      read(arg(3),*) Qmax
      do i=1,ma
         read(arg(i+3),*)a(i)
      enddo
D      write(6,*) 'Qmin Qmax', Qmin,Qmax
      do i=1,ma
         write(6,'(1x,a9,i2,2x,a2,f)')'Parameter',i,':',a(i)
      end do
c      stop 'test'

      OPEN(UNIT=1,file=FILENAME,status='OLD',FORM='FORMATTED',
     1ERR=1998)
      read(1,'(a)')ftitle
      read(1,'(a)')fhistry
      read(1,*)
	  do jj=1,Nsz
           READ(1,*,end=601) x1(jj),Y1(JJ)
           sig1(jj)=sqrt(abs(y1(jj)))
D           write(6,*) x1(jj),Y1(JJ),sig1(JJ)
	  end do
601	npts=jj-1
	close(1)

        call findQ(x1,1,npts,Qmin,m1)
        call findQ(x1,1,npts,Qmax,m2)
        ndata=m2-m1+1
        do i=1,ndata
           x(i)=x1(i+m1-1)
           y(i)=y1(i+m1-1)
           sig(i)=1.0
D           write(6,*) x(i),Y(i),sig(i)
        enddo

	write(6,*)' Qmin and Qmax for transformation:',x(1),x(ndata)         

      DO I=1,MA
         LISTA(I)=I
      ENDDO
      PPAR=0.00000000001
      ALAMDA=-1.
      CALL MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     &            COVAR,ALPHA,NCA,CHISQ,FJC,ALAMDA)
      DO I=1,1000
         CALL MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     &               COVAR,ALPHA,NCA,CHISQ,FJC,ALAMDA)
         PAR=2.*ABS(CHISQ-OCHISQ)/(CHISQ+OCHISQ)
         IF (PAR.LT.PPAR) GO TO 100
         OCHISQ=CHISQ
      ENDDO
 100  ALAMDA=0
      CALL MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,
     &            COVAR,ALPHA,NCA,CHISQ,FJC,ALAMDA)
      write(6,'(a,i)')' Number of iterations:',i-1
      DO i=1,MA
	write(6,*)' parameter',i,a(i)
      ENDDO
      write(6,'(a,f)')' Chisquare=',chisq
      write(fhistry,'(i2,3(1x,e14.6))')n_fit,(a(i),i=1,ma)
c      write(2,'(i5,3(2x,f6.2),2x,i2,2x,i2)') ndata,0.5,1.0,12.0,1,1
      do i=1,ndata
	 call fjc(x(i),a,yfit,dyda,ma)
	 y1(i)=yfit
      end do
      call fmax(1,ndata,yl,yh,x,y)
      write(line1,'(1x,i5,2(2x,e14.6))')ndata,yl,yh
      open(unit=12,file='fit.plot')
      write(12,'(a)') line1
      do i=1,ndata
         write(12,*)x(i),y(i),y1(i)
      end do
      close(12)
      write(fname,'(a,a4)')filename(:lench(filename)),'FIT1'
      write(6,*)fname
      open(unit=12,file=fname)
      write(12,'(a)') fhistry
      write(12,'(a)') ftitle
      write(12,'(i5,3(2x,f6.2),2x,i2,2x,i2)') ndata,0.5,Qmin,Qmax,1,0
      do i=1,ndata
         write(12,*) x(i),y(i)
      enddo
      close(12)
      j=system(fplot)
      stop
 1998 write(6,'(2a)')' Error open file',filename
      END
C-------------------------------------------------------------
        SUBROUTINE FJC(X,AC,YFIT,DYDA,MA)
        REAL AC(*),DYDA(*)
	x2=x*x
        YFIT=AC(1)+ac(2)*exp(-ac(3)*x2)
        DYDA(1)=1.
        DYDA(2)=exp(-ac(3)*x2)
	DYDA(3)=-ac(2)*exp(-ac(3)*x2)
        RETURN
        END

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
      SUBROUTINE FMAX(imin,imax,ymin,ymax,X,Y)
      DIMENSION X(*),Y(*)
      XMIN=X(imin)
      YMIN=Y(imin)
      XMAX=X(imin)
      YMAX=Y(imin)
      DO 10 I=imin,imax
         XMAX=AMAX1(XMAX,X(I))
         XMIN=AMIN1(XMIN,X(I))
         YMAX=AMAX1(YMAX,Y(I))
 10      YMIN=AMIN1(YMIN,Y(I))
         RETURN
         END
      
