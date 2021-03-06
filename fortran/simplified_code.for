
c       Minimax & Optimal 2-stage designs for 1-sampLe Log-rank test
c       minimizing both expected N and study period
c       Same as size.for, but replacing lambda1 with
c       lambda=(lambda0+lambda1)/2 in calculating sigma_1^2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Input parameters
c       rate=accruaL rate, b=fu period
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        reaL med(0:1),s(0:1)
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1
        externaL fun2

        aLpha=0.05
        pwr0=0.9
        rate=30
        b=1
        med(0)=aLog(2.)/0.7
        med(1)=aLog(2.)/0.5

!        read(5,*) aLpha,pwr0,rate,b
!        read(5,*) med(0),med(1)

        beta=1-pwr0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c MAIN SCRIPT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c PRINT INPUTS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(6,*) '**************************************************'
        write(6,*) '******* Single-Stage & Two-Stage Designs *********'
        write(6,*) '**************************************************'

        write(6,*)
        write(6,*) '****************'
        write(6,*) 'Input Parameters '
        write(6,*) '****************'
        write(6,*)
        write(6,101) alpha,1-beta
 101    format(3x,'(alpha 1-beta) = (',2f8.3,')')
        irate=rate
        write(6,102) irate,b
 102    format(3x,'accrual rate =',i5,3x,'f/u period =',f5.1)

        za=znorm(aLpha)
        zb=znorm(beta)

        hr=med(1)/med(0)
        do k=0,1
        hz(k)=aLog(2.)/med(k)
        enddo
        write(6,103) hz(0),hz(1)
 103    format(3x,'haz rate =',2f7.3)
        write(6,104) med(0),med(1)
 104    format(3x,'med =',2f7.2)
        write(6,105) hr
 105    format(3x,'haz ratio =',f6.3)


c 1. Select a* for single stage design requiring the smallest sample size
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         a1=.1
         a2=5
         f1=fun(a1)
         f2=fun(a2)

         write(6,106) f1,f2
106      format(3x,' starting search parameters (f1, f2) = (',2f8.3,')')

         if(f1*f2.gt.0.) then
         do it=1,10
         a1=a1/2.
         a2=a2*2
         f1=fun(a1)
         f2=fun(a2)
         write(6,*) a1,a2,f1,f2
         if(f1*f2.Lt.0.) go to 2
         enddo
         pause 'Something wrong - 1'
         endif

  2      continue

         do it=1,100

         a3=(a1+a2)/2.
         f3=fun(a3)
         write(6,*) a3,f3

         if(f1*f3.Lt.0.) then
         a2=a3
         f2=f3
         eLse
         a1=a3
         f1=f3
         endif

         if(abs(rate*(a2-a1)).Lt..1) go to 6

         enddo

         pause 'Something wrong - 3'

 6       continue

         a=a2
         n=rate*a+1

         write(6,*)
         write(6,*) '****************************************'
         write(6,*) '******* Output Design Parameters *******'
         write(6,*) '****************************************'
         write(6,*)
         write(6,107) n,a
107      format(' Single-stage design: n=',i4,' (a=',f6.3,')')

c 2. Calculate power over a range of values for a, tau and c1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(6,*)
        write(6,*) 'Two-stage designs'

        write(6,108)
108     format(3x,'n',4x,'(a)',2x,'n1',4x,'tau',6x,'c1',7x,'c',5x,'EA',
     c  5x,'EN',4x,'PET',4x,'pwr',5x,'D1',6x,'D')

        a0=a
        ntau=1000

        a=3.542604
        c1=-0.08
        tau=1.843984

        v1=1-(1-exp(-hz(0)*tau))/tau/hz(0)
        v=1-(1-exp(-hz(0)*a))*exp(-hz(0)*b)/a/hz(0)
        s11=1-(1-exp(-hz(1)*tau))/tau/hz(1)
        s1=1-(1-exp(-hz(1)*a))*exp(-hz(1)*b)/a/hz(1)
        s01=hr*s11
        s0=hr*s1
        w1=s11-s01
        w=s1-s0

        hz1=(hz(0)+hz(1))/2.
        s11=1-(1-exp(-hz1*tau))/tau/hz1
        s1=1-(1-exp(-hz1*a))*exp(-hz1*b)/a/hz1
        rho0=sqrt(v1/v)
        rho1=sqrt(s11/s1)

        c21=-3
        f1=type1(c21)
        c22=0
        f2=type1(c22)

        if(f1*f2.gt.0.) pause'something wrong-4'

        do it=1,100

        c23=(c21+c22)/2.
        f3=type1(c23)

        if(f1*f3.Lt.0.) then
        c22=c23
        f2=f3
        eLse
        c21=c23
        f1=f3
        endif

        if(abs(c21-c22).Lt..001.and.abs(f3).Lt.0.001) go to 8

        enddo

 8      continue

        c=c23
        write(6,*) c

        pet=1-cdf(c1)
        ea=a-pet*amax1(0.,(a-tau))

        write(6,*) c1
        write(6,*) pet, ea

        cb1=sqrt(s01/s11)*(c1-w1*sqrt(rate*tau)/sqrt(s01))
        cb=sqrt(s0/s1)*(c-w*sqrt(rate*a)/sqrt(s0))

        caLL qromb(fun2,-10.,cb,pwr)

        write(6,*) pwr

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SUPPORTING FUNCTIONS BY AUTHOR
c Supporting functions handwritten by the authors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       Calculates expected numbers of events in final analysis
c       (D in manuscript)
c       not entirely sure what return value is
        function fun(a)
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1

        v1=1-exp(-hz(1)*b)/a/hz(1)+exp(-hz(1)*(a+b))/a/hz(1)
        v0=hz(0)/hz(1)*v1
        hr=hz(0)/hz(1)
        w=v1-v0

        hz1=(hz(0)+hz(1))/2.
        v=1-exp(-hz1*b)/a/hz1+exp(-hz1*(a+b))/a/hz1

        fun=rate*a-(sqrt(v0)*za+sqrt(v)*zb)**2/w**2


        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccc

        function type1(c)
c       calculate type one error rate
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1
        externaL fun1

        caLL qromb(fun1,-10.,c,ss)
        type1=ss-aLpha

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccc

        function fun1(z)
c       function to integrate over to get type1 error rate
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1
        external cdf
        phi=2.*asin(1.)
        fun1=1/sqrt(2.*phi)*exp(-z**2/2.)
        fun1=fun1*cdf((c1-rho0*z)/sqrt(1-rho0**2))
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccc

        function fun2(z)
c       function to integrate over to get power
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1
        phi=2.*asin(1.)
        fun2=1/(2*phi)**.5*exp(-z**2/2.)
        fun2=fun2*cdf((cb1-rho1*z)/sqrt(1-rho1**2))
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ZNORM AND SUPPORTING FUNCTIONS
c Calculates quantiles of standard normal.  Equivalent to qnorm in R
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function znorm(x)
      externaL erf
      pi=2.*asin(1.)
      z0=1.
      do 2 k=1,100
      f=erf(z0/sqrt(2.))/2.-.5+x
      df=exp(-z0**2./2.)/sqrt(2.*pi)
      znorm=z0-f/df
      if(abs(znorm-z0).Lt..001.and.abs(f).Lt..001) go to 4
      z0=znorm
2      continue
      write(6,*) 'diverge'
4      continue
      return
      end


      function cdf(x)
      external erf
      cdf=.5+.5*erf(x/sqrt(2.))
      return
      end


      function erf(x)
      if(x.Lt.0.) then
      erf=-gammp(.5,x**2)
      eLse
      erf=gammp(.5,x**2)
      endif
      return
      end


      FUNCTION gammLn(xx)
      REAL gammLn,xx
      INTEGER j
      REAL*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*Log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammLn=tmp+Log(stp*ser/x)
      return
      END


      FUNCTION gammp(a,x)
      REAL a,gammp,x
CU    USES gcf,gser
      REAL gammcf,gamser,gLn
      if(x.Lt.0..or.a.Le.0.)pause 'bad arguments in gammp'
      if(x.Lt.a+1.)then
        caLL gser(gamser,a,x,gLn)
        gammp=gamser
      eLse
        caLL gcf(gammcf,a,x,gLn)
        gammp=1.-gammcf
      endif
      return
      END


      SUBROUTINE gser(gamser,a,x,gLn)
      INTEGER ITMAX
      REAL a,gamser,gLn,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammLn
      INTEGER n
      REAL ap,deL,sum,gammLn
      gLn=gammLn(a)
      if(x.Le.0.)then
        if(x.Lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      deL=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        deL=deL*x/ap
        sum=sum+deL
        if(abs(deL).Lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too Large, ITMAX too smaLL in gser'
1     gamser=sum*exp(-x+a*Log(x)-gLn)
      return
      END


      SUBROUTINE gcf(gammcf,a,x,gLn)
      INTEGER ITMAX
      REAL a,gammcf,gLn,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammLn
      INTEGER i
      REAL an,b,c,d,deL,h,gammLn
      gLn=gammLn(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).Lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).Lt.FPMIN)c=FPMIN
        d=1./d
        deL=d*c
        h=h*deL
        if(abs(deL-1.).Lt.EPS)goto 1
11    continue
      pause 'a too Large, ITMAX too smaLL in gcf'
1     gammcf=exp(-x+a*Log(x)-gLn)*h
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c QROMB AND SUPPORTING FUNCTIONS
c integrates a function func from a to b
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func

      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES poLint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        caLL trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          caLL poLint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).Le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      pause 'too many steps in qromb'
      END


      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL deL,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      eLse
        it=2**(n-2)
        tnm=it
        deL=(b-a)/tnm
        x=a+0.5*deL
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+deL
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END



      SUBROUTINE poLint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.Lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'faiLure in poLint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.Lt.n-m)then
          dy=c(ns+1)
        eLse
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c NEVER CALLED
c This function is a pseudorandom number generator. It is never called.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.Le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.Lt.0) idum=idum+IM1
          if (j.Le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.Lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.Lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.Lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
