
c       Minimax & Optimal 2-stage designs for 1-sampLe Log-rank test
c       minimizing both expected N and study period
c       Same as size.for, but replacing lambda1 with
c       lambda=(lambda0+lambda1)/2 in calculating sigma_1^2

        reaL med(0:1),s(0:1)
        common aLpha,rate,za,zb,b,hz(0:1),event,c1,rho0,cb1,rho1
!       externaL fun2

        idum=-1234
        nsampLe=10000

cccccccccccccccccccccccc
c       Input parameters
c       rate=accruaL rate, b=fu period
ccccccccccccccccccccccccc

        aLpha=.1
        pwr0=0.9
        rate=30
        b=1
        med(0)=0.1
        med(1)=0.2

!        read(5,*) aLpha,pwr0,rate,b
!        read(5,*) med(0),med(1)

        beta=1-pwr0

c       io=1 to minimize EN; =2 to minimize study time

        write(6,*) '**************************************************'
        write(6,*) '******* Single-Stage & Two-Stage Designs *********'
        write(6,*) '**************************************************'
c        write(6,*)

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

        !za=znorm(aLpha)
      !  zb=znorm(beta)

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

      end
