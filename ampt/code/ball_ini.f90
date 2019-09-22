      program box_ini
      implicit none
      real*8 :: nq,nua,nda,nsa,nuvua,nuvda,nuvsa,na, Gs,Ks, Gv, cut, mu0,md0,ms0,mu,md,ms, m,tem0, cp, trmu,fi, rad,a,V,&
                eua,eda,esa,euvua,euvda,euvsa,ea
      real*8, dimension(3) ::  x, p
      integer :: i, k, npart, ispc,ixj,k1
      integer, dimension(:), allocatable :: seed	  
      open(unit=20, file='uni_ball_small.dat', status='old',access='append')
      open(unit=30, file='NJLMUL.dat', status='unknown')	
      read(30,*)ixj,npart	
      !do i=1,npart	
      call random_seed(size=k1)
       allocate(seed(1:k1))	 
       seed(:) = ixj   
      call random_seed(put=seed)
      !enddo
      !rad=5.0
      !a=0.5
      mu0=0.0055
      md0=0.0055
      ms0=0.1407
      cut=0.6023
      Gs=3.67/cut/cut!2.44
      Ks=-12.36/cut/cut/cut/cut/cut
      Gv=0.0!0.6*Gs
      tem0=0.07
      nq=1.5
      call  MFE(nq,tem0,mu0,md0,ms0,Gs,Ks,Gv,cut,mu,md,ms,fi,cp,trmu)
      call rhoa(tem0,fi,cp,mu,cut,nua)
      call rhoa(tem0,fi,cp,md,cut,nda)
      call rhoa(tem0,fi,cp,ms,cut,nsa)
      call rhoauv(tem0,fi,trmu,mu0,cut,nuvua)
      call rhoauv(tem0,fi,trmu,md0,cut,nuvda)
      call rhoauv(tem0,fi,trmu,ms0,cut,nuvsa)
      na=nua+nda+nsa+nuvua+nuvda+nuvsa
      na=na/(0.197**3.0)
      call rhoek(tem0,fi,cp,mu,cut,eua)
      call rhoek(tem0,fi,cp,md,cut,eda)
      call rhoek(tem0,fi,cp,ms,cut,esa)
      call rhoeuv(tem0,fi,trmu,mu0,cut,euvua)
      call rhoeuv(tem0,fi,trmu,md0,cut,euvda)
      call rhoeuv(tem0,fi,trmu,ms0,cut,euvsa)
      ea=eua+eda+esa+euvua+euvda+euvsa
      ea=ea/(0.197**3.0)
      print *, "rhoek",ea
	  

	  
      rad=5.0
      a=0.5
      V = npart/na 
      rad = (V/3.1415/4.0*3.0)**(1.0/3.0)
      a = rad/10.0 

      call VWS(rad,a,V)
      print *, "V", V, na
      !npart=int(na*V*(1.+Exp(-rad/a)))
      print*, mu, md, ms,(1.+2.*Cos(fi))/3.0, cp, trmu
      do k=1, 50
        if(k.eq.50)then
         write(20,99) k, 1, npart, 0.000, 1, 1, 1, 1
         endif
         do i=1, npart
            call WS_ball(mu,md,ms,mu0,md0,ms0,tem0,cp,trmu,fi,cut,rad,a,x,p,m,ispc)
            if(k.eq.50)then
            write(20,100)  p(1), p(2), p(3), ispc, m, x(1), x(2), x(3), 0.05
             endif
         end do
      end do
99  format(3(1x,i6),1x,f10.4,4(1x,i6))
100 format(3(1x,f10.3),1x,i4,1x,f10.3,4(1x,f10.2))
      end program
     
      include 'njl_su3.h'
      
!*****uniform initial condition of Dirac distribution**********
      subroutine WS_ball(mu,md,ms,mu0,md0,ms0,T,cp,trmu,fi,cut,rad,a,x,p,m,ispc)
!**************************************************************
       implicit none
        real*8, intent(in) :: mu,md,ms, T, cp,trmu,fi, cut, mu0,md0,ms0
      real*8,  intent(in) :: rad, a
      integer*4, intent(out) :: ispc
      real*8, dimension(3), intent(out) :: x, p
      real*8, intent(out) :: m
      real*8 :: pr,ek,ytot,ypu,yau,ypd,yad,yps,yas,xx,fxx,yxx,yy,zz,m0,xi,cp0,psi,eku,ekd,eks
      ytot=1.0
      ypu=0.0
      yau=0.0
      ypd=0.0
      yad=0.0
      yps=0.0
      yas=0.0
      do while (ytot.gt.(ypu+yau+ypd+yad+yps+yas))
        call random_number(pr)
        pr=pr**(1./3.)
        pr=pr*(cp+6.0*T)
        if (pr.le.cut) then
            eku=sqrt(pr*pr+mu*mu)
            call f(eku,fi,T,cp,ypu)
            call f(eku,fi,T,-cp,yau)
            ekd=sqrt(pr*pr+md*md)
            call f(ekd,fi,T,cp,ypd)
            call f(ekd,fi,T,-cp,yad)
            eks=sqrt(pr*pr+ms*ms)
            call f(eks,fi,T,cp,yps)
            call f(eks,fi,T,-cp,yas)
        else
            eku=sqrt(pr*pr+mu0*mu0)
            call f(eku,fi,T,trmu,ypu)
            call f(eku,fi,T,-trmu,yau)
            ekd=sqrt(pr*pr+md0*md0)
            call f(ekd,fi,T,trmu,ypd)
            call f(ekd,fi,T,-trmu,yad)
            eks=sqrt(pr*pr+ms0*ms0)
            call f(eks,fi,T,trmu,yps)
            call f(eks,fi,T,-trmu,yas)
        end if
        call random_number(ytot)
        ytot=ytot*3.0
        if (ytot.lt.(ypu+yau+ypd+yad+yps+yas)) then
           call random_number(xi)
           xi=(xi-0.5)*2.0
           call random_number(psi)
           psi=psi*2.0*Acos(-1.0)
           p(1)=pr*sqrt(1.0-xi*xi)*cos(psi)
           p(2)=pr*sqrt(1.0-xi*xi)*sin(psi)
           p(3)=pr*xi
           if (ytot.lt.yau) then
              ispc=-1
              if (pr.le.cut) then
                  m=mu
              else 
                  m=mu0
              end if
           else if (ytot.lt.yau+yad) then
              ispc=-2
              if (pr.le.cut) then
                  m=md
              else 
                  m=md0
              end if
           else if (ytot.lt.yau+yad+yas) then
              ispc=-3
              if (pr.le.cut) then
                  m=ms
              else 
                  m=ms0
              end if
           else if (ytot.lt.yau+yad+yas+ypu) then
              ispc=1
              if (pr.le.cut) then
                  m=mu
              else 
                  m=mu0
              end if
           else if (ytot.lt.yau+yad+yas+ypu+ypd) then
              ispc=2
              if (pr.le.cut) then
                  m=md
              else 
                  m=md0
              end if
           else 
              ispc=3
              if (pr.le.cut) then
                  m=ms
              else 
                  m=ms0
              end if
           end if
        end if
      end do
      fxx=0.0
      yxx=1.0
      do while (fxx.lt.yxx)
         call random_number(xx)
         xx=(rad+5.0*a)*xx**(1.0/3.0)
         call WS(xx,rad,a,fxx)
         call random_number(yxx)
      end do
      call random_number(yy)
      yy=yy*2.0*3.1416
      call random_number(zz)
      zz=(zz-0.5)*2.0
      x(1)=sqrt(1-zz*zz)*cos(yy)*xx
      x(2)=sqrt(1-zz*zz)*sin(yy)*xx
      x(3)=zz*xx
      end subroutine WS_ball
!**************************************************************


!*************Wood Saxon***********************
      subroutine WS(r, rad, a, f)
!*********************************************
      implicit none
      real*8, intent(in) :: r, rad, a
      real*8, intent(out) :: f
      if(a.gt.0.0) then
        f=1.0/(Exp((r-rad)/a)+1.0)
      else
        if(r.le.rad) then
           f=1.0
        else
           f=0.0
        end if
      end if
      end subroutine WS
!*********************************************

!*********WS Volumn***************************
      subroutine VWS(rad,a,V)
!*********************************************
       implicit none
      real*8, intent(in) :: rad, a
      real*8, intent(out) :: V
      real*8 :: ra, rb, raa, rbb, da, db, x, y
      integer*4 :: i
      V=0.0
      ra=(rad-5.0*a)**3.0
      rb=(rad+5.0*a)**3.0
      if (ra.gt.0.0) then
         da=ra/100.0
         do i=1, 99
            x=dble(i)*da
            x=x**(1.0/3.0)
            call WS(x,rad,a,y)
            V=V+y*da
         end do
         raa=rad-5.0*a
         call WS(raa,rad,a,y)
         V=V+0.5*y*da
      end if
      if (ra.lt.0.0) ra=0.0
      db=(rb-ra)/100.0
      do i=1,99
         x=ra+dble(i)*db
         x=x**(1.0/3.0)
         call WS(x,rad,a,y)
         V=V+y*db
      end do
      raa=ra**(1.0/3.0)
      rbb=rad+5.0*a
      call WS(raa,rad,a,y) 
      V=V+0.5*db*y
      call WS(rbb,rad,a,y)
      V=V+0.5*db*y
      V=4.0*3.1416/3.0*V
      end subroutine VWS
!*********************************************


