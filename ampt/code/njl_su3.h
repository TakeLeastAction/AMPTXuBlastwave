!****************Dirac Distribution********************
      subroutine f(x,phi,t,cp,y)
!******************************************************
      implicit none
      real*8, intent(in) :: x,phi,t,cp
      real*8, intent(out) :: y
      real*8 :: a, b
      a=Exp((x-cp)/t)
      b=(1.0+2.0*Cos(phi))/3.0
      if (a.lt.10000000000.0) then 
      y=(b*(a*a+2.0*a)+1)/(a*a*a+3.0*b*(a*a+a)+1)
      else 
      y=0.0
      end if
      end subroutine f
!******************************************************

!****************Dirac Distribution********************
      subroutine Imf(x,phi,t,cp,y)
!******************************************************
      implicit none
      real*8, intent(in) :: x,phi,t,cp
      real*8, intent(out) :: y
      real*8 :: a
      a=Exp((x-cp)/t)
      if (a.lt.10000000000.0) then
      y=a*sin(phi)/(a*a+2.0*a*cos(phi)+1)
      else
      y=0.0
      end if
      end subroutine Imf
!******************************************************

!****************Dirac Distribution********************
      subroutine at(t,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t
      real*8, intent(out) :: y
      y=3.51-2.47*0.21/t+15.2*(0.21/t)**2.0
      end subroutine at
!******************************************************

!****************Dirac Distribution********************
      subroutine bt(t,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t
      real*8, intent(out) :: y
      y=1.75*(0.21/T)**3.0
      end subroutine bt
!******************************************************

!****************Dirac Distribution********************
      subroutine dat(t,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t
      real*8, intent(out) :: y
      y=2.47*0.21/t/t-2.0*15.2*(0.21/t)**2.0/t
      end subroutine dat
!******************************************************

!****************Dirac Distribution********************
      subroutine dbt(t,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t
      real*8, intent(out) :: y
      y=-3.0*1.75*(0.21/t)**3.0/t
      end subroutine dbt
!******************************************************


!****************Dirac Log Patition********************
      subroutine o(x,t,cp,y)
!******************************************************
      implicit none
      real*8, intent(in) :: x,t,cp
      real*8, intent(out) :: y
      y=t*Log(1.0+Exp((cp-x)/t))
      end subroutine o
!******************************************************

!****************Dirac Distribution********************
      subroutine U3(t,phi,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi
      real*8, intent(out) :: y
      real*8 :: aa, bb
      call at(t,aa)
      call bt(t,bb)
      y=1.0+cos(phi)*2.0
      y=y*(18.0*bb/sin(phi)-2.0*aa*sin(phi))
      y=-y*t*t*t/9.0
      end subroutine U3
!******************************************************

!****************Dirac Distribution********************
      subroutine U(t,phi,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi
      real*8, intent(out) :: y
      real*8 :: aa, bb
      call at(t,aa)
      call bt(t,bb)
      y=1.0+cos(phi)*2.0
      y=aa*y*y/18.0+bb*log(64.0*sin(phi/2)**4.0*sin(phi)**2.0/27.0)
      y=-y*t*t*t*t
      end subroutine U
!******************************************************

!****************Dirac Distribution********************
      subroutine tUt(t,phi,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi
      real*8, intent(out) :: y
      real*8 :: aa, bb, da, db
      call at(t,aa)
      call bt(t,bb)
      call dat(t,da)
      call dbt(t,db)
      y=1.0+cos(phi)*2.0
      y=(4.0*aa+da*t)*y*y/18.0+(4.0*bb+db*t)*log(64.0*sin(phi/2)**4.0*sin(phi)**2.0/27.0)
      y=-y*t*t*t*t
      end subroutine tUt
!******************************************************

!****************Dirac Distribution********************
      subroutine dU(t,phi,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi
      real*8, intent(out) :: y
      real*8 :: aa, bb, da, db
      call at(t,aa)
      call bt(t,bb)
      call dat(t,da)
      call dbt(t,db)
      y=1.0+cos(phi)*2.0
      y=(3.0*aa+da*t)*y*y/18.0+(3.0*bb+db*t)*log(64.0*sin(phi/2)**4.0*sin(phi)**2.0/27.0)
      y=y*t*t*t*t
      end subroutine dU
!******************************************************


!****************Charge Density************************
      subroutine Iu3(t,phi,cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, cp, m, cut
      real*8, intent(out) :: y
      real*8 :: lb, rb, cute, x, xp, dx, ans, fa, fb
      integer :: i
      cute=Sqrt(cut**2.0+m**2.0)
      lb=cp-6.0*t
      rb=cp+6.0*t
      if(lb.lt.m) lb=m
      if(lb.gt.cute) lb=cute
      if(rb.lt.lb) rb=lb
      if(rb.gt.cute) rb=cute
      ans=0.0
      if(lb.gt.m) then
         dx=(lb-m)/990.0
         do i=1,989,1
           x=m+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call Imf(x,phi,t,cp,fa)
           call Imf(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*xp*dx
         end do
         x=lb
         xp=Sqrt(x**2.0-m**2.0)
         call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
!      print *, "1st stage", ans
      if(lb.lt.cute.and.rb.gt.lb) then
      dx=(rb-lb)/1200.0
      do i=1,1199,1
         x=lb+i*dx
         xp=Sqrt(x**2.0-m**2.0)
          call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*x*xp*dx
      end do
      x=lb
      xp=Sqrt(x**2.0-m**2.0)
          call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      x=rb
      xp=Sqrt(x**2.0-m**2.0)
          call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
!      print *, "2nd stage", ans
      if(rb.lt.cute) then
         dx=(cute-rb)/990.0
         do i=1,989,1
           x=rb+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call Imf(x,phi,t,cp,fa)
           call Imf(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*xp*dx
         end do
         x=rb
         xp=Sqrt(x**2.0-m**2.0)
         call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
         x=cute
         xp=Sqrt(x**2.0-m**2.0)
         call Imf(x,phi,t,cp,fa)
         call Imf(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
!      print *, "3rd stage", ans
      y=ans*2.0/(Acos(-1.0)**2.0)
      end subroutine Iu3
!******************************************************
      
!****************Charge Density************************
      subroutine rho(t,phi,cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, cp, m, cut
      real*8, intent(out) :: y
      real*8 :: lb, rb, cute, x, xp, dx, ans,fa,fb
      integer :: i
      cute=Sqrt(cut**2.0+m**2.0)
      lb=cp-6.0*t
      rb=cp+6.0*t
      if(lb.lt.m) lb=m
      if(lb.gt.cute) lb=cute
      if(rb.lt.lb) rb=lb
      if(rb.gt.cute) rb=cute
      ans=0.0
      if(lb.gt.m) then
         dx=(lb-m)/990.0
         do i=1,989,1
           x=m+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa-fb)*x*xp*dx
         end do
         x=lb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa-fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "1st stage", ans
      if(lb.lt.cute.and.rb.gt.lb) then
      dx=(rb-lb)/1200.0
      do i=1,1199,1
         x=lb+i*dx
         xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa-fb)*x*xp*dx
      end do
      x=lb
      xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa-fb)*x*xp*dx
      x=rb
      xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa-fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "2nd stage", ans
      if(rb.lt.cute) then
         dx=(cute-rb)/990.0
         do i=1,989,1
           x=rb+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa-fb)*x*xp*dx
         end do
         x=rb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa-fb)*x*xp*dx
         x=cute
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa-fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "3rd stage", ans
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rho
!******************************************************

!****************Charge Density************************
      subroutine rhoa(t,phi,cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, cp, m, cut
      real*8, intent(out) :: y
      real*8 :: lb, rb, cute, x, xp, dx, ans,fa,fb
      integer :: i
      cute=Sqrt(cut**2.0+m**2.0)
      lb=cp-6.0*t
      rb=cp+6.0*t
      if(lb.lt.m) lb=m
      if(lb.gt.cute) lb=cute
      if(rb.lt.lb) rb=lb
      if(rb.gt.cute) rb=cute
      ans=0.0
      if(lb.gt.m) then
         dx=(lb-m)/990.0
         do i=1,989,1
           x=m+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*xp*dx
         end do
         x=lb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "1st stage", ans
      if(lb.lt.cute.and.rb.gt.lb) then
      dx=(rb-lb)/1200.0
      do i=1,1199,1
         x=lb+i*dx
         xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*x*xp*dx
      end do
      x=lb
      xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      x=rb
      xp=Sqrt(x**2.0-m**2.0)
          call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "2nd stage", ans
      if(rb.lt.cute) then
         dx=(cute-rb)/990.0
         do i=1,989,1
           x=rb+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*xp*dx
         end do
         x=rb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
         x=cute
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+0.5*(fa+fb)*x*xp*dx
      else
         ans=ans+0.0
      end if
     ! print *, "3rd stage", ans
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhoa
!******************************************************


!****************Scalar Density************************
      subroutine rhos(t,phi,cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, cp, m, cut
      real*8, intent(out) :: y
      real*8 :: lb, rb, cute, x, xp, dx, ans, fa, fb
      integer :: i
      cute=Sqrt(cut**2.0+m**2.0)
      lb=cp-6.0*t
      rb=cp+6.0*t
      if(lb.lt.m) lb=m
      if(lb.gt.cute) lb=cute
      if(rb.lt.lb) rb=lb
      if(rb.gt.cute) rb=cute
      ans=0.0
      if(lb.gt.m) then
         dx=(lb-m)/990.0
         do i=1,989,1
           x=m+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb-1)*xp*dx
         end do
         x=lb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb-1)*0.5*xp*dx
      else
         ans=ans+0.0
      end if
      if(lb.lt.cute.and.lb.lt.rb) then
      dx=(rb-lb)/1200.0
      do i=1,1199,1
         x=lb+i*dx
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb-1)*xp*dx
      end do
      x=lb
      xp=Sqrt(x**2.0-m**2.0)
      call f(x,phi,t,cp,fa)
      call f(x,phi,t,-cp,fb)
      ans=ans+(fa+fb-1)*0.5*xp*dx
      x=rb
      xp=Sqrt(x**2.0-m**2.0)
      call f(x,phi,t,cp,fa)
      call f(x,phi,t,-cp,fb)
      ans=ans+(fa+fb-1)*0.5*xp*dx
      else
         ans=ans+0.0
      end if
      if(rb.lt.cute) then
         dx=(cute-rb)/990.0
         do i=1,989,1
           x=rb+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb-1)*xp*dx
         end do
         x=rb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb-1)*0.5*xp*dx
         x=cute
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb-1)*0.5*xp*dx
      else
         ans=ans+0.0
      end if
      y=ans*3.0/(Acos(-1.0)**2.0)*m
      end subroutine rhos
!******************************************************

!****************Kinetic Energy Density****************
      subroutine rhoek(t,phi,cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi, cp, m, cut
      real*8, intent(out) :: y
      real*8 :: lb, rb, cute, x, xp, dx, ans, fa, fb
      integer :: i
      cute=Sqrt(cut**2.0+m**2.0)
      lb=cp-6.0*t
      rb=cp+6.0*t
      if(lb.lt.m) lb=m
      if(lb.gt.cute) lb=cute
      if(rb.lt.lb) rb=lb
      if(rb.gt.cute) rb=cute
      ans=0.0
      if(lb.gt.m) then
         dx=(lb-m)/990.0
         do i=1,989,1
           x=m+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*x*xp*dx
         end do
         x=lb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*0.5*x*x*xp*dx
      else
         ans=ans+0.0
      end if
      if(lb.lt.rb) then
      dx=(rb-lb)/1200.0
      do i=1,1199,1
         x=lb+i*dx
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*x*x*xp*dx
      end do
      x=lb
      xp=Sqrt(x**2.0-m**2.0)
      call f(x,phi,t,cp,fa)
      call f(x,phi,t,-cp,fb)
      ans=ans+(fa+fb)*0.5*x*x*xp*dx
      x=rb
      xp=Sqrt(x**2.0-m**2.0)
      call f(x,phi,t,cp,fa)
      call f(x,phi,t,-cp,fb)
      ans=ans+(fa+fb)*0.5*x*x*xp*dx
      else
         ans=ans+0.0
      end if
      if(rb.lt.cute) then
         dx=(cute-rb)/990.0
         do i=1,989,1
           x=rb+i*dx
           xp=Sqrt(x**2.0-m**2.0)
           call f(x,phi,t,cp,fa)
           call f(x,phi,t,-cp,fb)
           ans=ans+(fa+fb)*x*x*xp*dx
         end do
         x=rb
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*0.5*x*x*xp*dx
         x=cute
         xp=Sqrt(x**2.0-m**2.0)
         call f(x,phi,t,cp,fa)
         call f(x,phi,t,-cp,fb)
         ans=ans+(fa+fb)*0.5*x*x*xp*dx
      else
         ans=ans+0.0
      end if
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhoek
!******************************************************
      
      subroutine rho0(pf,cut,n0)
      implicit none
      real*8, intent(in) :: pf, cut
      real*8, intent(out) :: n0
      real*8 :: x
      x=pf
      if (x.gt.cut) x=cut
      n0=x*x*x/(Acos(-1.0)**2.0)
      end subroutine rho0
 
!**********Scalar Density @ zero temperature***********
      subroutine rhos0(cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: cp, m, cut
      real*8, intent(out) :: y
      real*8 :: cute,x0, x, xp, dx, ans
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+m**2.0)
      x0=cp
      if (m.gt.x0) x0=m
      dx=(cute-x0)/1200.0
      do i=1,1199,1
         x=x0+i*dx
         xp=Sqrt(x**2.0-m**2.0)
         ans=ans-xp*dx
      end do
      x=x0
      xp=Sqrt(x**2.0-m**2.0)
      ans=ans-xp*dx*0.5
      x=cute
      xp=Sqrt(x**2.0-m**2.0)
      ans=ans-xp*dx*0.5
      y=ans*3.0*m/(Acos(-1.0)**2.0)
      end subroutine rhos0
!******************************************************

!**********Energy Density @ zero temperature***********
      subroutine rhoe0(cp,m,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: cp, m, cut
      real*8, intent(out) :: y
      real*8 :: cute,x0, x, xp, dx, ans
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+m**2.0)
      x0=cp
      if (m.gt.x0) x0=m
      dx=(cute-x0)/1200.0
      do i=1,1199,1
         x=x0+i*dx
         xp=Sqrt(x**2.0-m**2.0)
         ans=ans-x*x*xp*dx
      end do
      x=x0
      xp=Sqrt(x**2.0-m**2.0)
      ans=ans-x*x*xp*dx*0.5
      x=cute
      xp=Sqrt(x**2.0-m**2.0)
      ans=ans-x*x*xp*dx*0.5
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhoe0
!******************************************************

!****************Density beyond cutoff****************
      subroutine Iu3uv(t,phi,trmu,nakem,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, trmu, nakem, cut
      real*8, intent(out) :: y
      real*8 :: cute, x, xp, dx, ans, yazi,fa,fb
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+nakem**2.0)
      yazi=cute
      if(trmu.gt.cute) then
        yazi=trmu
        dx=(yazi-cute)/490.0
        do i=1,489,1
           x=cute+i*dx
           xp=Sqrt(x**2.0-nakem**2.0)
           call Imf(x,phi,t,trmu,fa)
           call Imf(x,phi,t,-trmu,fb)
           ans=ans+(fa+fb)*x*xp*dx
        end do
        x=cute
        xp=Sqrt(x**2.0-nakem**2.0)
        call Imf(x,phi,t,trmu,fa)
        call Imf(x,phi,t,-trmu,fb)
        ans=ans+(fa+fb)*0.5*x*xp*dx
        x=yazi
        xp=Sqrt(x**2.0-nakem**2.0)
        call Imf(x,phi,t,trmu,fa)
        call Imf(x,phi,t,-trmu,fb)
        ans=ans+(fa+fb)*0.5*x*xp*dx
      else
        ans=ans+0.0
      end if
      dx=t/60.0
      do i=1,1199,1
         x=yazi+i*dx
         xp=Sqrt(x**2.0-nakem**2.0)
         call Imf(x,phi,t,trmu,fa)
         call Imf(x,phi,t,-trmu,fb)
         ans=ans+(fa+fb)*x*xp*dx
      end do
      x=yazi
      xp=Sqrt(x**2.0-nakem**2.0)
      call Imf(x,phi,t,trmu,fa)
      call Imf(x,phi,t,-trmu,fb)
      ans=ans+(fa+fb)*0.5*x*xp*dx
      x=yazi+20.0*t
      xp=Sqrt(x**2.0-nakem**2.0)
      call Imf(x,phi,t,trmu,fa)
      call Imf(x,phi,t,-trmu,fb)
      ans=ans+(fa+fb)*0.5*x*xp*dx
      y=ans*2.0/(Acos(-1.0)**2.0)
      end subroutine Iu3uv
!******************************************************

      subroutine find_phi(t, cp,trmu,mu,md,ms,mu0,md0,ms0,cut,fi)
      implicit none
      real*8, intent(in) :: t, cp,trmu,mu,md,ms,mu0,md0,ms0,cut
      real*8, intent(out) :: fi
      real*8 :: fia, fib, rb, y1a,y2ua,y2da,y2sa,y2a,y3ua,y3da,y3sa,y3a,ya, &
                y1b,y2ub,y2db,y2sb,y2b,y3ub,y3db,y3sb,y3b, yb
      integer*4 :: nlt, ngt 
      fia=0.3
      rb=Acos(-1.0)*2.0/3.0+0.1
      fib=1.0
      do while(abs(fia-fib).gt.0.0009)
         ngt=1
         nlt=1
         call U3(t,fia,y1a)
         call Iu3(t,fia,cp,mu,cut,y2ua)
         call Iu3(t,fia,cp,md,cut,y2da)
         call Iu3(t,fia,cp,ms,cut,y2sa)
         y2a=y2ua+y2da+y2sa
         call Iu3uv(t,fia,trmu,mu0,cut,y3ua)
         call Iu3uv(t,fia,trmu,md0,cut,y3da)
         call Iu3uv(t,fia,trmu,ms0,cut,y3sa)
         y3a=y3ua+y3da+y3sa
         ya=y1a+y2a+y3a
         call U3(t,fib,y1b)
         call Iu3(t,fib,cp,mu,cut,y2ub)
         call Iu3(t,fib,cp,md,cut,y2db)
         call Iu3(t,fib,cp,ms,cut,y2sb)
         y2b=y2ub+y2db+y2sb
         call Iu3uv(t,fib,trmu,mu0,cut,y3ub)
         call Iu3uv(t,fib,trmu,md0,cut,y3db)
         call Iu3uv(t,fib,trmu,ms0,cut,y3sb)
         y3b=y3ub+y3db+y3sb
         yb=y1b+y2b+y3b
         fi=(yb*fia-ya*fib)/(yb-ya)
         if (fi.gt.rb) then
             ngt=ngt+1
             fi=1.0-dble(ngt)*0.001
         end if
         if (fi.lt.0.0) then
             nlt=nlt+1
             fi=0.3+dble(nlt)*0.001
         end if
         fia=fib
         fib=fi
      end do
      end subroutine find_phi

!****************Density beyond cutoff****************
      subroutine rhouv(t,phi,trmu,nakem,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, trmu, nakem, cut
      real*8, intent(out) :: y
      real*8 :: cute, x, xp, dx, ans, yazi,fa, fb
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+nakem**2.0)
      yazi=cute
      if(trmu.gt.cute) then
        yazi=trmu
        dx=(yazi-cute)/490.0
        do i=1,489,1
           x=cute+i*dx
           xp=Sqrt(x**2.0-nakem**2.0)
           call f(x,phi,t,trmu,fa)
           call f(x,phi,t,-trmu,fb)
           ans=ans+(fa-fb)*x*xp*dx
        end do
        x=cute
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+(fa-fb)*0.5*x*xp*dx
        x=yazi
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+(fa-fb)*0.5*x*xp*dx
      else
        ans=ans+0.0
      end if
      dx=t/60.0
      do i=1,1199,1
         x=yazi+i*dx
         xp=Sqrt(x**2.0-nakem**2.0)
         call f(x,phi,t,trmu,fa)
         call f(x,phi,t,-trmu,fb)
         ans=ans+(fa-fb)*x*xp*dx
      end do
      x=yazi
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+(fa-fb)*0.5*x*xp*dx
      x=yazi+20.0*t
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+(fa-fb)*0.5*x*xp*dx
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhouv
!******************************************************

!****************Density beyond cutoff****************
      subroutine rhoauv(t,phi,trmu,nakem,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t,phi, trmu, nakem, cut
      real*8, intent(out) :: y
      real*8 :: cute, x, xp, dx, ans, yazi,fa, fb
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+nakem**2.0)
      yazi=cute
      if(trmu.gt.cute) then
        yazi=trmu
        dx=(yazi-cute)/490.0
        do i=1,489,1
           x=cute+i*dx
           xp=Sqrt(x**2.0-nakem**2.0)
           call f(x,phi,t,trmu,fa)
           call f(x,phi,t,-trmu,fb)
           ans=ans+(fa+fb)*x*xp*dx
        end do
        x=cute
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+(fa+fb)*0.5*x*xp*dx
        x=yazi
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+(fa+fb)*0.5*x*xp*dx
      else
        ans=ans+0.0
      end if
      dx=t/60.0
      do i=1,1199,1
         x=yazi+i*dx
         xp=Sqrt(x**2.0-nakem**2.0)
         call f(x,phi,t,trmu,fa)
         call f(x,phi,t,-trmu,fb)
         ans=ans+(fa+fb)*x*xp*dx
      end do
      x=yazi
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+(fa+fb)*0.5*x*xp*dx
      x=yazi+20.0*t
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+(fa+fb)*0.5*x*xp*dx
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhoauv
!******************************************************


!**********Energy Density beyond cutoff****************
      subroutine rhoeuv(t,phi,trmu,nakem,cut,y)
!******************************************************
      implicit none
      real*8, intent(in) :: t, phi, trmu, nakem, cut
      real*8, intent(out) :: y
      real*8 :: cute, x, xp, dx, ans, yazi,fa, fb
      integer :: i
      ans=0.0
      cute=Sqrt(cut**2.0+nakem**2.0)
      yazi=cute
      if(trmu.gt.cute) then
        yazi=trmu
        dx=(yazi-cute)/490.0
        do i=1,489,1
           x=cute+i*dx
           xp=Sqrt(x**2.0-nakem**2.0)
           call f(x,phi,t,trmu,fa)
           call f(x,phi,t,-trmu,fb)
           ans=ans+(fa+fb)*x*x*xp*dx
        end do
        x=cute
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+0.5*(fa+fb)*x*x*xp*dx
        x=yazi
        xp=Sqrt(x**2.0-nakem**2.0)
        call f(x,phi,t,trmu,fa)
        call f(x,phi,t,-trmu,fb)
        ans=ans+0.5*(fa+fb)*x*x*xp*dx
      else
        ans=ans+0.0
      end if
      dx=t/60.0
      do i=1,1199,1
         x=yazi+i*dx
         xp=Sqrt(x**2.0-nakem**2.0)
         call f(x,phi,t,trmu,fa)
         call f(x,phi,t,-trmu,fb)
         ans=ans+(fa+fb)*x*x*xp*dx
      end do
      x=yazi
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+0.5*(fa+fb)*x*x*xp*dx
      x=yazi+20.0*t
      xp=Sqrt(x**2.0-nakem**2.0)
      call f(x,phi,t,trmu,fa)
      call f(x,phi,t,-trmu,fb)
      ans=ans+0.5*(fa+fb)*x*x*xp*dx
      y=ans*3.0/(Acos(-1.0)**2.0)
      end subroutine rhoeuv
!******************************************************

      subroutine esea(m,cut,y)
      implicit none
      real*8, intent(in) :: m, cut
      real*8, intent(out) :: y
      y=-m**4.0*asinh(cut/m)
      y=y+cut*sqrt(m*m+cut*cut)*(m*m+2.0*cut*cut)
      y=y*3.0/8.0/(Acos(-1.0)**2.0)
      end subroutine esea 


!*****Eff Chemical Potential for Given Charge Density**
      subroutine find_mu(rh,t,fi,mu,md,ms,mu0,md0,ms0,Gv,cut,cp,trmu)
!******************************************************
      implicit none
      real*8, intent(in) :: rh, t,fi, mu,md,ms, mu0,md0,ms0, Gv, cut
      real*8, intent(out) :: cp, trmu
      real*8 :: x0, xa, xb, x, y1ua,y1da,y1sa,y1a, y2ua,y2da,y2sa,y2a, ya,&
               y1ub,y1db,y1sb,y1b, y2ub,y2db,y2sb,y2b, yb,y1u,y1d,y1s,y1
      integer :: ngt, nlt
      x0=Sqrt((rh*Acos(-1.0)**2.0/2.0)**(2.0/3.0)+mu**2.0)
      xb=x0
      xa=max(0.0,xb-t)
      ngt=1
      nlt=1
      do while(Abs(xa-xb).gt.0.0001*x0)
         call rho(t,fi,xa,mu,cut,y1ua)
         call rho(t,fi,xa,md,cut,y1da)
         call rho(t,fi,xa,ms,cut,y1sa)
         y1a=y1ua+y1da+y1sa
         call rhouv(t,fi,xa+2.0*Gv*y1a,mu0,cut,y2ua)
         call rhouv(t,fi,xa+2.0*Gv*y1a,md0,cut,y2da)
         call rhouv(t,fi,xa+2.0*Gv*y1a,ms0,cut,y2sa)
         y2a=y2ua+y2da+y2sa
         ya=y1a+y2a-rh
         call rho(t,fi,xb,mu,cut,y1ub)
         call rho(t,fi,xb,md,cut,y1db)
         call rho(t,fi,xb,ms,cut,y1sb)
         y1b=y1ub+y1db+y1sb
         call rhouv(t,fi,xb+2.0*Gv*y1b,mu0,cut,y2ub)
         call rhouv(t,fi,xb+2.0*Gv*y1b,md0,cut,y2db)
         call rhouv(t,fi,xb+2.0*Gv*y1b,ms0,cut,y2sb)
         y2b=y2ub+y2db+y2sb
         yb=y1b+y2b-rh
         x=(yb*xa-ya*xb)/(yb-ya)
         if (x.gt.x0) then
            ngt=ngt+1
            x=x0-x0*dble(ngt)*0.005
         end if
         if (x.lt.0.0) then
            nlt=nlt+1
            x=x0*dble(nlt)*0.005
         end if
         xa=xb
         xb=x        
      end do
      cp=x
      call rho(t,fi,x,mu,cut,y1u)
      call rho(t,fi,x,md,cut,y1d)
      call rho(t,fi,x,ms,cut,y1s)
      y1=y1u+y1d+y1s
      trmu=x+2.0*Gv*y1
      end subroutine find_mu
!******************************************************

!***Eff Mass and Chemical Potential under mean field***
      subroutine MFE(rhq,t,mu0,md0,ms0,Gs,K,Gv,cut,mu,md,ms,fi,cp,trmu)
!******************************************************
      implicit none
      real*8, intent(in) :: rhq, t, mu0,md0,ms0, Gs,K, Gv, cut
      real*8, intent(out) :: mu, md, ms, fi, cp, trmu
      real*8 :: rh, mua, mub, mda, mdb, msa, msb, fia, fib, uu,dd,ss
      integer :: nlt, ngt
      rh=rhq*(0.197)**3.0
      mu=mu0
      mua=-1.0
      mub=mu0
      md=md0
      mda=-1.0
      mdb=md0
      ms=ms0
      msa=-1.0
      msb=ms0
      fi=0.0
      fia=fi
      fib=fi
      nlt=0
      ngt=1
      do while((Abs(mua-mub).gt.0.001*mub).or.(Abs(mda-mdb).gt.0.001*mdb).or. &
               (Abs(msa-msb).gt.0.001*msb).or.(Abs(fib-fia).gt.0.001))
         mub=mu
         mdb=md
         msb=ms
         fib=fi 
         call find_mu(rh,t,fi,mu,md,ms,mu0,md0,ms0,Gv,cut,cp,trmu)
         call rhos(t,fi,cp,mu,cut,uu)
         call rhos(t,fi,cp,md,cut,dd)
         call rhos(t,fi,cp,ms,cut,ss)
         mua=mu0-2.0*Gs*uu-2.0*K*dd*ss
         mda=md0-2.0*Gs*dd-2.0*K*uu*ss
         msa=ms0-2.0*Gs*ss-2.0*K*dd*uu
         if (mua.lt.mu0) then
            nlt=nlt+1
            mua=mu0+dble(nlt)*0.1*mu0
         end if
         if (mua.gt.0.4) then
            ngt=ngt+1
            mua=0.4-dble(ngt)*0.1*mu0
         end if
         if (mda.lt.md0) then
            nlt=nlt+1
            mda=md0+dble(nlt)*0.1*md0
         end if
         if (mda.gt.0.4) then
            ngt=ngt+1
            mda=0.4-dble(ngt)*0.1*md0
         end if
         if (msa.lt.ms0) then
            nlt=nlt+1
            msa=ms0+dble(nlt)*0.1*ms0
         end if
         if (msa.gt.0.55) then
            ngt=ngt+1
            msa=0.4-dble(ngt)*0.1*ms0
         end if
         mu=mua
         md=mda
         ms=msa
      end do
      call find_mu(rh,t,fi,mu,md,ms,mu0,md0,ms0,Gv,cut,cp,trmu)
      end subroutine MFE
!******************************************************

!***Eff Mass and Chemical Potential under mean field***
      subroutine MFEmu(trmu,t,mu0,md0,ms0,Gs,K,Gv,cut,mu,md,ms,fi,cp)
!******************************************************
      implicit none
      real*8, intent(in) :: trmu, t, mu0,md0,ms0, Gs,K, Gv, cut
      real*8, intent(out) :: mu,md,ms, fi, cp
      real*8 ::  mua,mda,msa, mub,mdb,msb, cpa, cpb, fia, fib, uu,dd,ss,vu,vd,vs
      integer :: nlt, ngt
      mu=mu0
      md=md0
      ms=ms0
      cp=trmu
      fi=0.0
      mua=-1.0
      mda=-1.0
      msa=-1.0
      cpa=-1.0
      fia=fi
      mub=mu0
      mdb=md0
      msb=ms0
      cpb=trmu
      fib=fi
      nlt=0
      ngt=1
      do while((Abs(mua-mub).gt.0.001*mub).or.(Abs(mda-mdb).gt.0.001*mdb).or. &
               (Abs(msa-msb).gt.0.001*msb).or.Abs(cpa-cpb).gt.0.001*cpb.or.Abs(fib-fia).gt.0.001)
         mub=mu
         mdb=md
         msb=ms
         cpb=cp
         fib=fi
         call rhos(t,fi,cp,mu,cut,uu)
         call rhos(t,fi,cp,md,cut,dd)
         call rhos(t,fi,cp,ms,cut,ss)
         call rho(t,fi,cp,mu,cut,vu)
         call rho(t,fi,cp,md,cut,vd)
         call rho(t,fi,cp,ms,cut,vs)
         mua=mu0-2.0*Gs*uu-2.0*K*dd*ss
         mda=md0-2.0*Gs*dd-2.0*K*uu*ss
         msa=ms0-2.0*Gs*ss-2.0*K*dd*uu
         if (mua.lt.mu0) then
            nlt=nlt+1
            mua=mu0+dble(nlt)*0.1*mu0
         end if
         if (mua.gt.0.4) then
            ngt=ngt+1
            mua=0.4-dble(ngt)*0.1*mu0
         end if
         if (mda.lt.md0) then
            nlt=nlt+1
            mda=md0+dble(nlt)*0.1*md0
         end if
         if (mda.gt.0.4) then
            ngt=ngt+1
            mda=0.4-dble(ngt)*0.1*md0
         end if
         if (msa.lt.ms0) then
            nlt=nlt+1
            msa=ms0+dble(nlt)*0.1*ms0
         end if
         if (msa.gt.0.55) then
            ngt=ngt+1
            msa=0.4-dble(ngt)*0.1*ms0
         end if
         cpa=trmu-2.0*Gv*(vu+vd+vs)
         cp=cpa
         mu=mua
         md=mda
         ms=msa
         fi=fia
      end do
      end subroutine MFEmu
!******************************************************

!***Eff Mass and Chem Potential under MF @ 0 temp******
      subroutine MFE0(trmu,mu0,md0,ms0,Gs,K,Gv,cut,mu,md,ms,fi0)
!******************************************************
      implicit none
      real*8, intent(in) :: trmu, mu0,md0,ms0, Gs,K, Gv, cut
      real*8, intent(out) :: mu,md,ms, fi0
      real*8 :: mua,mda,msa,mub,mdb,msb,uu,dd,ss,vu,vd,vs,vua,vda,vsa,vub,vdb,vsb,cp,pf
      fi0=Acos(-1.0)/1.5
      mu=mu0
      md=md0
      ms=ms0
      mua=-100.0
      mda=-100.0
      msa=-100.0
      mub=mu0
      mdb=md0
      msb=ms0
      vu=0.0
      vd=0.0
      vs=0.0
      vua=-100.0
      vda=-100.0
      vsa=-100.0
      vub=0.0
      vdb=0.0
      vsb=0.0
      do while((Abs(mua-mub).gt.0.001*mub).or.(Abs(mda-mdb).gt.0.001*mdb).or.(Abs(msa-msb).gt.0.001*msb)&
               .or.(Abs(vua-vub).gt.0.001*mub).or.(Abs(vda-vdb).gt.0.001*mdb).or.(Abs(vsa-vsb).gt.0.001*msb))
         mub=mu
         mdb=md
         msb=ms
         vub=vu
         vdb=vd
         vsb=vs
         cp=trmu-vu-vd-vs
         call rhos0(cp,mu,cut,uu)
         call rhos0(cp,md,cut,dd)
         call rhos0(cp,ms,cut,ss)
         mua=mu0-2.0*Gs*uu-2.0*K*dd*ss
         mda=md0-2.0*Gs*dd-2.0*K*uu*ss
         msa=ms0-2.0*Gs*ss-2.0*K*dd*uu
         pf=cp*cp-mu*mu
         if (pf.lt.0.0) pf=0.0
         pf=sqrt(pf)
         call rho0(pf,cut,vua)
         vua=2.0*Gv*vua
         pf=cp*cp-md*md
         if (pf.lt.0.0) pf=0.0
         pf=sqrt(pf)
         call rho0(pf,cut,vda)
         vda=2.0*Gv*vda
         pf=cp*cp-ms*ms
         if (pf.lt.0.0) pf=0.0
         pf=sqrt(pf)
         call rho0(pf,cut,vsa)
         vua=2.0*Gv*vsa
         mu=mua
         md=mda
         ms=msa
         vu=vua
         vd=vda
         vs=vsa
      end do
      end subroutine MFE0
!******************************************************

!*************Zero Point Energy in Vacuum**************
      subroutine e0(cp,mu,md,ms,cut,Gs,K,y)
!******************************************************
      implicit none
      real*8, intent(in) :: cp, mu,md,ms, cut, Gs, K
      real*8, intent(out) :: y
      real*8 :: eku,ekd,eks,uu,dd,ss
      call rhoe0(cp,mu,cut,eku)
      call rhoe0(cp,md,cut,ekd)
      call rhoe0(cp,ms,cut,eks)
      call rhos0(cp,mu,cut,uu)
      call rhos0(cp,md,cut,dd)
      call rhos0(cp,ms,cut,ss)
      y=eku+ekd+eks+Gs*(uu*uu+dd*dd+ss*ss)+4.0*K*uu*dd*ss
      end subroutine e0
!******************************************************


!******************************************************
      subroutine WriteEosI(mu0,md0,ms0,Gs,K,Gv,cut)
!******************************************************
      implicit none
      real*8, intent(in) :: mu0,md0,ms0, Gs,K,Gv, cut
      real*8 :: t, nq, trmu, mu,md,ms, fi,cp, me0,fio,ekden,eden,ek1u,&
               ek1d,ek1s,ek2u,ek2d,ek2s,zero,uu,dd,ss,field,esu,esd,ess
      open (file='pNJL_su3_Gv_0.dat', unit=302, status='unknown')
      trmu=0.0
      call MFE0(trmu,mu0,md0,ms0,Gs,K,Gv,cut,mu,md,ms,fio)
      call e0(trmu,mu,md,ms,cut,Gs,K,zero)
      print*, mu,md,ms
      do nq=0.001, 2.6, 0.008
      do t=0.001, 0.4, 0.002
      call MFE(nq,t,mu0,md0,ms0,Gs,K,Gv,cut,mu,md,ms,fi,cp,trmu)
      call rhos(t,fi,cp,mu,cut,uu)
      call rhos(t,fi,cp,md,cut,dd)
      call rhos(t,fi,cp,ms,cut,ss)
      call dU(t,fi,field)
      call rhoek(t,fi,cp,mu,cut,ek1u)
      call rhoek(t,fi,cp,md,cut,ek1d)
      call rhoek(t,fi,cp,ms,cut,ek1s)
      call rhoeuv(t,fi,trmu,mu0,cut,ek2u)
      call rhoeuv(t,fi,trmu,md0,cut,ek2d)
      call rhoeuv(t,fi,trmu,ms0,cut,ek2s)
      call esea(mu,cut,esu)
      call esea(md,cut,esd)
      call esea(ms,cut,ess)
      ekden=ek1u+ek1d+ek1s+ek2u+ek2d+ek2s
      eden=ekden+Gs*(uu*uu+dd*dd+ss*ss)+4.0*K*uu*dd*ss-zero-esu-esd-ess
      ekden=ekden/(0.197**3.0)
      eden=eden/(0.197**3.0)
      fi=(1.0+2.0*Cos(fi))/3.0
      write(302,1111) nq,ekden, mu,md,ms,fi,t,cp,eden 
      write(*,1111) nq,t,ekden,eden, mu,md,ms,fi 
      end do
      end do
1111  format(9(1x,f10.6))
      end subroutine WriteEosI
!******************************************************

