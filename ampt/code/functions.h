!****************Distribution******************
      subroutine boltzm(pr,m,T,f)
!**********************************************
       implicit none
        real*8, intent(in) :: pr, m, T
      real*8, intent(out) :: f
      f=pr*pr*Exp(-Sqrt(m*m+pr*pr)/T)
      end subroutine boltzm
!**********************************************

!****************Distribution******************
      subroutine Dirac(pr,m,cp,T,f)
!**********************************************
       implicit none
        real*8, intent(in) :: pr, m, cp, T
      real*8, intent(out) :: f
      f=pr*pr/(Exp((Sqrt(m*m+pr*pr)-cp)/T)+1.0)
      end subroutine Dirac
!**********************************************

!****************Distribution******************
      subroutine DiracF(pr,m,cp,T,f)
!**********************************************
       implicit none
        real*8, intent(in) :: pr, m, cp, T
      real*8, intent(out) :: f
      f=1.0/(Exp((Sqrt(m*m+pr*pr)-cp)/T)+1.0)
      end subroutine DiracF
!**********************************************

!****************Distribution******************
      subroutine DiracE(ek,cp,T,f)
!**********************************************
       implicit none
        real*8, intent(in) :: ek, cp, T
      real*8, intent(out) :: f
      f=1.0/(Exp((ek-cp)/T)+1.0)
      end subroutine DiracE
!**********************************************


!****************Distribution******************
      subroutine sph_woods(r,Rm,a,rho)
!**********************************************
       implicit none
        real*8, intent(in) :: r, Rm, a
      real*8, intent(out) :: rho
      rho=r*r/(Exp((r-Rm)/a)+1.0)
      end subroutine sph_woods
!**********************************************

!*************************************************************
	  subroutine Normal(rr,width,Gauss)
!*************************************************************
      implicit none
      real*8, intent(in) :: rr,width
      real*8, intent(out) :: Gauss
      Gauss=Exp(-0.5*rr/width/width)/(width*sqrt(2.0*3.1416))**3.0
      end subroutine Normal
!**************************************************************

!*************************************************************
	  subroutine Normal_1D(r,width,Gauss)
!*************************************************************
      implicit none
      real*8, intent(in) :: r,width
      real*8, intent(out) :: Gauss
      Gauss=Exp(-0.5*r*r/width/width)/(width*sqrt(2.0*3.1416))
      end subroutine Normal_1D
!**************************************************************

!***********************************************************************
      subroutine sigma_lepton(srts,mq,ml,crs)
!***********************************************************************
      implicit none
      real*8, intent(in) :: srts, mq, ml
      real*8, intent(out) :: crs
      crs=4.0*3.1416/(3.0*137*137*srts*srts*sqrt(1-4.0*mq*mq/srts/srts))
      crs=crs*sqrt(1-4*ml*ml/srts/srts)*(1+2.0*(mq*mq+ml*ml)/srts/srts+4.0*mq*mq*ml*ml/srts/srts/srts/srts)
      crs=crs*0.197*0.197
      end subroutine
!***********************************************************************


