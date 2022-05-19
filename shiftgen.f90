! This is a general calculation of the force arising from an
! linearized oscillatory displacement (shift) perturbation of a
! localized structure that gives potential ENERGY phi. The potential,
! which is prescribed by two functions phigofz(z) and phigprimeofz(z),
! is either a hill or a valley and tends to zero at large distances,
! but its derivative passes through zero only once, at z=0, which is
! the only potential extremum (discounting infinity). When the extremum
! is a minimum, there are locally trapped particles, but when it is a
! maximum, reflected orbits are not locally trapped. There are
! therefore three types of orbit: Passing, Trapped, or Reflected,
! which must be treated differently.

! An orbit's equation of motion is dv/dt=-d\phi/dz, which means for
! species s of different mass, that the time is scaled differently: to
! omega_{ps}. Consequently, for a particular perturbation frequency
! omega, when there are multiple species the scaled value omegag must
! be set to omega/omega_{ps}, different for different mass
! species. For different charge sign, the hill peak psig must likewise
! be opposite. The length scale is normally the Debye length for some
! reference temperature, the default spatial extent is |zm|=15.
! But for energies near zero is sometimes automatically lengthened.

! The force contribution is the integral dz of -d\phi/dz times the
! non-adiabatic perturbed f, which is an integral over the past time
! d\tau of the orbit, integrated over velocity, using total (distant)
! density of unity, and given (linearized) for a unit perturbing
! z-shift.

! Both space and past time integrals can be expressed as time
! integrals.  However, equal intervals of neither time nor space are
! universally optimal choices. Equal intervals of space are suboptimal
! near a reflection. Equal intervals of time are suboptimal for orbits
! of passing energy close to zero (because they move slowly at large z
! which is a less important region).  The z array is z(i)= z1+
! z2/(1+2K)[2K(i/n)+(i/n)^2] when repelled, z1=zR, z2=(zm-zR), where
! zR is reflection point.  When attracted z1=0, z2=zR trapped, z2=zm
! passing.  And K=sqrt(max(0,W/psi-1)) repelled,
! K=-1-sqrt(max(0,-W/psi)) attracted.  Repelling hill psi>0, attracted
! valley psi<0.

! Version for inclusion of auxiliary coupled perturbation modes,
! and works only for sech^4 potential shape. 
! 1 gives just the other discrete mode, and 2 adds the continuum.
! Inner products should be only for Fattract (presumably electrons). 
! shift3mode version started 16 May 2022. nmd now controls nmodes.

module shiftgen
  integer, parameter :: ngz=100,nge=200,nhmax=60,nzd=200
  complex :: omegag=(1.,0.),Ftot,omegadiff,omegaonly,Vw
  real :: psig=.1,Wg,zm=29.9,Omegacg=5.,kg=0.
  real :: vshift=0.,vrshift=9999  ! The shape of ion distribution.
  real :: Tinf=1.,Tperpg=1.,Torepel=1.
! Tinf is really the reference (attracted). Torepel the repelled species tempr.
  real :: zR ! reflection position if any
  real :: rmime=1836. ! Ratio of mass of ion to mass of electron
  real, parameter :: pig=3.1415926
  complex, parameter :: sqm1g=(0.,1.)
  integer :: iwpowg=3,ippow=3,nharmonicsg,ivs,iws
  real :: kpar,f4norm
! Spatial Arrays variable
  real, dimension(-ngz:ngz) :: zg,vg,phig,phigprime,taug
! Fixed range  
  real, dimension(-nzd:nzd) :: zdent=0.,zdmid,vpsibyv,vinfbyv,phi0d
  complex, dimension(-nzd:nzd) :: CapPhid,dentpass,denttrap,CapPhidprev
  complex, dimension(-nzd:nzd) :: CapQd,dentq,auxzd,denqwint,dentqt
  complex, dimension(-nzd:nzd) :: ft4d,phipd,denstore,denqwintd
! Parallel energy arrays
  complex, dimension(0:nge) :: omegabg, Fnonresg
  complex, dimension(nge) :: Forcegp,Forcegt,Forcegr
  real, dimension(nge) :: Wgarray,Wgarrayp,Wgarrayr,vinfarrayp
  real, dimension(nge) :: vinfarrayr,tbr,tbp,Wgarrayu
!   Auxiliary Forces as a function of parallel energy/velocity.
  complex :: Fextqq,Fextqw,Fintqq,Fintqw,Fextqqwanal
! Perpendicular Harmonic force array.
  complex, dimension(-nhmax:nhmax) :: Fpg
! Total forces
  complex :: Ftotalrg,Ftotalpg,Ftotalsumg,FVwsumg
! Mode matrices
  integer, parameter :: nmdmax=3
  integer :: nmd=nmdmax ! Default
  complex, dimension(-ngz:ngz,nmdmax) :: pmds,CPmds,ftrmds
  complex, dimension(nmdmax,nmdmax) :: Fmdaccum
! Force totals Reflected, Passing, Sum, Attracted, Trapped, hill=Repelled  
  complex, dimension(nmdmax,nmdmax) :: Ftmdr,Ftmdp,Ftmdsum,Ftmda,Ftmdt,Ftmdh
  complex, dimension(nmdmax,nmdmax,0:nge) :: Fmdres,Fmdp,FmdofW,FmdpofW  
! Distribution function variables for mode dent densities 
  complex, dimension(-ngz:ngz,nmdmax-1) :: ftraux
  complex, dimension(-nzd:nzd,nmdmax-1) :: ftrauxd
! Whether to apply a correction to the trapped species
  logical :: lioncorrect=.true.,lbess=.false.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makezg(isigma)
    integer :: isigma
  real :: z0,z1,z2
! Calculate the zg mesh, and on it phig, vg, phigprime for an incoming
! orbit from z=isigma*zm (v sign -isigma) or trapped orbit from its
! isigma end.
    integer :: ncalls=0
    ncalls=ncalls+1
    z0=0.
    zmfac=1.
! Make sure we don't miss a shallow trapped orbit
1    if(Wg.lt.0..and.Wg.gt.phigofz(zmfac*isigma*zm))then
       zmfac=1.2*zmfac
       goto 1
    endif
! Find any reflection position zR. ! Put equal to 0 if W>max(phi).
    call orbitendg(Wg,z0,zmfac*isigma*zm)
    zR=z0
    ivs=-1
    if(psig.gt.0)then ! Repelling potential
       gK=sqrt(max(0.,Wg/psig-1.))
       z1=zR; z2=isigma*zm-zR
       if(Wg.lt.psig)ivs=1  ! Reflected orbit, all z are same sign.
    elseif(psig.lt.0)then !Attracted
       gK=-1.-10.*sqrt(max(0.,-Wg/psig))
       if(zR.eq.0.)then  ! Passing
!          gK=-100.  ! Uniform hack 
          z2=isigma*zm
       else
          z2=zR*(1.+1.e-5/ngz) ! Force trapped end vg to zero.
       endif
       z1=0.
    else  
       write(*,*)'ERROR makezg. psig is zero or worse',psig
       stop
   endif
    omegar=real(omegaonly)
! Set the kpar of the continuum mode.
    if(omegar.le.1.and.Omegacg.gt.omegar)then
       kpar=kg*real(omegaonly*sqrt((Omegacg**2+1-omegaonly**2)/&
            ((Omegacg**2-omegaonly**2)*(1-omegaonly**2))))
       Vw=1.+(kpar/omegaonly)**2+kg**2*Tperpg/(omegaonly**2-Omegacg**2)
       if(ncalls.eq.1)  write(*,'(a,f8.5,a,2f8.4)')' kpar=',kpar,' Vw=',Vw
    else
       write(*,*)'ERROR omegar >=1 is not allowed',omegar
       stop
    endif
! Construct the modes
    zg(0)=0.; phig(0)=psig; vg(0)=isigma*ivs*sqrt(2.*max(0.,Wg-phig(0)))
    do i=1,ngz
       zi=float(i)/ngz
       zg(i)=ivs*(z1+z2*((2.*gK+zi)*zi)/(1.+2.*gK))
       zg(-i)=ivs*zg(i)
       phig(i)=phigofz(zg(i))
       phig(-i)=phig(i)
       vg(i)=isigma*ivs*sqrt(2.*max(0.,Wg-phig(i)))
       vg(-i)=-ivs*vg(i)
       do j=1,nmd
          pmds(i,j)=mdofz(zg(i),j,kpar)
          pmds(-i,j)=ivs*pmds(i,j)
       enddo
    enddo
    f4norm=abs(16.*psig/(3.*sqrt(70.)))
    phigprime(:)=sign(1.,psig)*real(pmds(:,1)*f4norm)
  end subroutine makezg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine orbitendg(Wj,z0,zL)
    ! Find by bisection the turning point if any of the orbit whose
    ! energy is Wj lying between zL (normally negative) and z0 (=0).  
    ! That is where phi=-Wj, return z0 s.t. |z0| is just lower.
    real :: Wj,zm,z1,z0,enrgym,enrgy0,enrgy1
    nbi=25
    z0=0.
    z1=zL
    enrgy0=phigofz(z0)-Wj
    enrgy1=phigofz(z1)-Wj
    if(enrgy1*enrgy0.ge.0)then
       return             ! Improper starting z values.
    endif
    do i=1,nbi
       zm=(z1+z0)/2.
       enrgym=phigofz(zm)-Wj
! Which value to replace.
       if(sign(1.,enrgym).eq.sign(1.,enrgy0))then
          if(z0.eq.zm)exit ! should never happen
          enrgy0=enrgym
          z0=zm
       else
          if(z1.eq.zm)exit
          enrgy1=enrgym
          z1=zm
       endif
    enddo
  end subroutine orbitendg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Fdirect(Wgi,isigma,dForceg)
! Calculate the past integral of tau and force 
! vs z for orbit of energy W, starting at isigma side. Where
! for untrapped particles (W>=0), v_0=vinf=sqrt(2*W).
! For trapped particles (W<0), v_0=sqrt(2*(W-psig)).  
! The integration is done on a z-grid such that dz=v*dtau.
! isigma is the sign of zg at the start of orbit (v-sign=-isigma)
! The total needs to account for v-sign=+isigma as well. However,
! for a symmetric hole that just multiplies force by 2. 
! Force integration not using v-integration by parts
! dF/dv/xi/(omega df/dW)=(i)\int -dphi/dz|_t \int^t (-dphi/dz|_tau)
!          exp(-i omega(tau-t)) dtau |dz(t)|        where dz(t)=v(t)dt.
! But if we want the contribution per dvinf, dv/dvinf=vinf/v so that 
! dF/dvinf/xi/(omega df/dW)=(i)\int -dphi/dz|_t \int^t (-dphi/dz|_tau)
!          exp(-i omega(tau-t)) dtau |vinf dt|.
! Inner integral is CapPhi=\int^t (-dphi/dz|_tau) exp(-i omega(tau-t)) dtau 
! done as exp(-i omega(tau-t)) dtau = d[exp()]/(-i omega)
! The result needs to be multiplied by df/dWpar.
    complex :: dForceg,CPfactor,exptbb2,CPngz
    Wg=Wgi
    call makezg(isigma) ! Sets various time array values.
    vpsig=vg(-ngz)            ! Untrapped default ...
    if(Wg.lt.0.)vpsig=vg(0)   ! Trapped particle f(v) reference.
    ips=int(sign(1.,psig))
    iws=0                     ! dtau algorithm switch index (plotting only)
    taug(-ngz)=0.
    CPmds(-ngz,1:nmd)=0.
! Set the incoming CPmds for continuum to be the wave solution.
    if(Wg.ge.0)CPmds(-ngz,3)=pmds(-ngz,3)/(sqm1g*(-kpar*vg(ngz)-omegag))
    Fmdaccum=0.
    if(Wg.lt.phigofz(zm).and.psig.gt.0)return  ! Reflected Energy too small
    do i=-ngz+1,ngz
       phigp=0.5*(phigprime(i)+phigprime(i-1))
       vmean=0.5*(vg(i)+vg(i-1))
! Calculate dtau
       if(zR.ne.0.and.abs(vg(i)).le.0.5*sqrt(2.*max(Wg,Wg-psig)))then
          dtau=-((vg(i)-vg(i-1))/phigp) ! Use dtau=dv*dtau/dv
          if(ips.gt.0)iws=i               ! Reflected track
       else                               ! Use dtau=dx*dtau/dx
          dtau=((zg(i)-zg(i-1))*vmean/(vg(i-1)*vg(i)))
          if(ips.le.0..or.zR.eq.0)iws=i   ! Attracted or unreflected
       endif
       if(.not.dtau.lt.1e6)then           ! Test for NAN error
          write(*,*)'In Fdirect dtau error',Wg,' dtau=',dtau
          stop
       endif
! And past time integrals
       taug(i)=taug(i-1)+dtau
       CPfactor=exp(sqm1g*omegag*dtau) ! Current exponential
       CPmds(i,1:nmd)=CPfactor*CPmds(i-1,1:nmd)+0.5*&
            (pmds(i,1:nmd)+pmds(i-1,1:nmd))*(1.-CPfactor)*sqm1g/omegag
       do j=1,nmd
          Fmdaccum(1:nmd,j)=Fmdaccum(1:nmd,j)+&
               sqm1g*0.5*(conjg(pmds(i,1:nmd))*CPmds(i,j)+&
               conjg(pmds(i-1,1:nmd))*CPmds(i-1,j))*abs(vpsig*dtau)
       enddo
    enddo
    if(Wg.lt.0.)then     ! Trapped orbit. Add resonant term. 
       CPngz=CPmds(ngz,1)*f4norm
       exptbb2=exp(sqm1g*omegag*taug(ngz))
       vpsig=abs(vg(0))
! This form is to be divided later by (1-exptb) full resonant denominator.
! The division by the resonant denominator is done outside the routine
! because it involves complicated negotiation of the resonance to
! preserve accuracy for trapped particles.
       do j=1,nmd
          do i=1,nmd
             Fmdaccum(i,j)=Fmdaccum(i,j)*(1-exptbb2**2) &
                  +sqm1g*CPmds(ngz,i)*CPmds(ngz,j)*(1-exptbb2)*vpsig
          enddo
       enddo
    endif
! Put value into the passed argument. 
    dForceg=Fmdaccum(1,1)*f4norm**2  ! This is    <4V4>
  end subroutine Fdirect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgRepelEint(Ftotalg,isigma)
! Integrate the force over energy to obtain the full parallel distribution
! Just for positive (repelling) psig. isigma is the entering z-sign.
! For symmetric potentials and f(v), the returned Ftotalg can simply be
! doubled to give the total force since then it is symmetric in isigma.
    complex Ftotalg
    integer :: ncalls=1
    if(nmd.ne.1.and.ncalls.eq.1)then
       ncalls=ncalls+1
      write(*,*)'WARNING nmd=',nmd,'not implemented for repelling potential.'
    endif
    Emaxg=2.5*(sqrt(2*Tinf)+vshift)**2
    call FgPassingEint(Ftotalpg,isigma,Emaxg)
    call FgReflectedEint(Ftotalrg,isigma)
    Ftotalg=(Ftotalpg+Ftotalrg)*2. ! Both v-directions
    Ftotalg=Ftotalg*Tinf/Torepel   ! Correct for ion Temperature.
  end subroutine FgRepelEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FgPassingEint(Ftp,isigma,Emaxg)
    implicit none
    complex Ftp
    integer :: isigma,i
    real :: Emaxg
    real :: dvinf,fvinf,dfe,dfeperp
    real :: Wnext,dvnow,dvnext,vinfprior,vinfnow,vinfnext
    complex :: dfweight,ForceOfW
    omegadiff=omegag-omegaonly
    Ftp=0.
    Ftmdp=0.
    Fextqq=0.
    Fextqw=0.
    Fextqqwanal=0.
    dentpass=0.
    denqwint=0.
    
    vinfprior=sqrt(2.*max(psig,0.))  ! Speed of just passing at z=0.
    Wnext=max(psig,0.)+Emaxg*(1./float(nge))**ippow
    vinfnext=-isigma*sqrt(2.*Wnext)
    dvnext=abs(vinfnext-vinfprior )*2.  ! Approximate the start via dvnow*2.
    do i=1,nge  ! Passing, corrected for psig sign.
       Wgarray(i)=Wnext
       vinfnow=vinfnext
       dvnow=dvnext
       vinfarrayp(i)=vinfnow
       Wnext=max(psig,0.)+Emaxg*(min(i+1,nge)/float(nge))**ippow
       vinfnext=-isigma*sqrt(2.*Wnext)
       dvnext=abs(vinfnext-vinfnow)
       call Fdirect(Wgarray(i),isigma,ForceOfW)
       dfe=dfdWpar(vinfnow,fvinf) ! fparallel slope and value
       dfeperp=-fvinf/Tperpg
       dfweight=(omegag*dfe-omegadiff*dfeperp)
       dvinf=0.5*(dvnow+dvnext)
       Ftp=Ftp+ForceOfW*dvinf*dfweight
       Ftmdp=Ftmdp+Fmdaccum*dvinf*dfweight ! Increment mode force
       Forcegp(i)=ForceOfW*dfweight
       FmdpofW(:,:,i)=Fmdaccum*dfweight       ! Store contribution
       tbp(i)=taug(ngz)
       if(nmd.ge.3)then
! Now we must do or add the external continuum integration.
          call Fextern2(Wgarray(i),isigma,dvinf,dfweight)
       endif
       call dentadd(dfweight,dvinf) !Diagnostics
       vinfprior=vinfnow
    enddo
    Wgarrayp=Wgarray
  end subroutine FgPassingEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgReflectedEint(Ftr,isigma)
    implicit none
    complex :: Ftr,ForceOfW,FOfWprev
    complex :: dfweight,dfwprev
    integer :: isigma,i
    real :: dvinf,fvinf,dfe,dfeperp
    omegadiff=omegag-omegaonly
    do i=1,nge  ! Reflected
       Wgarray(i)=psig*(1.-(i/float(nge))**2)
       vinfarrayr(i)=-isigma*sqrt(2.*Wgarray(i))
       call Fdirect(Wgarray(i),isigma,ForceOfW)
       dfe=dfdWpar(vinfarrayr(i),fvinf) ! fparallel slope and value
       dfeperp=-fvinf/Tperpg
       dfweight=(omegag*dfe-omegadiff*dfeperp)
       if(i.eq.1)then
          dvinf=abs(vinfarrayr(i))-sqrt(2.*psig)
          Ftr=dvinf*ForceOfW*dfweight
       else
          dvinf=abs(vinfarrayr(i)-vinfarrayr(i-1))
          Ftr=Ftr+dvinf*0.5* (ForceOfW*dfweight +FOfWprev*dfwprev)
       endif
       Forcegr(i)=ForceOfW*dfweight
       FOfWprev=ForceOfW
       tbr(i)=taug(ngz)
       dfwprev=dfweight
    enddo
    Wgarrayr=Wgarray
  end subroutine FgReflectedEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgAttractEint(Ftotalg,isigma)
    complex Ftotalg
    Emaxg=6.*Tinf+vshift**2
    call FgPassingEint(Ftotalpg,isigma,Emaxg)
    call makezdent
    call qwint                         ! Find internal wave force.
    Ftmdp(3,3)=Ftmdp(3,3)-Fintqw       ! And subtract from qq
    call FgTrappedEint(Ftotalrg,-1/Tperpg,1.,isigma)
 ! Specifies dfperp, fperp, for Maxwellian perp distrib.    
    call qdenqint                      ! And densities when needed.
    Ftotalg=(Ftotalpg+Ftotalrg)*2. ! account for both v-directions 
    Ftmda=(Ftmdp+Ftmdt)*2   
  end subroutine FgAttractEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine FgTrappedEint(Ftotal,dfperpdWperp,fperp,isigma)
    ! Integrate over fe (trapped). Wj=vpsi^2/2-psi. So vpsi=sqrt(2(psi+Wj))
    ! Wj range 0 to -psi.
    implicit none
    complex :: Ftotal,exptb,cdvpsi,dob,dFdvpsig,resdenom,resdprev
    complex :: dfweight
    complex :: Fmdnr(1:nmdmax,1:nmdmax),Fmdnpr(1:nmdmax,1:nmdmax)
    real :: dfperpdWperp,fe,fperp,obi,vpsi,dvpsi,vpsiprev,Wj
    integer :: isigma,i

    omegadiff=omegag-omegaonly
    dentqt=0.
    denqwint=0.
    denttrap=0.
    Ftotal=0.
    vpsiprev=sqrt(-2.*psig)
    omegabg(0)=0.
    Fnonresg(0)=0.                !Don't add zero energy point.
    Ftmdt=0.
    Fmdnpr=0.
    resdprev=1.
    do i=1,nge       ! Trapped
       Wgarray(i)=psig*((float(i)/nge)**iwpowg)
       Wj=Wgarray(i)
       if(i.eq.nge)Wj=0.9999*psig   ! Prevent exact vpsi=0
       dfweight=(omegag*dfdWptrap(Wj,fe)*fperp-omegadiff*fe*dfperpdWperp)
       vpsi=sqrt(2.*(-psig+Wj))
       vinfarrayr(i)=vpsi ! reflected==trapped for attracting hill.
       dvpsi=vpsiprev-vpsi
! Calculate the force dFdvpsi for this vpsi and dvy element for one transit:
       call Fdirect(Wj,isigma,dFdvpsig)
       omegabg(i)=2.*pig/(2.*taug(ngz))
       ! Strictly to get dFdvpsi we need to multiply by the omega f' terms
       Fnonresg(i)=dFdvpsig*dfweight
       Fmdnr(1:nmd,1:nmd)=Fmdaccum(1:nmd,1:nmd)*dfweight
       if(.true.)then
          call pathshiftg(i,obi)
          omegabg(i)=omegabg(i)+sqm1g*obi
          dob=omegabg(i)-omegabg(i-1)
          cdvpsi=dvpsi*(1.+sqm1g*imag(dob)/real(dob))
       ! and correct for the imaginary shift of omegabg:
          Fnonresg(i)=Fnonresg(i)+sqm1g*real(Fnonresg(i)-Fnonresg(i-1))/dob*obi
          Fmdnr(1:nmd,1:nmd)=Fmdnr(1:nmd,1:nmd)+sqm1g*&
               &real(Fmdnr(1:nmd,1:nmd)-Fmdnpr(1:nmd,1:nmd))/dob*obi
       else
          obi=0.
          cdvpsi=dvpsi
       endif
       exptb=exp(sqm1g*omegag*pig/omegabg(i))   ! Half period
       if(.not.abs(exptb).lt.1.e10)exptb=1.e10  ! Enable divide by infinity
       resdenom=1.-exptb**2                        !Full period version
       if(i.eq.1)resdprev=resdenom    ! Why this? To avoid infinity.
       if(taug(ngz)*imag(omegag).lt.-2.)then!Hack fix giant dFdvpsi problem
          write(*,*)'Drop giant taug*omegai',taug(ngz),imag(omegag),dFdvpsig
          Fnonresg(i)=0.
          Fmdnr=0.
       endif

! Then divide by the resonance denominator and sum 
       Forcegr(i)=0.5*(Fnonresg(i)/resdenom + Fnonresg(i-1)/resdprev)
    ! Now Forcegr when multiplied by cdvpsi and summed over all ne energies
    ! and multiplied by 2 gives the total force. 
       Fmdres(1:nmd,1:nmd,i)=0.5*(Fmdnr(1:nmd,1:nmd)/resdenom &
            + Fmdnpr(1:nmd,1:nmd)/resdprev)
       if(.not.(abs(Forcegr(i)).ge.0))then
          write(*,*)'Forcegr NAN?',i
          write(*,*)Fnonresg(i),Forcegr(i)
          stop
       endif
! Add to Ftotal integral (doubled later).
       Ftotal=Ftotal+Forcegr(i)*cdvpsi
       Ftmdt=Ftmdt+Fmdres(1:nmd,1:nmd,i)*cdvpsi
! Save previous values       
       vpsiprev=vpsi
       Fmdnpr=Fmdnr
       resdprev=resdenom
       tbr(i)=taug(ngz)
       call dentaddtrap(dfweight,cdvpsi)
    enddo
    Wgarrayr=Wgarray
  end subroutine FgTrappedEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pathshiftg(i,obi)
! Calculate the required shift of the omegabg path below the real axis
! for this ith energy mesh value, assuming the previous omegabg(i-1)
! is known.
    real :: obi
    integer :: el,ielsign
    real :: doel,doel2,dob
!    real,parameter :: Sc=4,Rc=2./Sc
    real,parameter :: Sc=2,Rc=2./Sc   ! Slightly better.

    ! Select the closest odd resonance such that el*omegabg=omegar
    ielsign=int(sign(1.,real(omegag)))
    el=(2*int( (abs(real(omegag))/real(omegabg(i))-1)/2. )+1)*ielsign
    doel=real(omegabg(i))-real(omegag)/el
    doel2=real(omegabg(i))-real(omegag)/(el+2*ielsign)
    if(abs(doel2).lt.abs(doel))then
       el=el+2*ielsign
       doel=doel2
    endif
! Calculate the required omegabg imaginary part, which is that el*obi 
! must be at least |el|dob/Rc below imag(omegag) if omegabg is closer to
! the real resonance than Sc times the omegabg step size.
    dob=real(omegabg(i)-omegabg(i-1))
    obi=-ielsign*max(0.,dob/Rc-imag(omegag)/abs(el)) &
         *max(0.,1.-(doel/(Sc*dob))**2)
  end subroutine pathshiftg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FgEint(Ftotalg,isigma)
    complex :: Ftotalg
    if(psig.gt.0)then
       call FgRepelEint(Ftotalg,isigma)
    else
       call FgAttractEint(Ftotalg,isigma)
    endif
  end subroutine FgEint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Fextern2(Wgi,isigma,dvinf,dfweight)
! Revised Fextern that calculates analytically the external forces
    real :: Wgi
    complex :: Fextanal,Fext2anal,Fextqanal,Amp,dvdfw,dfweight
    vi=sqrt(2.*Wgi)
    dvdfw=dvinf*dfweight*sqm1g
    p=4.*kpar
    p2=p**2
    fn=sqrt(8.*3.1415926*(1+p2)*(4+p2)*(9+p2)*(16+p2)*(25+p2))
    Pk= -15*(p2**2 - 15*p2 + 8)
    Qk= p2**2 - 85*p2 + 274
    Amp=(Pk+complex(0.,p)*Qk)/fn
    Fextanal=conjg(Amp)*CPmds(ngz,1)*dvdfw*exp(-sqm1g*kpar*zg(ngz))&
         /(sqm1g*(kpar-omegag/vi))
    Fext2anal=conjg(Amp)*CPmds(ngz,2)*dvdfw*exp(-sqm1g*kpar*zg(ngz))&
         /(sqm1g*(kpar-omegag/vi))
    ! Analytic Difference <q|V-Vw|q>
    Fextqanal=conjg(Amp)*CPmds(ngz,3)*dvdfw*exp(-sqm1g*kpar*zg(ngz))&
         /(sqm1g*(kpar-omegag/vi))&
         +dvdfw*abs(Amp)**2/vi/(kpar-omegag/vi)**2
    Fextqqwanal=Fextqqwanal+Fextqanal ! Accumulated difference
    Ftmdp(3,1)=Ftmdp(3,1)+Fextanal
    Ftmdp(3,2)=Ftmdp(3,2)+Fext2anal
    Ftmdp(3,3)=Ftmdp(3,3)+Fextqanal
  end subroutine Fextern2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SumHarmonicsg(isigma)
    implicit none
    integer :: isigma
! Sum harmonic contributions to obtain Forces, equiv to v-perp integration. 
! Perpendicular k-vector kg matters only here, and limits minimum effective
! Omegacg to Oceff. 
    real :: EIm(0:nhmax),xit,vymax,hnum
    integer :: m,ncalc
    real :: Oceff   ! The effective Omegacg
! The maximum needed perp velocity, such that kg*vymax/Oc=nharmonicsg
    vymax=3.5*sqrt(2.*Tperpg)
    Oceff=max(1.e-6,max(Omegacg,kg*vymax/nhmax)) ! Don't allow zero Oceff.
    xit=kg*sqrt(Tperpg)/Oceff
! How many harmonics do we actually need? For large hnum:
    hnum=kg*vymax/Oceff                         ! This will not exceed nhmax
! If hnum is small, then use rather more for accuracy.
!    nharmonicsg=min(nhmax,int(hnum*(1.+3./(hnum+1.)))) ! Inadequate.
! Require the linear approx to Immin to be less than small.
    nharmonicsg=nint(max(hnum,min(4.,alog(.01)/alog(xit/2.+1.e-12))))
!    if(nharmonicsg.gt.nhmax)then
!       write(*,'(a,2f8.4,i3,4f7.3)')'kg,vymax,nharm,Oceff,Oc&
!            &,xit,psig',kg ,vymax,nharmonicsg,Oceff,Omegacg,xit,psig
!       stop 'Incorrect nharmonics exceeds nhmax'
!    endif
! Calculate the Integer[0.] exp*I[2] Bessel functions 0 to nharmonicsg
    ncalc=0
    call RIBESL(xit**2,0.,nharmonicsg+1,2,EIm,ncalc)
    if(.not.ncalc.eq.nharmonicsg+1)then ! All orders not calculated correctly.
       write(*,'(a,i3,a,i3,a)')'Bessel functions',ncalc,' of',nharmonicsg+1, &
       ' are precise. Others may be irrelevant.' 
    endif
    if(.not.(abs(EIm(0)).gt.-1.e-20))then   ! Insufficient test.
       write(*,*)'EIm NAN/Zero?',nharmonicsg+1,Oceff,xit**2,EIm(0)
       write(*,*)ncalc,EIm(0:nharmonicsg)
       stop
    endif
! m=0 always used.! fy is Maxwellian.
! But the fywy,fy need to be fixed in inner routines.
    omegaonly=omegag
    Ftotalsumg=0.
    Ftmdsum=0.
    FVwsumg=0.
    do m=1,nharmonicsg
       omegag=omegaonly+m*Oceff
       if(abs(omegag).lt.1.e-5)omegag=omegag+sqm1g*1.e-5 ! Prevent zero.
       call FgEint(Fpg(m),isigma)
       if(real(omegaonly).eq.0)then   !Short cut.
          Ftotalsumg=Ftotalsumg+2*real(Fpg(m))*EIm(m)
          Ftmdsum=Ftmdsum+2*Ftmda*EIm(m)
          FVwsumg=FVwsumg+(Fintqw+Fextqw)*2*EIm(m)
       else      ! Full sum over plus and minus m.
          Ftotalsumg=Ftotalsumg+Fpg(m)*EIm(m)
          Ftmdsum=Ftmdsum+Ftmda*EIm(m)
          FVwsumg=FVwsumg+(Fintqw+Fextqw)*2*EIm(m)
          omegag=omegaonly-m*Oceff
          if(abs(omegag).lt.1.e-5)omegag=omegag+sqm1g*1.e-5
          call FgEint(Fpg(-m),isigma)
          Ftotalsumg=Ftotalsumg+Fpg(-m)*EIm(m)
          Ftmdsum=Ftmdsum+Ftmda*EIm(m)
          FVwsumg=FVwsumg+(Fintqw+Fextqw)*2*EIm(m)
       endif
       if(lbess)write(*,'(a,i2,a,e11.4,''('',2e12.4,'')('',2e12.4,'')'')')&
            ' EI(',m,'),Fpg(+-m)=',EIm(m),Fpg(m),Fpg(-m)
    enddo
    omegag=omegaonly
    call FgEint(Fpg(0),isigma)
    if(lbess.and.nharmonicsg.gt.0)write(*,'(a,e11.4,''('',2e12.4&
         &,'')'',i4,f8.4)') ' EI(0),Ftt(0)  =',EIm(0),Fpg(0)&
         &,nharmonicsg,Oceff
    Ftotalsumg=Ftotalsumg+Fpg(0)*EIm(0)
    Ftmdsum=Ftmdsum+Ftmda*EIm(0)
    FVwsumg=FVwsumg+(Fintqw+Fextqw)*2*EIm(0)
  end subroutine SumHarmonicsg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function dfdWpar(vinf,fvinf)
! The derivative of the distant distribution func wrt energy W at 
! velocity vinf, in the rest frame of the hole. fvinf is returned the f(v). 
! Example. Maxwellian of temperature T. f=exp(-vinv^2/2T)/sqrt(2 pi T)
!    dfdWpar=-exp(-vinf**2/(2.*T))/T/sqrt(2.*3.1415926*T)
! Two half-density Maxwellians symmetrically shifted by vshift
! f=0.5*(exp(-(vinv-vshift)^2/2T)+exp(-(vinv+vshift)^2/2T))/sqrt(2 pi T))
! Even though for ions the temperature Torepel might be different from
! Tinf, we ignore that fact here, and correct for it in FgRepelEint
    Tp=Tinf
    e1=exp(-(vinf-vshift)**2/(2.*Tp))
    e2=exp(-(vinf+vshift)**2/(2.*Tp))
    fvinf=0.5*(e1+e2)/sqrt(2.*3.1415926*Tp)
    dfdWpar=0.5*(-e1*(vinf-vshift)-e2*(vinf+vshift)) &
         &/sign(max(abs(vinf*Tp),1.e-6),vinf)/sqrt(2.*3.1415926*Tp)
  end function dfdWpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hole potential form and trapped distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If some different form of phi is required. Replace these functions,
! e.g. with forms that call functions outside the module. 
  real function phigofz(zval)
    doubleprecision :: zby4
    zby4=zval/4.
    phigofz=real(psig/cosh(zby4)**4)
  end function phigofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function phigprimeofz(zval)
    phigprimeofz=-psig*sinh(zval/4.)/cosh(zval/4.)**5
  end function phigprimeofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex function mdofz(zval,imd,kpar)
! Eigenmodes of the potential, assumed to be sech^4(z/4) (z-normalized)
    real :: kpar,p2
    x=zval/4.
    S=1./cosh(x)
    T=tanh(x)
    if(imd.eq.1)then     ! 4
       mdofz=-3.*T*S**4*sqrt(70.)/16.
    elseif(imd.eq.2)then ! 2
       mdofz=T*S**2*(3*S**2-2)*sqrt(105.)/8
    elseif(imd.eq.3)then ! q
       p=4.*kpar
       p2=p**2
       fnorm=sqrt(8.*3.1415926*(1+p2)*(4+p2)*(9+p2)*(16+p2)*(25+p2))
       Pk= -15*(p2**2 + (28*S**2 - 15)*p2 + 63*S**4 - 56*S**2 + 8)
       Qk= p2**2 + (105*S**2 - 85)*p2 + 945*S**4 - 1155*S**2 + 274
! complex version, antisymmetric
       mdofz=(T*Pk+sign(1.,x)*complex(0.,p)*Qk)&
            *exp(complex(0.,p*abs(x)))/fnorm
       if(x.eq.0)mdofz=0.
    else
       mdofz=0.
       write(*,*)'ERROR: Incorrect mode number',imd
    endif    
  end function mdofz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real function dfdWptrap(Wj,fe)
! Return the energy derivative and the parallel distribution function
! Possibly including contribution from repelled species.
    real, parameter :: ri=1.4
    logical :: lfirst=.true.
    sqWj=sqrt(-Wj)
    fe=((2./pig)*sqWj+(15./16.)*Wj/sqrt(-psig)+experfcc(sqWj)&
       /sqrt(pig))/sqrt(2.)        ! New exact sech^4 hole form.
    dfdWptrap=((15./16.)/sqrt(-psig)-experfcc(sqWj)/sqrt(pig))/sqrt(2.)
! This is the approximate correction for ion charge applied only when
! we are evaluating an attracted species.
    if(lioncorrect.and.vrshift.ne.9999)then
       if(lfirst)then
          write(*,*)'WARNING ############ lioncorrect is TRUE ions are active.'
          lfirst=.false.
       endif
       vsx=1.3+0.2*psig
       vsa=vrshift     ! The vshift of the reflected species.
       denem1=(-1.+(vsa/vsx)**ri) /(1.+0.25*psig+vsa**2*(vsa/(vsx +(3.3&
            &/vsa)**1.5))**ri)*Tinf/Torepel
!       fdfac=(1-denem1*psig)   ! error.
       fdfac=(1-denem1)
       dfdWptrap=dfdWptrap*fdfac
! ftilde=fe-fflat is multiplied by fdfac. fflat=1/sqrt(2pi)
       fflat=1/sqrt(2.*pig)
       fe=fflat+fdfac*(fe-fflat)
!       write(*,'(a,4f10.4)')'lioncorrect',vrshift,denem1,fdfac
    endif
end function dfdWptrap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fvinfplot
! vshift, Tinf, defined in module shiftgen,     
    integer, parameter :: ngpl=100
    real, dimension(ngpl) :: vplinf,fplinf
    character*10 string
    vpmax=4.*sqrt(Tinf)+vshift
    do i=1,ngpl
       vplinf(i)=-vpmax+2.*i/ngpl*vpmax
       blah=dfdWpar(vplinf(i),fvinf)
       fplinf(i)=fvinf   ! This is really f, not v.
    enddo
    call autoplot(vplinf,fplinf,ngpl)
    call axlabels('v','f!di!A;!@!d')
    call fwrite(vshift,iwidth,2,string)
    call legendline(0.1,0.9,258,'v!ds!d='//string(1:iwidth))
    call pltend
  end subroutine fvinfplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for density calculation for diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dentadd(dfweight,dvinf)
    complex :: dfweight
! Add to dentpass the weighted average of CPmds over the zg
! points that reside in each of its (equally spaced) intervals.
! This is the ~V|4> density. When appropriately rescaled.
! Scaling is that using for |4> instead phigprime gives an unnormalized
! value |4>=8\psi/3sqrt(70)*|^4> (the normalized mode)
! The densities here accumulated are from normalized modes (including 4).

! Need to correct *vinf/v.
    vinfbyv=sqrt((Wg)/(-phi0d+Wg))
!    vinfbyv=1.
    call remesh(zg,CPmds(:,1)*f4norm,2*ngz+1,zdent,CapPhid,2*nzd+1)
    dentpass=dentpass+dvinf*sqm1g*dfweight*CapPhid *vinfbyv
! Similarly for the ~V|q>  passing density.
    call remesh(zg,CPmds(:,3),2*ngz+1,zdent,CapQd,2*nzd+1)
    dentq=dentq+dvinf*sqm1g*dfweight*CapQd *vinfbyv
! And the wave-generated internal passing density V_w|q>.
    denqwint=denqwint+dvinf*sqm1g*dfweight*(auxzd/&
         (sqm1g*(sign(1.,zdent)*kpar*vg(ngz)-omegag))) *vinfbyv
  end subroutine dentadd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dentaddtrap(dfweight,cdvinf)
! Make denttrap equal to the weighted average of \Phi over the zg
! points that reside in each of its (equally spaced) intervals.
! dfweight is omegag*dfe etc. ftilde is i*dfweight*\Phi (summed over
! the perpendicular harmonics). Correct for trapped resonance.
! The densities here accumulated are from normalized modes (including 4).
    complex :: dfweight,cdvinf
    complex :: resfac
    complex, dimension(-ngz:ngz) :: ft4
! Need to correct *vpsi/v. Limit how close to zero v can be.
    vpsibyv=sqrt((-psig+Wg)/max(-phi0d+Wg,2.e-5))
! Maybe we need to add in the form including prior weighting? dfwtprev=dfweight
    resfac=1+exp(sqm1g*omegag*taug(ngz)) ! Unshifted Half period version.
    ft4=sqm1g*dfweight*(CPmds(:,1)&
         -exp(sqm1g*omegag*taug)*CPmds(ngz,1)/resfac)
    call remesh(zg,ft4,2*ngz+1,zdent,ft4d,2*nzd+1)
    do j=1,nmd-1 ! Similarly for the modes trapped f
       ftraux(:,j)=sqm1g*dfweight*(CPmds(:,j+1)&
            -exp(sqm1g*omegag*taug)*CPmds(ngz,j+1)/resfac)
       call remesh(zg,ftraux(:,j),2*ngz+1,zdent,ftrauxd(:,j),2*nzd+1)
    enddo    
    denttrap=denttrap+cdvinf*ft4d *vpsibyv
    dentqt=dentqt+cdvinf*ftrauxd(:,2) *vpsibyv
  end subroutine dentaddtrap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine remesh(x1,y1,n1,x2,y2,n2)
! Interpolate arrays x1,y1, of length n1, where x1 is monotonic
! increasing, onto another monotonic increasing array x2 of length n2.
    integer n1,n2
    real :: x1(n1)
    complex :: y1(n1)
    real :: x2(n2)
    complex :: y2(n2)
! Assume both monotonic increasing for simplicity.    
    nx1=1
    do j=1,n2
       if(x2(j).lt.x1(nx1).or.x2(j).gt.x1(n1))then
          y2(j)=0. ! Make values outside the x1 range zero.
!          write(*,'(2i4,7f9.4)')nx1,j,x1(nx1),x2(j),ff,y1(nx1),y2(j)
       else
          do i=nx1,n1
             if(x1(i).ge.x2(j))then
                ff=(x2(j)-x1(i-1))/(x1(i)-x1(i-1))
                y2(j)=y1(i-1)+ff*(y1(i)-y1(i-1))
                nx1=max(i-1,1)
!                write(*,'(2i4,7f9.4)')i,j,x1(i),x2(j),ff,y1(i),y2(j)
                exit
             endif
          enddo          
       endif
    enddo
  end subroutine remesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qdenqint
! Integrate the continuum mode force using auxzd and dentq in the
! hole region -zm to zm.
    dz=zm/nzd
    Fintqq=-conjg(auxzd(-nzd))*dentq(-nzd)*dz/2.
    do i=-nzd,nzd
       Fintqq=Fintqq+conjg(auxzd(i))*dentq(i)*dz
    enddo
    Fintqq=Fintqq-conjg(auxzd(nzd))*dentq(nzd)*dz/2.
  end subroutine qdenqint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qwint
! Integrate the wave mode force in the hole region.
    dz=zm/nzd
    Fintqw=-conjg(auxzd(-nzd))*Vw*auxzd(-nzd)/2.*dz/2.
    do i=-nzd,nzd
       Fintqw=Fintqw+conjg(auxzd(i))*Vw*auxzd(i)/2.*dz
    enddo
    Fintqw=Fintqw-conjg(auxzd(nzd))*Vw*auxzd(nzd)/2.*dz/2.
  end subroutine qwint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shiftgen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines that use shiftgen and are thus the interface.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ionforce(Fi,omega,kin,Omegacin,psiin,vsin,isigma)
! Input parameters are in electron units, and produce no permanent
! changes to shiftgen internal parameters.
  use shiftgen
  complex :: Fi,omega,omg
  real :: psiin,vsin,Omegacin,kin
  if(vsin.eq.9999)then
     Fi=0.
     return
  endif
  kg=kin
!  write(*,*)'ionforce kg=',kg
  omg=omegag;omc=Omegacg;pg=psig;vs=vshift
  omegag=omega*sqrt(rmime)
  Omegacg=Omegacin/sqrt(rmime)
  psig=psiin
  vshift=vsin
  if(abs(isigma).ne.1) then
     write(*,*)'ERROR in ionforce. isigma=',isigma
     stop
  endif
  call SumHarmonicsg(isigma)
  Fi=Ftotalsumg
! Undo changes
  omegag=omg;Omegacg=omc;psig=pg;vshift=vs
end subroutine ionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine electronforce(Felec,omegain,kin,Omegacp,psiin,vsin,isigma)
! Inputs in electron units except for vsin in ion units. No permanent
! changes to shiftgen parameters. vshift is set temporarily to zero,
! and vrshift to ion shift vsin. 
! Turn off ion density correction by passing vsin=9999.
  use shiftgen
  real :: kin
  complex :: omegain,Felec
  omegag=omegain
  kg=kin
!  write(*,*)'electronforce kg=',kg
  Omegacg=Omegacp
  vr=vrshift
  vrshift=vsin   ! lioncorrect value
  vs=vshift
  vshift=0.
  psig=-psiin; call SumHarmonicsg(isigma); psig=-psig
  Felec=Ftotalsumg
  vshift=vs      ! Restore vshift and vrshift.
  vrshift=vr
end subroutine electronforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LBessSet
  use shiftgen
  lbess=.true.
end subroutine LBessSet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Tset(Te,Ti)
  use shiftgen
  Tinf=Te
  Torepel=Ti
end subroutine Tset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotfv(vsin)
  use shiftgen
  vshift=vsin
  call fvinfplot
end subroutine plotfv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function inharm()
  use shiftgen
  inharm=nharmonicsg
end function inharm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine makezdent
! Initialize the uniform zdent, triggering dent saving during FtEint.
! Make sure that it is never seen as being outside zm.
  use shiftgen
  do i=-nzd,nzd
     zdent(i)=.999999*zm*i/float(nzd)
     ! And the continuum auxzd internally on the zd mesh.  
     auxzd(i)=mdofz(zdent(i),3,kpar)
     phi0d(i)=phigofz(zdent(i)) ! and phi0d
     phipd(i)=phigprimeofz(zdent(i)) ! and phi0prime.
     zdmid(i)=zdent(i)
  enddo
  CapPhidprev=0.
end subroutine makezdent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Routines not dependent on shiftgen.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!****************************************************************
! This is exp(X^2)*erfc(X)
  FUNCTION expERFCC(X)
    Z=ABS(X)      
    T=1./(1.+0.5*Z)
    expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+                  &
         &    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+       &
         &    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
    IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
  END FUNCTION expERFCC
!********************************************************************
