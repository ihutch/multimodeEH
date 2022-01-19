module iterfind
! Provides iterfindroot for solving shiftmode dispersion relation by Newton iteration. 
! External calls are to electronforce, ionforce.
  integer, parameter :: niterfind=12
  real, dimension(0:niterfind) :: Frit,Fiit  !The sequence of complex roots for access.
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine iterfindroot(psip,vsin,Omegacp,omegap,kp,isigma,nit)
  real, intent(in) :: psip,vsin,Omegacp,kp
  integer, intent(in) :: isigma
  complex, intent(inout) :: omegap    ! The found root on out.
  integer, intent(out) :: nit         ! zero is failure, else n-iterations.
  integer :: nunconv=3
  complex :: Fec,Fic,Fsum
  logical :: litw=.false.
  zoif=.001  ! Iteration minimum oi limit factor.
  nzo=0
  omegap=complex(0.9*sqrt(psip)/8.,.9*sqrt(psip)/8./(1.+vsin))
  Frit(0)=max(real(omegap),-2e-3)
  Fiit(0)=imag(omegap)
  call electronforce(Fec,omegap,kp,Omegacp,psip,vsin,isigma)
  ienharm=inharm()
  call      ionforce(Fic,omegap,kp,Omegacp,psip,vsin,isigma)
  iinharm=inharm()
  FE=kp**2*psip**2*128./315.
  Fsum=Fec+Fic-FE
  err=0
  do i=1,niterfind
     if(litw)write(*,'(a,3i3,5f9.6)')'i,nharm,omegap,Fsum,err=',i&
          &-1,ienharm,iinharm,omegap,Fsum,err
     nit=i
     call complexnewton(Fsum,omegap,kp,err,psip,isigma,vsin,Omegacp)
     Frit(i)=max(real(omegap),-2e-3)
     Fiit(i)=imag(omegap)
     if(.not.abs(omegap).lt.1.e6)write(*,*)'Iterfindroot',i,psip,vsin,omegap
     if(imag(omegap).lt.zoif*sqrt(psip))then
        nzo=nzo+1
        if(nzo.ge.nunconv)then
           write(*,'(a,i2,a,g10.3)')'Uncoverged after',nzo,'&
                & omegai less than',zoif*sqrt(psip)
           err=1.
           nit=0
           goto 1
        endif
        omegap=complex(real(omegap),zoif*sqrt(psip))
     endif
     if(.not.abs(omegap).lt.1.e3)then
        write(*,*)'Iterfind diverging',omegap,err
        nit=0
        goto 1
     endif
     call electronforce(Fec,omegap,kp,Omegacp,psip,vsin,isigma)
     call      ionforce(Fic,omegap,kp,Omegacp,psip,vsin,isigma)
     Fsum=Fec+Fic-FE
     if(err.lt..5e-4)goto 1
     if(err*abs(omegap).lt.1.e-6)then
        write(*,*)'Apparent convergence at low omegap'
        goto 1
     endif
  enddo
  nit=0
  i=i-1
  write(*,*)'Unconverged after',i,' iterations'
1 continue
  if(litw)write(*,'(a,3i3,5f9.6)')'i,nharm,omegap,Fsum,err=',i&
       &,ienharm,iinharm,omegap,Fsum,err
end subroutine iterfindroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine complexnewton(Fsum,omegap,kp,err,psip,isigma,vsin,Omegacp)
  ! Take a Newton Step in finding the complex omega root of
  ! complex force by inverting the 2x2 matrix from adjacent evaluations.
  ! Return the Fsum,omegap,err new values.
  complex :: Fsum,omegap
  real :: err,psip,vsin,Omegacp,kp
  real :: eps1=.05,J11,J12,J21,J22,domega1,domega2,det
  complex :: om1,om2,om3,f1,f2,f3,domega,Fec,Fic
  det=1.
  FE=kp**2*psip**2*128./315.
  f1=Fsum
  om1=omegap
  ! Calculate the Jacobian's coefficients.
  om2=om1+eps1*abs(real(om1))
!write(*,*)'compnewt call, kp=',kp
  call electronforce(Fec,om2,kp,Omegacp,psip,vsin,isigma)
  call      ionforce(Fic,om2,kp,Omegacp,psip,vsin,isigma)
  f2=Fec+Fic-FE
!  write(*,'(6g12.4)')om1,om2,f2
  if(.not.om2-om1.ne.0)stop 'om2-om1=0'
  J11=real(f2-f1)/real(om2-om1)
  J21=imag(f2-f1)/real(om2-om1)
  om3=om1+complex(0.,eps1*abs(imag(om1)))
  call electronforce(Fec,om3,kp,Omegacp,psip,vsin,isigma)
  call      ionforce(Fic,om3,kp,Omegacp,psip,vsin,isigma)
  f3=Fec+Fic-FE
!  write(*,'(6g12.4)')om1,om3,f3
  if(.not.om3-om1.ne.0)stop 'om3-om1=0'
  J12=real(f3-f1)/imag(om3-om1)
  J22=imag(f3-f1)/imag(om3-om1)
  ! Solve J domega = -f to find domega
  det=J11*J22-J21*J12
  if(det.eq.0)then
     write(*,*)'Det=0'
     write(*,*)om1,om2,om3,f1,f2,f3
     write(*,'(2f10.4)')J11,J12,J21,J22
     ! Try to recover.
     det=1.
  endif
  domega1=-(J22*real(f1)-J12*imag(f1))/det
  domega2=-(-J21*real(f1)+J11*imag(f1))/det
!  write(*,*)J11*domega1+J12*domega2,J21*domega1+J22*domega2
  domega=complex(domega1,domega2)
  omegap=om1+domega
  err=abs(domega/omegap)
end subroutine complexnewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module iterfind
