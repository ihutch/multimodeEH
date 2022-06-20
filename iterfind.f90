module iterfind
! Provides iterfindroot for solving shiftmode dispersion relation by
! Newton iteration.  Modified to iterate on the dispersion matrix
! determinant when lgetdet is true.
! External calls are to electronforce, ionforce.
  integer, parameter :: niterfind=12
  real, dimension(0:niterfind) :: Frit,Fiit  !The sequence of complex roots.
  logical :: lgetdet=.true.,litw=.true.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine iterfindroot(psip,vsin,Omegacp,omegap,kp,isigma,nit)
  real, intent(in) :: psip,vsin,Omegacp,kp
  integer, intent(in) :: isigma
  complex, intent(inout) :: omegap    ! The found root on out.
  integer, intent(out) :: nit         ! zero is failure, else n-iterations.
  integer :: nunconv=5
  complex :: Fec,Fic,Fm1,Fm2,om1
  complex, external :: getdet
  ! kp is k not k/sqrt(psip).
  zoif=.0001  ! Iteration minimum oi limit factor.
  nzo=0
  FE=kp**2*psip**2*128./315.
  Fic=0.
  Frit(0)=max(real(omegap),-2e-3)
  Fiit(0)=imag(omegap)
  om1=omegap*1.05
  call      ionforce(Fic,om1,kp,Omegacp,psip,vsin,isigma)
  iinharm=inharm()
  call electronforce(Fec,om1,kp,Omegacp,psip,vsin,isigma)
  ienharm=inharm()
  Fm1=Fec+Fic-FE
  if(lgetdet)Fm1=getdet()                !Iterate on |M| not Fsum
  call      ionforce(Fic,omegap,kp,Omegacp,psip,vsin,isigma)
  iinharm=inharm()
  call electronforce(Fec,omegap,kp,Omegacp,psip,vsin,isigma)
  ienharm=inharm()
  if(litw)write(*,'(a)')'  i  nharm      omega       oerr   |Fe+Fi=FE|  ||M|| '
  Fm2=Fec+Fic-FE
  if(lgetdet)Fm2=getdet()
  err=abs((omegap-om1)/om1)
  do i=1,niterfind
     if(litw)then
        write(*,'(3i3,8f9.6)')i-1 ,ienharm,iinharm,omegap,err,abs(Fec+Fic-FE),abs(getdet())
     endif
     nit=i
     call complexdnewton(Fm1,Fm2,om1,omegap)
! The following limits how big a step can be taken when far from convergence.
     omegap=om1+min(abs(omegap-om1),0.03*sqrt(abs(psip)))/abs(omegap-om1)*(omegap-om1)
     Frit(i)=max(real(omegap),-2e-3)
     Fiit(i)=imag(omegap)
     if(.not.abs(omegap).lt.1.e6)write(*,*)'Iterfindroot',i,psip,vsin,omegap
     if(imag(omegap).lt.zoif*sqrt(psip))then
        zoif=zoif*max(1,nzo)
        nzo=nzo+1
!        write(*,*)'nzo=',nzo,' zoif=',zoif,' oi=',imag(omegap)
        if(nzo.ge.nunconv.and.imag(omegap).le.0.)then
           write(*,'(a,i2,a,g10.3,a,g10.3)')'Unconverged after',nzo,'&
                & omegai (',imag(omegap),') less than',zoif*sqrt(psip)
           err=1.
           nit=0
           goto 1
        endif
        zoif=zoif/nzo
        omegap=complex(real(omegap),zoif*sqrt(psip))
     endif
     if(.not.abs(omegap).lt.1.e3)then
        write(*,*)'Iterfind diverging',omegap,err
        nit=0
        goto 1
     endif
     err=abs((omegap-om1)/om1)
     call      ionforce(Fic,omegap,kp,Omegacp,psip,vsin,isigma)
     call electronforce(Fec,omegap,kp,Omegacp,psip,vsin,isigma)
     Fm2=Fec+Fic-FE
     if(lgetdet)Fm2=getdet()
     if(err.lt..5e-4)goto 1
     if(err*abs(omegap).lt.1.e-6)then
        write(*,*)'Apparent convergence at low omega'
        goto 1
     endif
  enddo
  nit=0
  i=i-1
  write(*,*)'Unconverged after',i,' iterations'
1 continue
  if(litw)write(*,'(3i3,5f9.6)')i,ienharm,iinharm,omegap,err
end subroutine iterfindroot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine complexdnewton(F1,F2,x1,x2)
! Take a Newton iteration step based on the previous 2 complex values.
  complex :: F1,F2,x1,x2,dfdx
  dfdx=(F2-F1)/(x2-x1)
  x1=x2
  F1=F2
  x2=x2-F2/dfdx
end subroutine complexdnewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module iterfind
