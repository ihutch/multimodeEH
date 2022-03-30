! New version of omegacont.f90 to use for slow hole stability including
! the effects of ion force Fi.
! This version uses only shiftgen.
! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.
! k-variation is not yet actually implemented. 
! The variation of interest here is more about psip and vshift.

program fomegasolve
  use iterfind
  real :: kin,kmid
  logical :: lcont=.true.,lplot=.false.,lerase=.false.,lTiscan=.false.
  integer, parameter :: nk=1,noc=1,npsi=8,nvsin=21,nv0=0
  real :: Omegacarr(noc),karr(nk),psiparray(npsi),vsinarray(nv0:nvsin)
  real :: Tiarray(nvsin)
  real :: ormax,oimax
  complex :: omegasolve(nv0:nvsin,npsi),omegap
  character*30 string
  
! Default parameters
  kmid=0.
  Omegacmax=10
  rangek=.5    !fractional k-range
  Typ=1.
  

  isigma=-1
  psip=.25
  vsin=1.25         ! Maxwellian ion component velocity shift
  ormax=0.
  oimax=0.
  call parsefoarguments(psip,vsin,ormax,oimax,kmid,Omegacmax,lerase,lcont,lplot,Ti,lTiscan)
  if(ormax.eq.0)then
     ormax=.3*sqrt(psip)
  endif
  if(oimax.eq.0)then
     oimax=.1*sqrt(psip)
     if(vsin.le.1.2)oimax=.3*sqrt(psip)
     if(vsin.le..9)oimax=.4*sqrt(psip)
  endif
  omegap=complex(ormax,oimax)
  if(lplot)then
! contouring of psimax and vsmax case.       
     do ik=1,nk
        if(nk.gt.1)then
           kin=kmid*(1.+rangek*(2.*(ik-1.)/max(1.,nk-1.)-1.))
        else
           kin=kmid
        endif
        karr(ik)=kin/sqrt(psip)
!     write(*,*)'kin=',kin,' karr=',karr(ik)
        do ioc=1,noc
           Omegacp=ioc*Omegacmax/noc
           Omegacarr(ioc)=Omegacp/sqrt(psip)
           call fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot,err&
                &,ormax,oimax,lerase)
        enddo
     enddo
  elseif(lTiscan)then
! omega versus Ti at several vsin
     kin=0.
     Timax=Ti
     Timin=.5
     psip=.05
     omegasolve=0.
     nvs=4
     write(*,*)' psip   vshift   Ti    it     omega'
     do ip=1,nvs ! Really this is vs iterations here
        vs=vsin*(ip-1)/(nvs-1.)
        vsinarray(ip)=vs
        do iv=1,nvsin  ! Really this is Ti iterations
           Ti=Timin+(Timax-Timin)*(iv-1.)/(nvsin-.99999)
           if(ip.eq.1)Tiarray(iv)=Ti
           omegap=complex(0.7*sqrt(psip)/8.,.7*sqrt(psip)/8./(1.+vsin))
           call Tset(1.,Ti)
           call iterfindroot(psip,vs,Omegacp,omegap,kin,isigma,ires)
           if(ires.ne.0)omegasolve(iv,ip)=omegap/sqrt(psip)
           write(*,'(3f8.4,i3,$)')psip,vs,Ti,ires
           write(*,*)omegap
        enddo
     enddo
     call pltinit(Tiarray(1),Tiarray(nvsin),0.,.3+psip)
     call charsize(.019,.019)
     call axis; call axis2
     call axlabels('T!di!d','!Aw!@/!Ay!@!u1/2!u')
     call winset(.true.)
     call legendline(.55,.95,258,'v!ds!d=')
     iline=1
     do ip=1,nvs
        call color(ip)
        call fwrite(vsinarray(ip),iwidth,2,string)
        call dashset(2*ip-2)
        call polyline(Tiarray(1),real(omegasolve(1:,ip)),nvsin)
        if(real(omegasolve(1,ip)).gt.0.001)then
           call legendline(.5,.95-.05*iline,0,' '//string(1:iwidth)&
                &//' !Aw!@!dr!d')
           iline=iline+1
        endif
        call dashset(2*ip-1)
        call polyline(Tiarray(1),imag(omegasolve(1:,ip)),nvsin)
        call legendline(.5,.95-.05*iline,0,' '//string(1:iwidth)//' !Aw!@!di!d')
        iline=iline+1
        call dashset(0)
     enddo
     call pltend
     
  else
! omega versus vsin plots at various psi.
     kin=0.
     psipmax=psip
     vsmax=vsin
     Omegacp=Omegacmax
     lplot=.false.
     if(nvsin.lt.15)then
        vsmin=0.5
     else
        vsmin=0.1
     endif
     omegasolve=0.
     do ip=1,npsi
        psip=ip*psipmax/npsi
        psiparray(ip)=psip
        do iv=nv0,nvsin
           vsin=vsmin+(vsmax-vsmin)*(iv-1.)/(nvsin-.99999)
           if(iv.eq.0)vsin=0.
           if(ip.eq.1)vsinarray(iv)=vsin
           omegap=complex(0.7*sqrt(psip)/8.,.7*sqrt(psip)/8./(1.+vsin))
           call iterfindroot(psip,vsin,Omegacp,omegap,kin,isigma,ires)
           if(ires.ne.0)omegasolve(iv,ip)=omegap/sqrt(psip)
           write(*,'(2f8.4,i3,$)')psip,vsin,ires
           write(*,*)omegap
        enddo
     enddo
!  call minmax(omegasolve,2*npsi*nvsin,ommin,ommax)
!  call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,ommax*1.1)
! Now scaled version  
! Uncorrected  call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,.34)
     call pltinit(vsinarray(nv0),vsinarray(nvsin),0.,.3+psipmax)
     call charsize(.019,.019)
     call axis; call axis2
     call axlabels('!Bv!ds!d!@','!Aw!@/!Ay!@!u1/2!u')
     if(Ti.ne.1.)then
        call fwrite(Ti,iwidth,1,string)
        call legendline(.05,.07,258,'T!di!d='//string(1:iwidth))
     endif
     call winset(.true.)
     call legendline(.7,.95,258,'!Ay!@=')
     do ip=1,npsi
        call color(15)
        call fwrite(psiparray(ip),iwidth,3,string)
        call dashset(ip-1)
        call legendline(.65,.95-.05*ip,0,' '//string(1:iwidth))
        call color(1)
        call polyline(vsinarray,real(omegasolve(:,ip)),(nvsin+1-nv0))
        if(ip.eq.1)call legendline(.65,.5,258,'real(!Aw!@)/!Ay!@!u1/2!u')
        call color(4)
        call polyline(vsinarray,imag(omegasolve(:,ip)),(nvsin+1-nv0))
        if(ip.eq.1)call legendline(.1,.9,258,'imag(!Aw!@)/!Ay!@!u1/2!u')
        call dashset(0)
     enddo
     call color(0)
     call pltend
  endif
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parsefoarguments(psip,vsin,ormax,oimax,kin,Omegacmax,lerase,lcont,lplot,Ti,lTiscan)
  character*20 argument
  real :: kin
  logical :: lerase,lcont,lplot,lTiscan
  ipfset=3 ! default
  Ti=1.
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:2).eq.'-k')read(argument(3:),*)kin
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsin
     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-oc')read(argument(4:),*)Omegacmax
     if(argument(1:3).eq.'-Ti')then
        read(argument(4:),*)Ti
        write(*,'(a,f8.2)')'Setting Ti to',Ti
        call Tset(1.,Ti)
     endif
     if(argument(1:3).eq.'-lc')lcont=.not.lcont
     if(argument(1:3).eq.'-lp')lplot=.not.lplot
     if(argument(1:3).eq.'-lT')lTiscan=.not.lTiscan
     if(argument(1:2).eq.'-e')lerase=.not.lerase
     if(argument(1:2).eq.'-c')ipfset=-3
     if(argument(1:2).eq.'-h')goto 1
  enddo
  call pfset(ipfset)
  return
1 write(*,*)'-p psi, -vs vshift, -or -oi real, imag omega,',&
       ' -c no-stopping, -e erase file'
  write(*,*)'-Ti ion tempr, -oc Omegac -lp toggle on contours, -lt toggle Tiscan' 
  stop
end subroutine parsefoarguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot&
     &,err,ormax,oimax,lerase)
  use iterfind
  use mpiloops
  logical :: lerase
  complex :: omegap      ! omegag surrogate maybe read from file
  integer, parameter ::   nor=41,noi=41
  real :: or(nor),oi(noi),kin
  complex ::  omegacomplex(nor,noi),forcecomplex(nor,noi),Fi
  integer :: ienharm,iinharm
  complex ::  Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lcont,lplot
  logical :: lplot3=.false.,lplot2=.true.
  logical :: lions=.true.,lreadit=.false.
  character*30 string,filename,argument
  real ormax,oimax

  nproc=1
  id=0
  call mpilprep(id,nproc)
!  write(*,*)'fomegacont k=',kin
  ! Create filename in accordance with passed parameters
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a)')   &
       'F',nor,noi, &
       'Oc',min(999,abs(nint(100*Omegacp))),   &
       'k',min(999,abs(nint(100*kin))),   &
       'v',min(999,abs(nint(100*vsin))),  &
       'p',min(999,abs(nint(100*psip))),'.arr'
  do i=1,iargc()   ! Check cmdline for filenames.
     call getarg(i,argument)
     if(.not.argument(1:1).eq.'-')then
        filename=argument
        write(*,'(a,$)')'Specified file: '
     endif
  enddo
  ! Try to open the file.
  open(12,file=filename,status='old',form='unformatted',err=101)
  if(.not.lerase)then
     read(12,err=101)norf,noif
     if(norf.ne.nor.or.noif.ne.noi)then
        write(*,*)'Reading from ',filename
        write(*,*)'File array dimensions',norf,nori,' not compatible'
        write(*,*)'Adjust allocation or delete file'
        stop
     endif
     read(12,err=101)or,oi
     read(12,err=101)psip,Omegacp,Typ,omegap,kin,ormax,oimax
     read(12,err=101)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
     read(12,end=100)lions
     close(12)
100  lreadit=.true.
     if(id.eq.0)write(*,*)'Read forcecomplex from file ',filename
     if(id.eq.0)write(*,*)'kin=',kin,' omegap=',omegap
     call mpilstopslaves
  else
     if(id.eq.0)write(*,*)'Overwriting forcecomplex file ',filename
     close(12,status='delete')
  endif
  goto 102
101 if(id.eq.0)write(*,*)'Failed to open or read from file: ', filename
102 continue


  if(lplot.and.id.eq.0)call plotfv(vsin)  
  isigma=-1
  FE=kin**2*psip**2*128./315.

  if(.not.lreadit.and.lcont)then    ! Failed to read from file so calculate
     impi=0
     dioi=-0.04   ! Offset of oi(1) from zero.
     call mpilprep(id,nproc)
     ! Contruct the forcecomplex matrix
     if(id.eq.0)write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr &
          & Ftotali     k   nharms'
     do ioi=1,noi
        oi(ioi)=(ioi-1+dioi)*oimax/(noi-1+dioi)
        do ior=1,nor
           or(ior)=(ior-1)*ormax/(nor-1.)
           omegap=complex(or(ior),oi(ioi))
           omegacomplex(ior,ioi)=omegap
           
           impi=impi+1
           iactiv=mod(impi,nproc)  ! Decide the active rank process iactiv
           if(iactiv.eq.id)then    ! If I am active I do the work needed ...
              
              call electronforce(Ftcomplex(ior,ioi),omegacomplex(ior&
                   &,ioi),kin,Omegacp,psip,vsin,isigma)
              ienharm=inharm()
              if(lions)then
                 call ionforce(Fi,omegacomplex(ior,ioi) ,kin,Omegacp,&
                      & psip,vsin ,isigma)
                 iinharm=inharm()
              endif
              Ficomplex(ior,ioi)=Fi
              forcecomplex(ior,ioi)=Ftcomplex(ior,ioi)+Fi-FE
!              write(*,'(2i4,2f8.4,2f10.6,f7.3,3i4)')ior,ioi,omegacomplex(ior&
!                   &,ioi),forcecomplex(ior,ioi),kin,ienharm,iinharm,iactiv
              
           endif
           call mpilcommscomplex(forcecomplex(ior,ioi),iactiv,1,impi)
           call mpilcommscomplex(Ftcomplex(ior,ioi),iactiv,1,impi)
           call mpilcommscomplex(Ficomplex(ior,ioi),iactiv,1,impi)
! We probably don't need to pass ienharm etc.. Values remain what they
! were last time master did the work.
           if(id.eq.0)write(*,'(2i4,2f8.4,2f10.6,f7.3,2i4)')ior&
                &,ioi,omegacomplex(ior ,ioi),forcecomplex(ior,ioi)&
                &,kin,ienharm,iinharm
        enddo
     enddo
     call mpilstopslaves        ! Prevent slaves from continuing.
     write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali     k   nharms'
     write(*,*)'Omegacp,k,psip',Omegacp,kin,psip
     
        ! Attempt to write but skip if file exists.
     open(12,file=filename,status='new',form='unformatted',err=103)
     write(*,*)'Opened new file: ',filename,' and writing'
     write(12)nor,noi
     write(12)or,oi
     write(12)psip,Omegacp,Typ,omegap,kin,ormax,oimax
     write(12)omegacomplex,forcecomplex,Ftcomplex,Ficomplex
     write(12)lions
     goto 104
103  write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104  close(12)
  endif

  write(*,*)'Omegacp/omegab=',2*Omegacp/sqrt(psip)
  if(lplot)then
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ftcomplex/psip**2)
     call legendline(.1,.9,258,'!p!o~!o!qF!de!d/!Ay!@!u2!u')
     call pltend
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ficomplex/psip**2)
     call legendline(.1,.9,258,'!p!o~!o!qF!di!d/!Ay!@!u2!u')
     call pltend

     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,forcecomplex/psip**2)
     if(FE.eq.0.)then
        call legendline(.1,.9,258,'!p!o~!o!qF/!Ay!@!u2!u')
     else
        call legendline(.1,.9,258,'(!p!o~!o!qF-F!dE!d)/!Ay!@!u2!u')
     endif
  endif
     
! Find root and plot it converging (omegag is set to found omegap implicitly)  
!  omegap=complex(0.7*sqrt(psip)/8.,1.*sqrt(psip)/8./(1.+vsin))
  write(*,*)'calling iterfindroot',omegap,kin
  call iterfindroot(psip,vsin,Omegacp,omegap,kin,isigma,ires)
  if(ires.ne.0)then
     do j=0,ires
        call color(6)
        call polymark(max(Frit(j),-2e-3),Fiit(j),1,ichar('0')+j)
     enddo
  endif
  write(*,*)'Eigenfrequency=',omegap
  if(lplot)     call pltend
  if(lplot3)then
     call multiframe(2,1,3)
     call ocomplot(or,nor,vsin,omegacp,psip,(Ftcomplex(:,1))/psip**2)
     call legendline(.1,.9,258,'F!de!d at imag(!Aw!@)=0')
!  call orealplot(or,nor,vsin,omegacp,psi,real(Ficomplex(:,1)))
     call ocomplot(or,nor,vsin,omegacp,psip,(Ficomplex(:,1)/psip**2))
     call legendline(.1,.9,258,'F!di!d')
     call multiframe(0,0,0)
     call pltend()
  endif
  if(lplot2)then
     sqpsi=sqrt(psip)
     tqpsi=psip**0.75
     uqpsi=psip**1.5
     qqpsi=psip**.25
     call pltinit(0.,ormax/tqpsi,0.,oimax/uqpsi)
     call charsize(0.02,0.02)
     call axis
     call axis2
     call axlabels('!Aw!B!dr!d/!Ay!@!u3/4!u','!Aw!B!di!d!@/!Ay!@!u3/2!u')
     icsw=1
     icl=0
     zclv(1)=10
     call color(1)
     call dashset(2)
     call contourl(real(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(0.1,-.1,0,'real')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(real(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     icl=0
     zclv(1)=20
     call color(2)
     call dashset(3)
     call contourl(imag(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call legendline(.65,-.1,0,'imag')
     icl=-1
     zclv(1)=0.
     call dashset(0)
     call contourl(imag(forcecomplex)/psip**2.5,cworka,nor,nor,noi,zclv,icl, &
          or/tqpsi,oi/uqpsi,icsw)
     call color(15)
     call fwrite(kin/qqpsi,iwidth,3,string)
     call legendline(0.05,1.04,258,'k/!Ay!@!u1/4!u='//string)
     call fwrite(omegacp,iwidth,2,string)
     call legendline(.45,1.04,258,'!AW!@='//string)
     call fwrite(psip,iwidth,2,string)
     call legendline(.8,1.04,258,'!Ay!@='//string)
     call color(6)
!     write(*,*)omegap,.001*sqrt(psip)
     if(imag(omegap).gt.0.0011*sqrt(psip))& 
          call polymark(real(omegap/tqpsi),imag(omegap/uqpsi),1,3)
     call color(15)
     call pltend
  endif
end subroutine fomegacont
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!include 'iterfindold.f90'  ! Obsolete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plotting routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lplot1(or,oi,nor,noi,vsin,omegacp,kin,psi,forcecomplex)
  real :: or(nor),oi(noi),kin
  complex ::  forcecomplex(nor,noi)
!  complex ::  Fpcomplex(nor,noi),Ftcomplex(nor,noi),Ficomplex(nor,noi)
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  character*30 string
  call pltinit(0.,or(nor),0.,oi(noi))
  call charsize(0.02,0.02)
  call axis
  call axis2
  call axlabels('!Aw!B!dr!d','!Aw!B!di!d!@')
  
  icsw=1
  call color(1)
  call dashset(4)
  icl=0
  zclv(1)=20
  call contourl(real(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
  call legendline(0.1,-.1,0,'real')
  icl=-1
  zclv(1)=0.
  call color(1)
  call dashset(0)
  call contourl(real(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
  icl=0
  zclv(1)=20
  call color(2)
  call dashset(2)
  call contourl(imag(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
  call legendline(.6,-.1,0,'imag')
  call dashset(0)
  icl=-1
  zclv(1)=0.0000
  call color(2)
  call dashset(0)
  call contourl(imag(forcecomplex),cworka,nor,nor,noi,zclv,icl,or,oi,icsw)
  call color(15)
!        call fwrite(k,iwidth,3,string)
!        call legendline(0.1,1.04,258,'k='//string)
  call fwrite(vsin,iwidth,2,string)
  call legendline(0.02,1.04,258,'v!ds!d='//string)
  call fwrite(omegacp,iwidth,2,string)
  call legendline(.27,1.04,258,'!AW!@='//string)
  call fwrite(kin,iwidth,2,string)
  call legendline(.52,1.04,258,'k='//string)
  call fwrite(psi,iwidth,3,string)
  call legendline(.77,1.04,258,'!Ay!@='//string)
  call dashset(0)
end subroutine lplot1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orealplot(or,nor,vsin,omegacp,psi,orealforce)
  real ::  or(nor)
  real ::  orealforce(nor)

  call minmax(orealforce,nor,fmin,fmax)
  call pltinit(0.,or(nor),fmin,fmax)
  call axis; call axis2
  call axlabels('real(!Aw!@)','real(Force)')
  call polyline(or,orealforce,nor)
end subroutine orealplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ocomplot(or,nor,vsin,omegacp,psi,ocomforce)
  real ::  or(nor)
  complex ::  ocomforce(nor)
  call minmax(ocomforce,2*nor,fmin,fmax)
  call pltinit(0.,or(nor),fmin,fmax)
  call axis; call axis2
  call axlabels('real(!Aw!@)','Force')
  call polyline(or,real(ocomforce),nor)
  call legendline(.1,.8,0,'real')
  call dashset(1)
  call polyline(or,imag(ocomforce),nor)
  call legendline(.1,.7,0,'imag')
  call dashset(0)
end subroutine ocomplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotionforce(psi,Typ,vsin,Omegacin)
! Repurposed for Ftot and Fe.
  real :: psi,Typ,vsin
  integer, parameter :: nfi=100
  complex, dimension(nfi) :: Ftotarray,Fiarray,Fearray,omegaFi
  omegamax=3.
  write(*,*)'psi=',psi,' vsin=',vsin
  do i=1,nfi
     omegaFi(i)=omegamax*(float(i)/nfi)/sqrt(1836.)+complex(0.,.0001)
     call ionforce(Fiarray(i),omegaFi(i),kin,Omegacin,psi,vsin,isigma)
     call electronforce(Fearray(i),omegaFi(i),kin,Omegacin,psi,vsin,-1)
     Ftotarray(i)=Fearray(i)+Fiarray(i)
  enddo
  call minmax(Ftotarray,2*nfi,fmin,fmax)
  call pltinit(0.,omegaFi(nfi),fmin,fmax)
  call axis; call axis2
  call axlabels('omega','Ftot, Fe')
  call color(1)
  call polyline(real(omegaFi),real(Ftotarray),nfi)
  call color(2)
  call dashset(1)
  call polyline(real(omegaFi),imag(Ftotarray),nfi)
  if(.true.)then
  call color(3)
  call dashset(2)
  call polyline(real(omegaFi),real(Fearray),nfi)
  call color(4)
  call dashset(3)
  call polyline(real(omegaFi),imag(Fearray),nfi)
  endif
  call dashset(0)
  call color(15)
  call pltend
end subroutine plotionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
