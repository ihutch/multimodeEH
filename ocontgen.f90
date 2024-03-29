! New version of omegacont.f90 to use for slow hole stability including
! the effects of ion force Fi.
! This version uses only shiftgen.
! Contour the real and imaginary force(s) over a complex omega domain.
! For k, psi and Omegac having some values.
! Solve for the complex omega that makes the complex force zero.
! Stripped down versoin of slowstabgen.
! Change kin to refer to k/sqrt(psi). Change Omegacp to refer to O/sqrt(psi)

program fomegasolve
  use iterfind
  real :: kin,kmid
  logical :: lcont=.true.,lplot=.true.,lerase=.false.,lTiscan=.false.
  integer, parameter :: nk=1,noc=1,npsi=8,nvsin=21,nv0=0
  real :: Omegacarr(noc),karr(nk)!,vsinarray(nv0:nvsin),psiparray(npsi)
!  real :: Tiarray(nvsin)
  real :: ormax,oimax
  complex :: omegap!,omegasolve(nv0:nvsin,npsi)
!  character*30 string
  
! Default parameters
  kmid=0.
  Omegacmax=10
  rangek=.5    !fractional k-range
  Typ=1.
  call setnmd(3,0)  ! nmds, skip2?
  
  isigma=-1
  psip=.25
!  vsin=1.25         ! Maxwellian ion component velocity shift
  vsin=9999          ! No ions. ! Also lions=.false.
  ormax=0.
  oimax=0.
  call parsefoarguments(psip,vsin,ormax,oimax,kmid,Omegacmax,lerase,lcont,lplot,Ti,lTiscan,lgetdet)
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
        karr(ik)=kin
!     write(*,*)'kin=',kin,' karr=',karr(ik)
        do ioc=1,noc
           Omegacp=ioc*Omegacmax/noc
           Omegacarr(ioc)=Omegacp
           call fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot,err&
                &,ormax,oimax,lerase)
        enddo
     enddo
  endif
end program fomegasolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parsefoarguments(psip,vsin,ormax,oimax,kin,Omegacmax,lerase,lcont,lplot,Ti,lTiscan,lgetdet)
  character*20 argument
  real :: kin
  logical :: lerase,lcont,lplot,lTiscan,lgetdet
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
     if(argument(1:2).eq.'-M')lgetdet=.not.lgetdet
     if(argument(1:2).eq.'-h')goto 1
  enddo
  call pfset(ipfset)
  return
1 write(*,*)'-p psi, -vs vshift, -or -oi real, imag omega,',&
       ' -c no-stopping, -e erase file'
  write(*,*)'-Ti ion tempr, -oc Omegac, -k kp. -lp toggle on contours, -lt toggle Tiscan'
  write(*,*)'-M toggle iterate determinant -e overwrite file'
  stop
end subroutine parsefoarguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fomegacont(psip,Omegacp,Typ,omegap,kin,vsin,lcont,lplot&
     &,err,ormax,oimax,lerase)
  use iterfind
  use mpiloops
  logical :: lerase
  complex :: omegap      ! omegag surrogate maybe read from file
  integer, parameter ::   nor=21,noi=21
  real :: or(nor),oi(noi),kin
  complex ::  omegacomplex(nor,noi),forcecomplex(nor,noi),Fi
  integer :: ienharm,iinharm
  complex, dimension(nor,noi) ::  Ftcomplex,Ficomplex,DispDet,DetCap
  real, dimension(nor,noi) :: cworka
  integer :: icl
  real :: zclv(20)
  logical :: lcont,lplot
  logical :: lplot3=.false.,lplot2=.true.
  logical :: lions=.false.,lreadit=.false.
  character*30 string,filename,argument
  real ormax,oimax
  complex, external :: getdet
  
  id=0
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
  open(11,file=filename,status='old',form='unformatted',err=101)
  write(*,*)'Opened file:',filename
  if(.not.lerase)then
     read(11,err=105)norf,noif
     if(norf.ne.nor.or.noif.ne.noi)then
        write(*,*)'Reading from ',filename
        write(*,*)'File array dimensions',norf,nori,' not compatible'
        write(*,*)'Adjust allocation or delete file'
        stop
     endif
     read(11,err=105)or,oi
     read(11,err=105)psip,Omegacp,Typ,omegap,kin,ormax,oimax
     read(11,err=105)omegacomplex,forcecomplex,Ftcomplex,Ficomplex,DispDet
     read(11,end=100)lions
     close(11)
     goto 100
105  write(*,*)'Read error from ',filename,' Delete and recalculate.'
     close(11,status='delete')
     goto 101
100  lreadit=.true.
     call mpilserial  ! Remove for mpiexec operation else multiple plots.
     call mpilprep(id,nproc) ! Only needed if mpiexec used.
     call mpilkillslaves
     if(id.eq.0)write(*,*)'Read forcecomplex from file ',filename,' id=',id
     if(id.eq.0)write(*,'(a,f8.4,a,2f10.5,a,f8.4)')'kin=',kin,' omegmaxes=',omegap,' psip',psip
  else
     call mpilprep(id,nproc) ! Only needed if mpiexec used.
     if(id.eq.0)then
        write(*,*)'Overwriting forcecomplex file ',filename
        close(11,status='delete')
     endif
  endif
  goto 102
  
101 if(id.eq.0)write(*,*)'Failed to open or read from file: ', filename
102 continue


!  if(lplot.and.id.eq.0)call plotfv(vsin)  
  isigma=-1
  FE=kin**2*psip**3*128./315.  ! kin refers to k/sqrt(psi).

  if(.not.lreadit.and.lcont)then    ! Failed to read from file so calculate
     impi=0
     dioi=min(1.5*kin,.8)  ! Offset of oi(1) from zero.
     call mpilprep(id,nproc)
     ! Construct the forcecomplex matrix
     if(id.eq.0)write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr &
          & Ftotali     k   nharms'
     DetPhase=0.
     do ioi=1,noi
        oi(ioi)=(ioi-.98+dioi)*oimax/(noi-1+dioi)
        do ior=1,nor
           or(ior)=(ior-.98)*ormax/(nor-1.)  ! just avoid zero.
           omegap=complex(or(ior),oi(ioi))
           omegacomplex(ior,ioi)=omegap
           
           impi=impi+1
           iactiv=mod(impi,nproc)  ! Decide the active rank process iactiv
           if(iactiv.eq.id)then    ! If I am active I do the work needed ...
              
              if(lions)then
                 call ionforce(Fi,omegacomplex(ior,ioi),kin*sqrt(psip)&
                      ,Omegacp*sqrt(psip), psip,vsin ,isigma)
                 iinharm=inharm()
              endif
              Ficomplex(ior,ioi)=Fi
              call electronforce(Ftcomplex(ior,ioi),omegacomplex(ior&
                   &,ioi),kin*sqrt(psip),Omegacp*sqrt(psip),psip,vsin,isigma)
              ienharm=inharm()
              forcecomplex(ior,ioi)=Ftcomplex(ior,ioi)+Fi-FE
              DispDet(ior,ioi)=getdet()
              if(.false..and.ioi.eq.1)then
                 write(*,*)ior,or(ior),ioi,oi(ioi)
                 call printmatrix(phip)
              endif
           endif
           call mpilcommscomplex(forcecomplex(ior,ioi),iactiv,1,impi)
           call mpilcommscomplex(Ftcomplex(ior,ioi),iactiv,1,impi)
           call mpilcommscomplex(Ficomplex(ior,ioi),iactiv,1,impi)
           call mpilcommscomplex(DispDet(ior,ioi),iactiv,1,impi)
! We probably don't need to pass ienharm etc.. Values remain what they
! were last time master did the work.
           if(id.eq.0)write(*,'(2i4,2f8.4,2f10.6,f7.3,2i4)')ior&
                &,ioi,omegacomplex(ior ,ioi),forcecomplex(ior,ioi)&
                &,kin,ienharm,iinharm
        enddo
     enddo
     call mpilkillslaves        ! Prevent slaves from continuing.
     write(*,'(a)')'ior,  ioi  omegar omegai   Ftotalr  Ftotali     k   nharms'
     write(*,*)'Omegacp/sqrt(psi),k/sqrt(psip),psip',Omegacp,kin,psip
     
        ! Attempt to write but skip if file exists.
     open(12,file=filename,status='new',form='unformatted',err=103)
     write(*,*)'Opened new file: ',filename,' and writing'
     write(12)nor,noi
     write(12)or,oi
     write(12)psip,Omegacp,Typ,omegap,kin,ormax,oimax
     write(12)omegacomplex,forcecomplex,Ftcomplex,Ficomplex,DispDet
     write(12)lions
     goto 104
103  write(*,*)'New File: ',filename,' cannot be opened; not rewriting.'
104  close(12)
  endif

! Often confusing:  
!  write(*,'(10f8.4)')imag(DispDet(:,1))
!  write(*,'(10f8.4)')imag(DispDet(:,2))
  Fscale=(9*70.)/16**2/psip**2
  write(*,*)'Omegacp/omegab=',2*Omegacp
  if(lplot)then
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ftcomplex*Fscale)
     call legendline(.1,.9,258,'!p!o~!o!qF!de!d/!Ay!@!u2!u')
     call pltend
     if(lions)then
        call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,Ficomplex*Fscale)
        call legendline(.1,.9,258,'!p!o~!o!qF!di!d/!Ay!@!u2!u')
        call pltend
     endif
     call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,forcecomplex*Fscale)
     if(FE.eq.0.)then
        call legendline(.05,.9,258,'!p!o~!o!qF/!Ay!@!u2!u')
     else
        call legendline(.05,.95,258,'M!d11!d')
     endif
  endif
     
! Find root and plot it converging (omegag is set to found omegap implicitly)  
  omegap=complex(min(.5*kin*sqrt(psip),ormax*.7),min(kin*sqrt(psip)/2,oimax/2.))
  if(.false.)then
  lgetdet=.false.
  write(*,*)'calling M11 iterfindroot',omegap,kin,lgetdet
  call iterfindroot(psip,vsin,Omegacp*sqrt(psip),omegap,kin&
       &*sqrt(psip),isigma,ires)
  call printmatrix(psip)
  lgetdet=.true.
  endif
  write(*,*)'calling iterfindroot',omegap,kin,lgetdet
  call iterfindroot(psip,vsin,Omegacp*sqrt(psip),omegap,kin&
       &*sqrt(psip),isigma,ires)
  call dashset(1)
  call color(5)
  call contourl(real(DispDet(:,1:noi)),cworka,nor,nor,noi-1,0.,1,or,oi(1:noi),1)
  call jdrwstr(wx2nx(real(omegap)),wy2ny(imag(omegap))+.06,'real(D)=0',1.1)
  call color(7)
  call contourl(imag(DispDet(:,1:noi)),cworka,nor,nor,noi-1,0.,1,or,oi(1:noi),1)
  if(real(omegap).gt.1.e-5)&
      call jdrwstr(wx2nx(real(omegap)),wy2ny(imag(omegap))+.03,'imag(D)=0',-1.)
  call dashset(0)
  call pltend
  call printmatrix(psip)
!  write(*,'(10f8.4)')real(DispDet(1:10,1:10))

! truncate range of DispDet probably no longer needed.
  dmax=abs(DetCap(int(nor/2.),int(noi/2.))-DetCap(2,int(noi/2.)))
!  write(*,*)'dmax=',dmax,DetCap(int(nor/2.),int(noi/2.)),DetCap(2,int(noi/2.))
!  dmax=0.1
  DetCap=DispDet
  do i=1,nor
     do j=1,noi
        if(abs(real(DispDet(i,j))).gt.dmax)DetCap(i,j)&
             =complex(sign(dmax,real(DispDet(i,j))),imag(DetCap(i,j)))
        if(abs(imag(DispDet(i,j))).gt.dmax)DetCap(i,j)&
             =complex(real(DetCap(i,j)),sign(dmax,imag(DispDet(i,j))))
     enddo
  enddo
!  call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,DetCap) !Maybe unneeded.
  if(.false.)then
  call autoplot(or,imag(forcecomplex(:,1)*Fscale),noi)
  call axlabels('omegar','forcecomplex,DispDet')
  call dashset(1)
  call polyline(or,imag(DispDet(:,1)),noi)
  call dashset(0)
  call pltend

  call pltinit(0.,or(nor),-3.,3.)
  call axis
  call axlabels('omegar','DispDet/Forcecomplex')
  call polyline(or,imag(DispDet(:,1))/imag(forcecomplex(:,1)*Fscale),nor)
  call pltend
  endif
  
  call lplot1(or,oi,nor,noi,vsin,omegacp,kin,psip,DispDet)
  call legendline(.05,.95,258,'D')
  if(ires.ne.0)then
     do j=0,ires
        call color(6)
        call polymark(max(Frit(j),-2e-3),Fiit(j),1,ichar('0')+j)
     enddo
  endif
  call color(5)
  call dashset(1)
  call contourl(real(forcecomplex(:,1:noi)*Fscale),&
       cworka,nor,nor,noi-1,0.,1,or,oi(1:noi),1)
  call jdrwstr(wx2nx(real(omegap)),wy2ny(imag(omegap))+.06,&
       'real(M!d11!d)=0',1.1)
  call color(7)
  call contourl(imag(forcecomplex(:,1:noi)*Fscale),&
       cworka,nor,nor,noi-1,0.,1,or,oi(1:noi),1)
  if(real(omegap).gt.1.e-5)&
       call jdrwstr(wx2nx(real(omegap)),wy2ny(imag(omegap))+.03,&
       'imag(M!d11!d)=0',-1.1)
  call dashset(0)
  call pltend
  write(*,*)'Eigenfrequency=',omegap
!  if(lplot)     call pltend
  if(lplot3)then
     call multiframe(2,1,3)
     call ocomplot(or,nor,vsin,omegacp,psip,(Ftcomplex(:,1))*Fscale)
     call legendline(.1,.9,258,'F!de!d at imag(!Aw!@)=0')
!  call orealplot(or,nor,vsin,omegacp,psi,real(Ficomplex(:,1)))
     call ocomplot(or,nor,vsin,omegacp,psip,(Ficomplex(:,1)*Fscale))
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
     call fwrite(kin*qqpsi,iwidth,3,string)
     call legendline(0.05,1.04,258,'k/!Ay!@!u1/4!u='//string)
     call fwrite(omegacp,iwidth,2,string)
     call legendline(.45,1.04,258,'!AW!@='//string)
     call fwrite(psip,iwidth,2,string)
     call legendline(.8,1.04,258,'!Ay!@='//string)
     call color(6)
!     write(*,*)omegap,.001*sqrt(psip)
!     if(imag(omegap).gt.0.0011*sqrt(psip))& 
          call polymark(real(omegap/tqpsi),imag(omegap/uqpsi),1,3)
!     call color(15)
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
  if(vsin.eq.9999.)then
!     call legendline(-0.15,1.05,258,'No ions')
  else
     call legendline(-.2,1.05,258,'v!ds!d='//string)
  endif
  call fwrite(omegacp,iwidth,2,string)
  call legendline(.1,1.04,258,'!AW!@/!A)y!@='//string)
  call fwrite(kin,iwidth,2,string)
  call legendline(.43,1.04,258,'k/!A)y!@='//string)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine printmatrix(psip)
  use shiftgen
  complex :: Ftotalg,Cfactor
  complex :: w2byw4,w2byw41,wqbyw4       !,w2byw42,w2byw43
  complex :: mdofzd(-nzd:nzd,nmdmax),mdtotal(-nzd:nzd)
  character, dimension(nmdmax) :: mname=(['4','2','q'])
  w2byw4=0.
  write(*,'(a,f7.4,a,2f8.5,a,f7.4,a,f8.3)')&
       ' kg=',kg,' omega=(',omegag,') psi=',psip,' Omegac/sqrt(psi)&
       &=',Omegacg/sqrt(psip) ,' kpar=',kpar
  Ftotalg=Fpg(0)
  qresdenom=4.*kpar*(1/real(omegag)**2-1.)
 ! These are total forces integrated dW.
!  Cfactor=1.+sqm1g*3.1415926*2*(Fintqq-Fintqw+Fextqqwanal)/(qresdenom/16.)
  Cfactor=1.+sqm1g*3.1415926*2*(Fintqq-Fintqw+Fextqqwanal)/(qresdenom/16.)
!  if(nharmonicsg.gt.0)then
     write(*,*)'Nharmonics=',nharmonicsg,'  Harmonic sum values:'
     write(*,'(a,8f8.5)')' <4|V|4>  (Ftotalsum)=',Ftotalsumg/f4norm**2
     write(*,*)'Ftmdsum  <4|                  <2|                   <q|'
     write(*,'('' ('',2f11.6,'') ('',2f11.6,'') ('',2f11.6,'')'')')Ftmdsum
     write(*,*)'dispM    <4|                  <2|                   <q|'
     write(*,'('' ('',2f11.6,'') ('',2f11.6,'') ('',2f11.6,'')'')')dispM
!  endif
  write(*,'(a,2g14.6)')'Determinant=',dispMdet
  w2byw41=(dispM(1,3)*dispM(3,1)-dispM(1,1)*dispM(3,3))/&
       (dispM(1,2)*dispM(3,3)-dispM(1,3)*dispM(3,2))
  wqbyw4=-(dispM(1,1)*dispM(3,2)-dispM(3,1)*dispM(1,2))/&
       (dispM(1,3)*dispM(3,2)-dispM(3,3)*dispM(1,2))
  write(*,'(a,2f10.6,a,2f10.6)')' a2/a4=',amd(2),'   aq/a4=',amd(3)
  write(*,*)'M13*M31/M33',dispM(1,3)*dispM(3,1)/dispM(3,3)
  
  mdtotal=0.
  mdofzd=0.
  do i=-nzd,nzd
     do j=1,nmd
        mdofzd(i,j)=amd(j)*mdofz(zdent(i),j,kpar)
        mdtotal(i)=mdtotal(i)+mdofzd(i,j)
     enddo     
  enddo
  call minmax(real(mdtotal),2*nzd+1,rmin,rmax)
  rmax=.446/rmax  ! Scale to make height same as |4>
  open(15,file='totalmode.plt',status='unknown')
  write(15,*)'Output of ocontgen, z, real, imag'
!      write(15,'(a)')':-rx-20,20'
!      write(15,'(a)')':-lo.05,.05'
      write(15,'(a)')'legend:real(total)'
      write(15,'(a)')'legend:imag(total)'
!      write(15,'(a)')'legend:shiftmode |4>'
      write(15,*)2*nzd+1,2
      write(15,'(3f10.5)')(zdent(i),real(mdtotal(i)*rmax)&
           &,imag(mdtotal(i)*rmax) ,i=-nzd,nzd)
      close(15)

  is=-20
  ie=int(.5*nzd)
!  call autoinit(zdent(is:ie),real(mdofzd(is:ie,2)),ie+1-is)
  call multiframe(2,1,0)
  call dcharsize(.02,.02)
  call minmax(mdofzd(is:ie,2),2*(ie+1-is),amin,amax)
  call minmax(mdofzd(is:ie,3),2*(ie+1-is),bmin,bmax)
  call pltinit(zdent(is),zdent(ie),min(amin,bmin),max(amax,bmax))
  call axis
  call axis2
  call axlabels('','|u>')
  call legendline(.7,.2,0,' Real')
  call dashset(2)
  call legendline(.7,.1,0,' Imag')
  call dashset(0)
  call winset(.true.)
  do j=1,nmd
     call color(2*j-1)
     call labeline(zdent(is:ie),real(mdofzd(is:ie,j)),ie+1-is,mname(j),1)
     call dashset(2)
     call polyline(zdent(is:ie),imag(mdofzd(is:ie,j)),ie+1-is)
     call dashset(0)
  enddo
  call color(15)
  call minmax(mdofzd(is:ie,1),2*(ie+1-is),amin,amax)
  call pltinit(zdent(is),zdent(ie),amin,amax)
  call axis
  call axis2
  call axlabels('z','|u>')
  call color(4)
  call labeline(zdent(is:ie),real(mdtotal(is:ie)),ie+1-is,'total',5)
  call dashset(4)
  call color(1)
  call labeline(zdent(is:ie),real(mdofzd(is:ie,1)),ie+1-is,mname(1),1)
  call dashset(0)
  call color(15)
  call pltend()
  call multiframe(0,0,0)

end subroutine printmatrix
