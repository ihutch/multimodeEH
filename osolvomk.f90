! Solve for omega as a function of Omegac for various k,
! and specified psi.
program osolvomk
  use iterfind
  use mpiloops
  integer, parameter :: nomax=100,nkmax=20,nmds=3
  real, dimension(nomax,nkmax) :: or,oi,dor,doi
  complex, dimension(nmds,nomax,nkmax) :: amds
  complex, dimension(nomax,nkmax) :: rw1
  real, dimension(nomax) :: oc
  real, dimension(nkmax) :: kparray
  complex :: omegap
  complex, external :: rowone
  character*100 filename
  character*20 string
  integer :: isigma=-1,no=10,nk=1
  real :: kp=.1,psip=0.09,Ocmax=2.,vs=9999
  logical :: lg
  
  call parseos(psip,vs,kp,Ocmax,no,nk,lgetdet)
  call mpilprep(id,nproc)
  if(id.ne.0)litw=.false.
  write(filename,'(a,2i2.2,a,i3.3,a,i3.3,a,i3.3,a,i3.3,l1,a)')   &
       'F',no,nk, &
       'B',min(999,abs(nint(100*Ocmax))),   &
       'k',min(999,abs(nint(100*kp))),   &
       'v',min(999,abs(nint(100*vs))),  &
       'p',min(999,abs(nint(100*psip))),lgetdet,'.omk'

  if(id.eq.0)write(*,'(a,$)')'Attempt opening  '&
       &//filename(1:lentrim(filename))//' '
  open(12,file=filename,status='old',err=1)
  read(12,*)nof,ncol,nkf
  if(nof.ne.no.or.nkf.ne.nk)goto 1
  do i=1,nof
     read(12,*)oc(i),(or(i,j),oi(i,j),j=1,nk)
     read(12,*)oc(i),(dor(i,j),doi(i,j),j=1,nk)
     read(12,*)oc(i),(amds(2,i,j),amds(3,i,j),rw1(i,j),j=1,nk)
  enddo
  read(12,*)(kparray(j),j=1,nk)
  read(12,*,end=10)lg
  if(lg.neqv.lgetdet)goto 1
10 close(12)
  call mpilkillslaves
  write(*,*)' Succeeded.'
  goto 2        ! Just plot the read-in data.
1 continue     ! No or bad old file. Calculate and save.
  if(id.eq.0)close(12,status='delete')
  if(id.ne.0)close(12)
  if(id.eq.0)write(*,*)'  Failed. Calculating ...'

  nit=0
  impi=0
  oi=0;or=0.;rw1=0;dor=0;doi=0
  do j=1,nk
     kparray(j)=kp*j/nk
     thek=kparray(j)*sqrt(psip)
     do i=1,no
        if(j.eq.1)oc(i)=Ocmax*(i-.7)/(no-.7)
        Omegacp=sqrt(psip)*oc(i)
        if(nit.eq.0.or.oc(i).gt.0.3)omegap=complex(min(thek,oc(i)),thek/2)
        impi=impi+1
        iactiv=mod(impi,nproc)
        if(iactiv.eq.id)then
           call iterfindroot(psip,vs,Omegacp,omegap,thek,isigma,nit)
           if(nit.gt.1.and.nit.le.niterfind)then
              or(i,j)=real(omegap)/sqrt(psip)
              oi(i,j)=imag(omegap)/sqrt(psip)
              if(lgetdet)then
                 call getamd(nmds,amds(1:nmds,i,j))
                 rw1(i,j)=rowone()
              endif
              write(*,'(a,2i3,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a&
                &,i3)')'Found:',i,j,' Oc=' ,oc(i),' k/sqpsi=',kparray(j),'&
                & or=',or(i,j),' oi=',oi(i ,j),' id=',id
              if(.true.)then
                 lgetdet=.not.lgetdet
                 nit=0
                 call iterfindroot(psip,vs,Omegacp,omegap,thek,isigma,nit)
                 if(nit.gt.1.and.nit.le.niterfind)then
                    dor(i,j)=real(omegap)/sqrt(psip)-or(i,j)
                    doi(i,j)=imag(omegap)/sqrt(psip)-oi(i,j)
                    if(lgetdet)then
                       call getamd(nmds,amds(1:nmds,i,j))
                       rw1(i,j)=rowone()
                    endif
                    write(*,'(a,2i3,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a&
                         &,l1)')'Second:',i,j,' Oc=' ,oc(i),' k/sqpsi=',kparray(j),'&
                         & dor=',dor(i,j),' doi=',doi(i ,j),' lgetdet=',lgetdet
                 else
                    write(*,*)'Failed second iteration'
                 endif
                 lgetdet=.not.lgetdet
              endif
           endif
        endif
        call mpilcommsreal(or(i,j),iactiv,1,impi)
        call mpilcommsreal(oi(i,j),iactiv,1,impi)
        call mpilcommsreal(dor(i,j),iactiv,1,impi)
        call mpilcommsreal(doi(i,j),iactiv,1,impi)
        call mpilcommscomplex(amds(:,i,j),iactiv,nmds,impi)
        call mpilcommscomplex(rw1(i,j),iactiv,1,impi)
     enddo
  enddo
  call mpilkillslaves
  open(12,file=filename,status='new')
  write(12,*)no,2*nk,nk
  do i=1,no
     write(12,*)oc(i),(or(i,j),oi(i,j),j=1,nk)
     write(12,*)oc(i),(dor(i,j),doi(i,j),j=1,nk)
     write(12,*)oc(i),(amds(2,i,j),amds(3,i,j),rw1(i,j),j=1,nk)
  enddo
  write(12,*)(kparray(j),j=1,nk)
  write(12,*)lgetdet  
  close(12)

  
2 continue                  !Plot
!  do i=1,no
!     write(*,*)oc(i),(rw1(i,j),j=1,nk)
!  enddo

  call multiframe(2,1,3)
  call dcharsize(.02,.02)
  call minmax(or(1:no,1:nk),no*nk,rmin,rmax)
  call minmax(oi(1:no,1:nk),no*nk,gmin,gmax)
  omax=1.1*max(rmax,gmax)
  call pltinit(0.,oc(no),0.,omax)
  call axis;call axis2
!  call axlabels('!AW!@/!A)y!@','!Aw!@/!A)y!@')
  call axlabels('','!Aw!@/!A)y!@')
  call fwrite(psip,iwidth,2,string)
  call legendline(.2,1.05,258,'!Ay!@='//string(1:iwidth))
  call fwrite(vs,iwidth,2,string)
  if(vs.lt.9999)call legendline(.6,1.05,258,'v!ds!d='//string(1:iwidth))
  if(lgetdet)call legendline(.6,1.05,258,'Multimode')
  if(.not.lgetdet)call legendline(.6,1.05,258,'Shiftmode')
  call color(1)
  call jdrwstr(wx2nx(oc(no)*.7),wy2ny(max(.05*omax,or(no,1)-0.05&
       &*omax)),'real',0.)
  call color(4)
  call jdrwstr(wx2nx(oc(no)*.5),wy2ny(max(.05*omax,oi(no,1)+0.05&
       &*omax)),'imag',0.)
!  call winset(.true.) triggers dashed line plotting bug.
  do j=1,nk
     call dashset(j-1)
     call color(1)
     call polyline(oc,or(1,j),no)
     call color(4)
     call polyline(oc,oi(1,j),no)
     call color(15)
     call fwrite(kparray(j),iwidth,2,string)
     call legendline(.73,.97-j*.09,0,' '//string(1:iwidth))
  enddo
  call legendline(.73,.95,258,'k/!A)y!@=')
!  call pltend
  call minmax(doi(1:no,1:nk),no*nk,gmin,gmax)
  call pltinit(0.,oc(no),gmin,gmax)
  call axis;call axis2
  call axlabels('!AW!@/!A)y!@','!Adw!@/!A)y!@')
  if(lgetdet)call legendline(.6,1.05,258,'Shiftmode-Multimode')
  if(.not.lgetdet)call legendline(.6,1.05,258,'Multimode-Shiftmode')
  do j=1,nk
     call dashset(j-1)
     call color(1)
     call polyline(oc,dor(1,j),no)
     call color(4)
     call polyline(oc,doi(1,j),no)
     call color(15)
     call fwrite(kparray(j),iwidth,2,string)
!     call legendline(.73,.6-j*.09,0,' '//string(1:iwidth))
  enddo
  call dashset(0)
  call pltend
  
  call minmax(amds(2,1:no,1:nk),2*no*nk,gmin,gmax)
  call pltinit(0.,oc(no),gmin,gmax)
  call axis;call axis2
  call axlabels('','a!d2!d/a!d4!d')
  call fwrite(psip,iwidth,2,string)
  call legendline(.2,1.05,258,'!Ay!@='//string(1:iwidth))
  call color(1)
  call legendline(.1,.75,258,'real')
  call color(4)
  call legendline(.1,.9,258,'imag')
  do j=1,nk
     call dashset(j-1)
     call color(1)
     call polyline(oc,real(amds(2,1:no,j)),no)
     call color(4)
     call polyline(oc,imag(amds(2,1:no,j)),no)
     call color(15)
  enddo
  call dashset(0)
  
  call minmax(amds(3,1:no,1:nk),2*no*nk,gmin,gmax)
  call pltinit(0.,oc(no),gmin,gmax)
  call axis;call axis2
  call axlabels('!AW!@/!A)y!@','a!dq!d/a!d4!d')
  call legendline(.1,.8,258,'k/!A)y!@=')
  do j=1,nk
     call dashset(j-1)
     call color(1)
     call polyline(oc,real(amds(3,1:no,j)),no)
     call color(4)
     call polyline(oc,imag(amds(3,1:no,j)),no)
     call color(15)
     call fwrite(kparray(j),iwidth,2,string)
     call legendline(.1,.8-j*.09,0,' '//string(1:iwidth))
  enddo
  call pltend
end program osolvomk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parseos(psip,vsp,kp,Ocmax,no,nk,lgetdet)
  real :: kp
  character*50 :: argument
  integer :: ipset=0
  logical :: lgetdet
  do i=1,iargc()
     call getarg(i,argument)
     if(argument(1:2).eq.'-p')read(argument(3:),*)psip
     if(argument(1:2).eq.'-k')read(argument(3:),*)kp
     if(argument(1:2).eq.'-M')lgetdet=.not.lgetdet
     if(argument(1:3).eq.'-vs')read(argument(4:),*)vsp
!     if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
!     if(argument(1:3).eq.'-oi')read(argument(4:),*)oimax
     if(argument(1:3).eq.'-oc')read(argument(4:),*)Ocmax
     if(argument(1:3).eq.'-no')read(argument(4:),*)no
     if(argument(1:3).eq.'-nk')read(argument(4:),*)nk
     if(argument(1:3).eq.'-Ti')then
        read(argument(4:),*)Ti
        write(*,'(a,f8.2)')'Setting Ti to',Ti
        call Tset(1.,Ti)
     endif
     if(argument(1:2).eq.'-w')then
        if(lentrim(argument).gt.2)then
           read(argument(3:),*,err=2)ipset
           call pfset(ipset)
        endif
     endif
2    if(argument(1:2).eq.'-h')goto 1
  enddo
  
  return
1 continue
  write(*,'(a,f7.3)')'-p... set psip     [',psip
  write(*,'(a,f7.3)')'-k... set kp       [',kp
  write(*,'(a,l3  )')'-M toggle matrix   [',lgetdet
  write(*,'(a,f9.3)')'-vs.. set vs       [',vsp
  write(*,'(a,f7.3)')'-oc.. set Ocmax    [',Ocmax
  write(*,'(a,i5)')  '-no.. set no       [',no
  write(*,'(a,i5)')  '-nk.. set nk       [',nk
  write(*,'(a,i5)')  '-w... write plot?  [',ipset
  stop
end subroutine parseos
