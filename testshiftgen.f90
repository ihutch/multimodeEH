! Version of testshiftgen with all shiftmode dependence purged out. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testFrepel
    use shiftgen
    complex :: Ftotalg
    character*100 annote,tban
    write(*,*)'testFrepel'
    if(real(omegag).eq.0.)omegag=(.1,0.00000)
    write(*,*)omegag
    omegaonly=omegag
    if(psig.le.0)psig=.5
    isigma=-1
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'', '//&
         ' ''  v!ds!d='',f5.3)')psig,real(omegag),imag(omegag),vshift
    call FgRepelEint(Ftotalg,isigma)
    tbmax=tbr(nge-2)*0.9
    tbfac=5.*nint(tbmax/Wgarrayp(nge)/5.)
    write(*,*)'Repelling Ftotalpg',Ftotalpg
    write(*,*)'Repelling Ftotalrg',Ftotalrg
    write(*,*)'2*Sum             ',2*(Ftotalpg+Ftotalrg)
    write(*,*)'Repelling Ftotalg ',Ftotalg
    write(*,*)'R(Ftotalg)/psi^2  ',real(Ftotalg)/psig**2,128./315,'=128/315'
    call dcharsize(.03,.03)
    call multiframe(2,2,0)
    call pltinit(vinfarrayr(nge),vinfarrayr(1),0.0005,Wgarrayp(nge))
!       call axlabels('v!d!A;!@!d','W')
    call axis
    call axlabels('','W!d!A|!@!d')
    call legendline(0.5,1.08,258,annote(1:lentrim(annote)))
    call winset(.true.)
    call polyline(vinfarrayr,Wgarrayr,nge)
!    call polymark(vinfarrayr,Wgarrayr,nge,ichar('|'))
    call dashset(2)
    call color(5)
    call fwrite(tbfac,iwidth,0,tban)
    call polyline(vinfarrayr,tbr/tbfac,nge)
    call legendline(.5,.3,0,' !Bt!@!dorbit!d/'//tban(1:iwidth))
    call dashset(0)
    call color(15)
    call polymark(vinfarrayr,Wgarrayp(nge)*.98+1.e-6*Wgarrayp,nge,ichar('|'))
!    call pltend
    call minmax(Forcegp,2*nge,pmin,pmax)
    call minmax(Forcegr,2*nge,rmin,rmax)
    call pltinit(vinfarrayr(nge),vinfarrayr(1),min(pmin,rmin),max(pmax,rmax))
    call axis
    call axlabels('v!d!A;!@!d','dF/dv!d!a;!@!d')
    call legendline(.3,.9,258,'Reflected')
    call color(1)
    call polyline(vinfarrayr,real(Forcegr),nge)
    call legendline(.2,.16,0,' real')
    call color(2)
    call dashset(2)
    call polyline(vinfarrayr,imag(Forcegr),nge)
    call legendline(.2,.08,0,' imag')
    call dashset(0)
    call color(15)
    call pltinit(vinfarrayp(1),vinfarrayp(nge),0.,Wgarrayp(nge))
    call axis
    call axis2
    call winset(.true.)
    call polymark(vinfarrayp,Wgarrayp(nge)*.98+1.e-6*Wgarrayp,nge,ichar('|'))
    call polyline(vinfarrayp,Wgarrayp,nge)
    call dashset(2)
    call color(5)
    call polyline(vinfarrayp,tbp/tbfac,nge)
    call dashset(0)
    call color(15)
    call pltinit(vinfarrayp(1),vinfarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call legendline(.3,.9,258,'Passing')
    call axis
    call axis2
    call axlabels('v!d!A;!@!d','')
    call color(1)
    call polyline(vinfarrayp,real(Forcegp),nge)
    call color(2)
    call dashset(2)
    call polyline(vinfarrayp,imag(Forcegp),nge)
    call dashset(0)
    call pltend
    call multiframe(0,0,0)
  end subroutine testFrepel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testAttract(isw)
    use shiftgen
    real, dimension(nge) :: vpsiarrayp
    complex :: Ftotalg
    character*40 annote,ffan,cij
    character*20, dimension(2,3) :: modelabel
! Give access to pf switch pfsw    
    integer pfsw
    integer pfilno,pfnextsw,pfPS,psini
    common/pltfil/pfsw,pfilno,pfnextsw,pfPS,psini
    data modelabel/'<2|!p!o~!o!qV|4>','<q|!p!o~!o!qV|4>'&
         ,'<4|!p!o~!o!qV|2>','<4!p!o~!o!q|V|q>'&
         ,'<2|!p!o~!o!qV|2>','<q|!p!o~!o!qV|q>'/

    write(*,*)'testAttract'
    lioncorrect=.false.
    if(real(omegag).eq.0)omegag=(.1,0.001000)
    omegaonly=omegag
    kg=real(omegag)
    if(psig.ge.0)psig=-.5
    isigma=-1    
    if(vshift.ne.0.)write(*,*)'WARNING testAttract with vshift=',vshift
!    write(*,*)'Entered testattract'
    write(annote,'(''!Ay!@='',f5.3,'' !Aw!@=('',f5.3'','',f5.3,'')'')')&
         psig,real(omegag),imag(omegag)
    if(pfnextsw.ge.0)then
       call dcharsize(.02,.02)
    else
       call dcharsize(.025,.025)
    endif
    call FgEint(Ftotalg,isigma)  ! Generic call is the same.
    if(isw.eq.2.and..false.)then
       write(*,'(a,3f8.5)')'     W               tb/2       for  psi=',-psig
       write(*,'(i4,g12.4,f12.4)')(i,Wgarrayr(i),tbr(i),i=1,nge)
    endif
    vpsiarrayp=sqrt(2.*(Wgarrayp(1:nge)-psig))
    psi=-psig                     ! psi is the positive depth
    omegac=10.
    write(*,*)'omega=',omegag
    write(*,*)'Ftotalg        ',Ftotalg
    write(*,*)'R(Ftotalg)/psi^2  ',real(Ftotalg)/psig**2,-.14776!,-0.5*.14776
    call multiframe(1,2,0)
    call minmax(Forcegp,2*nge,pmin,pmax)
    call minmax(Forcegr,2*nge,rmin,rmax)
!    fpfac=max(5.*int(min(abs(rmax/pmax),abs(rmin/pmin))/5.),1.)
    fpfac=max(10.*int(min(abs(rmax/pmax),abs(rmin/pmin))/10.),1.)
    call fwrite(fpfac,iwidth,0,ffan)
    call pltinit(vinfarrayr(nge),vinfarrayr(1)*1.01,min(pmin,rmin),max(pmax,rmax))
    call axis; call axis2
    call axlabels('v!d!Ay!@!d','dF/dv!d!Ay!@!d')
    call legendline(0.5,1.03,258,annote(1:lentrim(annote)))
    call legendline(0.3,.9,258,'Trapped')
    call winset(.true.)
!    call polymark(vinfarrayr,(max(pmax,rmax)*.97+vinfarrayr*1.e-7),nge&
!         &,ichar('|'))
    call color(1)
    call polyline(vinfarrayr,real(Forcegr),nge)
    call legendline(.05,.1,0,' real')
    call color(2)
    call dashset(2)
    call polyline(vinfarrayr,imag(Forcegr),nge)
    call legendline(.05,.05,0,' imag')
    call color(15)
    call dashset(0)
    call pltinit(vpsiarrayp(1),vpsiarrayp(nge),min(pmin,rmin),max(pmax,rmax))
    call axis; call axis2 
    call axlabels('v!d!Ay!@!d','')
    call legendline(0.3,.9,258,'Passing')
    call winset(.true.)
!    call polymark(vpsiarrayp,(max(pmax,rmax)*.97+vpsiarrayp*1.e-7),nge&
!         &,ichar('|'))
    call color(3)
    call polyline(vpsiarrayp,fpfac*imag(Forcegp),nge)
    call legendline(.4,.1,0,' realx'//ffan(1:iwidth))
    call color(4)
    call dashset(2)
    call polyline(vpsiarrayp,fpfac*real(Forcegp),nge)
    call legendline(.4,.05,0,' imagx'//ffan(1:iwidth))
    call dashset(0)
    call pltend
    call multiframe(0,0,0)

    if(.true.)then
    call charsize(.018,.018)
    call pltinit(0.,sqrt(-psig),min(pmin,rmin),max(pmax,rmax))
    call axis
    call axlabels('(-W)!u1/2!u','dF/dv!d!Ay!@!d')
    call polymark(sqrt(-psig-vinfarrayr**2/2),(max(pmax,rmax)*.97&
         &+vpsiarrayp*1.e-7),nge ,ichar('|'))
    call polyline(sqrt(-psig-vinfarrayr**2/2),real(Forcegr),nge)
    call legendline(.7,.1,0,' real')
    call dashset(2)
    call polyline(sqrt(-psig-vinfarrayr**2/2),imag(Forcegr),nge)
    call legendline(.7,.2,0,' imag')
    call dashset(0)
    call pltend
    endif

    if(nmd.ge.2)then  ! Plot of auxforces.
       write(*,*)'f4norm=',f4norm,' kpar=',kpar
       write(*,'(a,8f8.4)')'Complex Ftotalg <4|V|4> =',Ftotalg
       write(*,'(a,8f8.4)')'Self-Adjoint verification ratios (2,q)',&
          abs(Ftmda(2,1)/Ftmda(1,2)),abs(Ftmda(3,1)/Ftmda(1,3))
      do j=1,3
         f4h=f4norm
         if(j.eq.3)f4h=1.
         do i=1,nmd-1
            k1=mod(j,2)*i+1
            k2=min(j-1,1)*i+1
       call multiframe(1,2,0)
       call minmax(FmdpofW(k1,k2,:)*f4h,2*nge,pmin,pmax)
       call minmax(Fmdres(k1,k2,:)*f4h,2*nge,rmin,rmax)
       ratio=min(abs(rmax/pmax),abs(rmin/pmin))
       step=min(10,int(ratio))
       fpfac=max(step*int(ratio/step),1.)
!       write(*,'(a,5f8.4)')'rmax/pmax,rmin/pmin,fpfac',rmax/pmax,rmin/pmin,fpfac
       call fwrite(fpfac,iwidth,0,ffan)
       call pltinit(vinfarrayr(nge),vinfarrayr(1)*1.01,rmin,rmax)
!       write(cij,'(''i='',i1,'' j='',i1)')i,j
!       call legendline(.1,.95,258,cij(1:lentrim(cij)))
       call legendline(.1,.95,258,modelabel(i,j))
       write(cij,'(''('',2f8.4,'')'')')Ftmdt(k1,k2)*f4h
       call legendline(.35,.95,258,cij(1:lentrim(cij)))
       call axis; call axis2
       call axlabels('v!d!Ay!@!d','dFaux/dv!d!Ay!@!d')
       call legendline(0.5,1.03,258,annote(1:lentrim(annote)))
       call legendline(0.3,.9,258,'Trapped')
       call winset(.true.)
!    call polymark(vinfarrayr,(max(pmax,rmax)*.97+vinfarrayr*1.e-7),nge&
!         &,ichar('|'))
       call color(1)
       call polyline(vinfarrayr,real(Fmdres(k1,k2,:)*f4h),nge)
       call legendline(.05,.1,0,' real')
       call color(2)
       call dashset(2)
       call polyline(vinfarrayr,imag(Fmdres(k1,k2,:)*f4h),nge)
       call legendline(.05,.05,0,' imag')
       call color(15)
       call dashset(0)
       call pltinit(vpsiarrayp(1),vpsiarrayp(nge),rmin,rmax)
       call axis; call axis2 
       call axlabels('v!d!Ay!@!d','')
       call legendline(0.3,.9,258,'Passing')
       write(cij,'(''('',2f8.4,'')'')')Ftmdp(k1,k2)*f4h
       call legendline(.35,.95,258,cij(1:lentrim(cij)))
       call winset(.true.)
       call color(3)
       call polyline(vpsiarrayp,imag(FmdpofW(k1,k2,:))*f4h*fpfac,nge)
       call legendline(.4,.1,0,' real*'//ffan(1:iwidth))
       call color(4)
       call dashset(2)
       call polyline(vpsiarrayp,real(FmdpofW(k1,k2,:))*f4h*fpfac,nge)
       call legendline(.4,.05,0,' imag*'//ffan(1:iwidth))
       call dashset(0)
       call pltend
       enddo
      enddo
              
    endif
    
    lioncorrect=.true.
  end subroutine testAttract
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Frepelofomega
    use shiftgen
    integer, parameter :: nor=100
    real, dimension(nor) :: or
    complex, dimension(nor) :: frcomplex,fion
    character*30 string
    write(*,*)'Frepelofomega'
!    complex :: Ftotalg
    kg=0.
    ormax=10.   !defaults
    oi=0.01
    psig=0.01
    vs=0.
    isigma=-1
    Omegacin=10.
    Fimmobile=(128./315.)       ! Now normalized *psig**2 
! Because frcomplex includes only one velocity direction.
    nvs=1
    write(*,*)'vshift,psig,nvs',vshift,psig,nvs
    psig=abs(psig)
    vsmax=vshift
    ol=.4
    nl=int(nor*min(ol/ormax,1.))
    do j=1,nvs
       if(nvs.gt.1)vshift=vsmax*(j-1.)/(nvs-1.)
       do i=1,nor
!       or(i)=ormax*(i-1.)/(nor-1.)
          or(i)=ormax*float(i)/(nor)
          omegag=complex(or(i),oi)
          omegaonly=omegag
          call FgRepelEint(frcomplex(i),isigma)
          if(.not.real(frcomplex(i)).lt.1.e20)then
             write(*,*)'Frepelofomega Force Nan?',&
                  i,omegag,omegaonly,frcomplex(i),psig,isigma
             stop
          endif
! Test of ionforce, which takes electron omega arguments, whereas omega
! here has been specified in ion units.
          vs=vshift
!          write(*,*)'ionforce call',omegag,rmime,psig,vs,isigma
          call ionforce(fion(i),omegag/sqrt(rmime),0.,Omegacin,psig,vs,isigma)
          fion(i)=fion(i)/psig**2
          frcomplex(i)=frcomplex(i)/psig**2
          diffmax=max(diffmax,abs(fion(i)-frcomplex(i)))
       enddo
       if(j.eq.1)then
          write(*,*)'diffmax=',diffmax
          write(*,*)'Fimmobile/2=',Fimmobile,' Fdirect=',frcomplex(nor)
          write(*,'(a,f9.5,a,f9.5)')' vshift=',vshift,' psig=',psig
          call pltinit(0.,or(nor),-0.8*Fimmobile,1.2*Fimmobile)
          call charsize(.02,.02)
          call axis
          call axptset(0.,1.)
          call ticrev
          call altxaxis(1./sqrt(rmime),1./sqrt(rmime))
          call legendline(.35,1.13,258,'real(!Aw!@)/!Aw!@!dpe!d')
          call axptset(1.,0.)
          call ticlabtog
          call altyaxis(1.,1.)
          call ticlabtog
          call axptset(0.,0.)
          call ticrev
          call axlabels('real(!Aw!@)/!Aw!@!dpi!d','!p!o~!o!qF!di!d/!Ay!@!u2!u')
          call polymark(ormax,Fimmobile,1,1)
       endif
       call color(mod(j-1,15)+1)
       call polyline(or,real(frcomplex),nor)
       call fwrite(vshift,iwidth,2,string)
       call jdrwstr(wx2nx(ormax*.95),wy2ny(real(frcomplex(nor))),string(1:iwidth),-1.)
       call jdrwstr(wx2nx(ol),wy2ny(imag(frcomplex(nl))),string(1:iwidth),0.)
       if(j.eq.1)call legendline(.5,.1,0,' real')
       call dashset(2)
       call polyline(or,imag(frcomplex),nor)
       if(j.eq.1)call legendline(.5,.15,0,' imag')
       call dashset(0)
    enddo
    call pltend
  end subroutine Frepelofomega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Frepelofimagomega(rmime)
! Version that does not "use" shiftgen, only encapsulated forces.
! Hence one does not have to worry about resetting things.
    integer, parameter :: nor=100
    real, dimension(nor) :: or
    complex, dimension(nor) :: fion,fion2,felec,felec2
    character*30 string
    complex :: omegain,Fe
!    real :: rmime=1836.
    write(*,*)'Frepelofiomagomega'
    ormax=10.   !defaults
    oi=0.01
    Omegacin=10.
    isigma=-1
    Fimmobile=(128./315.)       ! Now normalized *psig**2 
    nvs=1
    ps=.1
    vs=1.
    call tsparse(ormax,oi,nvs,isw,vs,ps) ! Get values from cmdline
    vsmax=vs
    ol=.4
    nl=int(nor*min(ol/ormax,1.))
    write(*,'(a,f9.5,a,f9.5)')' vsmax=',vs,' ps=',ps
    do j=1,nvs
       if(nvs.gt.1)vs=vsmax*(j-1.)/(nvs-1.)
       do i=1,nor
          or(i)=ormax*(i-.9)/(nor-.9)
          omegain=complex(oi,or(i))/sqrt(rmime)
          call ionforce(fion(i),omegain,0.,Omegacin,ps,vs,isigma)
          fion(i)=fion(i)/ps**2
          call ionforce(fion2(i),omegain,0.,Omegacin,ps/10.,vs,isigma)
          fion2(i)=fion2(i)/(ps/10.)**2
          if(j.eq.1)then
             call electronforce(Fe,omegain,0.,Omegacin,ps/10.,vs,isigma)
             felec2(i)=Fe/(ps/10.)**2
             call electronforce(Fe,omegain,0.,Omegacin,ps,vs,isigma)
             felec(i)=Fe/ps**2
          endif
       enddo
       if(j.eq.1)then
          call pltinit(0.,or(nor),-0.8*Fimmobile,1.2*Fimmobile)
          call charsize(.02,.02)
          call axis
          call axptset(0.,1.)
          call ticrev
          call altxaxis(1./sqrt(rmime),1./sqrt(rmime))
          call legendline(.35,1.13,258,'imag(!Aw!@)/!Aw!@!dpe!d')
          call axptset(1.,0.)
          call ticlabtog
          call altyaxis(1.,1.)
          call ticlabtog
          call axptset(0.,0.)
          call ticrev
          call axlabels('imag(!Aw!@)/!Aw!@!dpi!d','real(!p!o~!o!qF)&
               &/!Ay!@!u2!u')
          call legendline(.7,.05,0,' +!p!o~!o!qF!di!d')
          call dashset(2)
          call winset(.true.)
          call polyline(or,-real(felec),nor)
          call fwrite(ps,iwidth,3,string)
          call legendline(.1,.12,0,' -!p!o~!o!qF!de!d, !Ay!@='&
               &//string(1:iwidth))
          call dashset(3)
          call polyline(or,-real(felec2),nor)
          call fwrite(ps/10.,iwidth,3,string)
          call legendline(.1,.05,0,' -!p!o~!o!qF!de!d, !Ay!@='&
               &//string(1:iwidth))
          call dashset(0)
       endif
       call color(mod(j-1,15)+1)
       call polyline(or,real(fion),nor)
       call polyline(or,real(fion2),nor)
       call fwrite(vs,iwidth,2,string)
       call jdrwstr(wx2nx(ormax*.95),wy2ny(real(fion(nor)))&
            &,string(1:iwidth),-1.)
    enddo
    call pltend
  end subroutine Frepelofimagomega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testSumHarm
    use shiftgen
!    complex :: Ftotalg
! Defaults
    write(*,*)'testSumHarm'
    nvs=1
    kg=0.
    ormax=.1    
    ps=-.1
    vs=0.
    lbess=.true.
    call tsparse(ormax,oi,nvs,isw,vs,ps)
    lioncorrect=.false.
    if(oi.lt.0.00001)oi=.00001
    omegag=complex(ormax,oi)
    omegaonly=omegag
    write(*,*)'      psig                      omegag,             &
         & Omegacg','   lioncorrect'
    write(*,*)psig,omegag,Omegacg,lioncorrect
    if(kg.ne.0.)write(*,*)'kg=',kg
    isigma=-1
    call SumHarmonicsg(isigma)
    write(*,*)'FtotalSumg=',Ftotalsumg
    lioncorrect=.true.
  end subroutine testSumHarm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tsparse(ormax,oi,nvs,isw,vs,ps)
    use shiftgen
    character*30 argument
    do i=1,iargc()
       call getarg(i,argument)
       if(argument(1:2).eq.'-p')read(argument(3:),*)ps
       if(argument(1:2).eq.'-v')read(argument(3:),*)vs
       if(argument(1:2).eq.'-i')read(argument(3:),*)isigma
       if(argument(1:3).eq.'-zm')read(argument(4:),*)zm
       if(argument(1:3).eq.'-or')read(argument(4:),*)ormax
       if(argument(1:3).eq.'-oi')read(argument(4:),*)oi
       if(argument(1:3).eq.'-oc')read(argument(4:),*)Omegacg
       if(argument(1:3).eq.'-kg')read(argument(4:),*)kg
       if(argument(1:2).eq.'-n')read(argument(3:),*)nvs
       if(argument(1:2).eq.'-a')read(argument(3:),*)nmd
       if(argument(1:2).eq.'-c')then
          ipfset=-3
          read(argument(3:),*,err=201,end=201)ipfset
201       call pfset(ipfset)
       endif
       if(argument(1:2).eq.'-s')then
          if(lentrim(argument).eq.2)then
             isw=0
          else
             read(argument(3:),*)j
             isw=isw+j
          endif
       endif
       if(argument(1:2).eq.'-h')goto 1
    enddo
    if(isw.eq.0)isw=1
    vshift=vs
    psig=ps
    return
1   continue
    write(*,*)' Usage: testshiftgen [-p,-v,-i,-zm,-or,-oi,-oc,-kg,-s,-n -h]'
    write(*,101) '-p..         set psi         [',psig
    write(*,101) '-v..         set vshift(max) [',vshift
    write(*,101) '-zm..        set zmax        [',zm
    write(*,101) '-or.. -oi..  set omega       [',ormax,oi
    write(*,101) '-oc.. -kg..  set O_c, set k  [',Omegacg,kg
    write(*,102) '-n..         set number of vs[',nvs
    write(*,102) '-a..         set nmodes      [',nmd
    write(*,102) '-s..         set switches    [',isw
    write(*,*)' s=1 Frepel, 2 Fattr, 4 Fi(omega), 8 PlotForce, 16&
         & denem,',' 32 Modes, 64 SumHarm'
    write(*,*)' s=128 fvinfplot, 256 Frepelofimagomega, 512 plotdent'
!    write(*,'(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,i3,a,f8.3,a,f8.3)') 'psi&
!         &=',psig,' zm=',zm,' omax=',ormax,' v=',vshift,' nv=',nvs ,'&
!         & oc=',Omegacg,' kg=',
    101 format(a,6f8.4)
    102 format(a,i5)
  end subroutine tsparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotionforce(psi,Typ,vsin)
  use shiftgen
  real :: psi,Typ,vsin
  integer, parameter :: nfi=10
  complex, dimension(nfi) :: Fiarray,omegaFi
  write(*,*)'plotionforce'
  omegamax=10
  nmd=1
  Omegacg=20
  isigma=1
  write(*,*)'psi=',psi,' vsin=',vsin,' isigma=',isigma
  write(*,*)'  i           omega                       Fi'
  do i=1,nfi
     omegaFi(i)=omegamax*(float(i)/nfi)/sqrt(rmime)+complex(0.,.01)/sqrt(rmime)
!         ionforce(Fi,omega,kin,Omegacin,psiin,vsin,isigma)
     call ionforce(Fiarray(i),omegaFi(i),0.1,5.,psi,vsin,isigma)
     Fiarray(i)=Fiarray(i)/psi
     write(*,'(i4,''  ('',2e11.3,'')'','' ('',2e11.3,'')'')')&
          i,omegaFi(i),Fiarray(i)
  enddo
  call minmax(Fiarray,2*nfi,fmin,fmax)
  call pltinit(0.,omegamax/sqrt(rmime),fmin,fmax)
  call axis
  call axlabels('omega','Fi')
  call polyline(real(omegaFi),real(Fiarray),nfi)
  call polyline(real(omegaFi),imag(Fiarray),nfi)
  call pltend
end subroutine plotionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testdenem
  use shiftgen
  integer, parameter :: nvs=100
  real, dimension(nvs) :: vs,denem
  denem=1.
  psig=.1
  vsmax=2.
  write(*,*)'testdenem'
  do i=1,nvs
     vrshift=i*vsmax/nvs
     vs(i)=vrshift
     call dfefac(denem(i),psig,vrshift)
  enddo
  call autoplot(vs,denem,nvs)
  call axlabels('vshift','fefac')
  write(*,'(''Factor by which dfe trapped is multiplied for psi='',f8.3)')psig
  call pltend
end subroutine testdenem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dfefac(dfeval,psig,vrshift)  ! Testing only. Inline in module.
! Correct the dfedW value approximately for the alteration of trapped
! equilibrium distribution caused by ion density perturbation.
    real, parameter :: ri=1.4
    vsx=1.3+0.2*psig
    vsa=vrshift
    denem1=(-1.+(vsa/vsx)**ri) /(1.+0.25*psig+vsa**2*(vsa/(vsx +(3.3&
         &/vsa)**1.5))**ri)
    dfeval=dfeval*(1-denem1)
  end subroutine dfefac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotmodes
  use shiftgen
!  kpar=0.01
  omegag=complex(.05,.001)
  call tsparse(ormax,oi,nvs,isw,vs,ps)
  kg=.1
  Omegacg=20
  zm=25
  psig=-.1
  isigma=-1
  Wg=.1*abs(psig)
  omegaonly=omegag   ! Needed for evaluation of kpar
  call makezg(isigma)
  write(*,'(i3,3f10.5)')(i,zg(i),pmds(i,3),i=-1,1)
  write(*,*)'kpar=',kpar,'omegag=',omegag
  dot0=0;dot1=0;dot2=0;dot12=0;dot01=0;dot02=0
  do i=-ngz+1,ngz
     dot0=dot0+(phigprime(i-1)**2+phigprime(i)**2)*(zg(i)-zg(i-1))/2
     ! Fix the z/x normalization by dividing by 2*4
     dot1=dot1+real(pmds(i-1,2)**2+pmds(i,2)**2)*(zg(i)-zg(i-1))/2
     dot2=dot2+(abs(pmds(i-1,3))**2+abs(pmds(i,3))**2)*(zg(i)-zg(i-1))/8
     dot12=dot12+real(pmds(i-1,2)*pmds(i-1,3)&
          +pmds(i,2)*pmds(i,3))*(zg(i)-zg(i-1))/8
     dot01=dot01+real(pmds(i-1,2)*phigprime(i-1)&
          +pmds(i,2)*phigprime(i))*(zg(i)-zg(i-1))/8
     dot02=dot02+real(phigprime(i-1)*pmds(i-1,3)&
          +phigprime(i)*pmds(i,3))*(zg(i)-zg(i-1))/8
  enddo

  ! Normalization of phigprime.
  ppnorm=3*sqrt(70.)/16
  dot0=dot0/(psig/ppnorm)**2
  dot01=dot01/(psig/ppnorm)
  dot02=dot02/(psig/ppnorm)
  write(*,*)'Integral of mode^2: should = 1,    1,   Uncomplete'
  write(*,*)'dot0 =',dot0,'  dot1 =',dot1,'  dot2 =',dot2
  write(*,*)'Overlap integrals: for zm=big,  should be zero'
  write(*,*)'dot12=',dot12,'  dot01=',dot01,'  dot02=',dot02
  write(*,*)'zm=',zm,' kpar=',kpar
  write(*,*)'psig/ppnorm',psig/ppnorm
  if(nmd.gt.1)call multiframe(nmd,1,1)
! call autoplot(zg,phigprime*ppnorm/psig,2*ngz+1) ! Normalized
  call autoplot(zg,phigprime,2*ngz+1)
  call axlabels('','d!Af!@/dz')
  do j=1,nmd-1
     call autoplot(zg,real(pmds(:,j+1)),2*ngz+1)
  enddo
  call axlabels('z','')
  call pltend

end subroutine plotmodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plotdent
  use shiftgen
  complex :: Ftotalg,Cfactor ! ,dfweight=999,cdummy
  complex :: F44t=0,F44p=0,w2byw4,wqbyw4
  character*10 string
  character(len=20), external :: ffwrite
  lioncorrect=.false.
  write(*,*)'To get density plots check out branch AuxmodeVersion'
  psidef=-.1
  if(psig.ge.0)psig=psidef
  omegacg=20
!  omegag=complex(.01,.001)
  omegag=complex(.047,.00119)
  omegaonly=omegag
  isigma=-1
  dentpass=0.
  kg=real(omegag)
  kg=.01
! Cross compare test
  omegacg=100
!  omegag=complex(1.814554073E-02,9.004140884E-05)
  omegag=complex(1.764957793E-02,8.126287139E-05)
  psig=-.25
  kg=.1
!  kg=.143
!  kpar=kg*real(omegag)/sqrt(1.-real(omegag))  ! Needed for makezdent.
  kpar=kg*real(omegaonly*sqrt((Omegacg**2+1-omegaonly**2)/&
       ((Omegacg**2-omegaonly**2)*(1-omegaonly**2))))

  write(*,'(a,f7.4,a,2f8.5,a,f7.4,a,f8.3)')&
       ' kg=',kg,' omega=(',omegag,') psi=',-psig,' Omegac=',Omegacg
  call makezdent
!  call FgEint(Ftotalg,isigma)
  call SumHarmonicsg(isigma) 
  Ftotalg=Fpg(0)
  qresdenom=4.*kpar*(1/real(omegag)**2-1.)
 ! These are total forces integrated dW.
!  Cfactor=1.+sqm1g*3.1415926*2*(Fintqq-Fintqw+Fextqqwanal)/(qresdenom/16.)
  Cfactor=1.+sqm1g*3.1415926*2*(Fintqq-Fintqw+Fextqqwanal)/(qresdenom/16.)

  if(nharmonicsg.gt.0)then
     write(*,*)'Nharmonics=',nharmonicsg,'  Harmonic sum values:'
     write(*,'(a,8f8.4)')' <4|V|4>  (Ftotalsum)=',Ftotalsumg/f4norm**2
     write(*,*)'Harm Sum <4|               <2|                <q|'
     write(*,'('' ('',2f8.4,'') ('',2f8.4,'') ('',2f8.4,'')'')')Ftmdsum
  endif
! Form the density versions of inner products, compensating for symmetry
! \int_-^+ phipd*(n4(+)-n4(-)) dz etc.
  dzd=zm/nzd
  do i=-nzd+1,nzd
     F44t=F44t-0.5*&
          (conjg(phipd(i))*(denttrap(i)-denttrap(-i))&
          +conjg(phipd(i-1))*(denttrap(i-1)-denttrap(1-i)))*dzd
     F44p=F44p-0.5*&
          (conjg(phipd(i))*(dentpass(i)-dentpass(-i))&
          +conjg(phipd(i-1))*(dentpass(i-1)-dentpass(1-i)))*dzd
  enddo
  write(*,*)
  write(*,*)'Details from just the m=0 contribution:'
  write(*,'(a,f8.4,a)')'Testshiftgen: Normalizing factor for |4>',f4norm,&
       ' applied to all forces.'
  write(*,*)'Force check:        Total             Passing          Trapped'
  write(*,'(a,8f8.4)')'Complex <4|V|4>  =',&
       Ftotalg/f4norm**2,2.*Ftotalpg/f4norm**2,2.*Ftotalrg/f4norm**2
  write(*,'(a,8f8.4)')'Complex <4|n_4   =',&
       (F44t+F44p)/f4norm,F44p/f4norm,F44t/f4norm
  write(*,'(a,2f8.4)')'Force nt4 verification ratio',&
       f4norm*F44t/(2.*Ftotalrg)
  write(*,*)'Inner product matrices:'
  write(*,*)'Passing <4|               <2|                <q|   [V-Vw|q>]'
  write(*,'('' ('',2f8.4,'') ('',2f8.4,'') ('',2f8.4,'')'')')2.*Ftmdp
  write(*,*)'Trapped <4|               <2|                <q|'
  write(*,'('' ('',2f8.4,'') ('',2f8.4,'') ('',2f8.4,'')'')')2.*Ftmdt     
  write(*,*)'Attract <4|               <2|                <q|'
  write(*,'('' ('',2f8.4,'') ('',2f8.4,'') ('',2f8.4,'')'')')Ftmda
  write(*,'(a,8f8.4)')'Self-Adjoint verification ratios (2,q)',&
       abs(Ftmda(2,1)/Ftmda(1,2)),abs(Ftmda(3,1)/Ftmda(1,3))
  write(*,*)
  w2byw4=Ftmda(2,1)/(kg**2+12/16.)
  write(*,*)'Amplitude         w_2/w_4=',w2byw4
  wqbyw4=-sqm1g*16*3.1415926*Ftmda(3,1)/(qresdenom*Cfactor)
  write(*,*)'Amplitude \int w_q dq/w_4=',wqbyw4          
  write(*,*)'         Coefficient magnitudes'
  write(*,'(a,2f8.4,a,f8.4)')'C= (',Cfactor,&
       ')    qresdenom=q0*(1./real(omegag)**2-1.)=',qresdenom
  write(*,'(a,$)')'Size 2 v 4 <2|V|4><4|V|2>/<4|V|4>(k^2-l2)='
  write(*,*)Ftmda(2,1)*Ftmda(1,2)/Ftotalg/(kg**2+12/16.)
  write(*,'(2a)')'Size of q term relative to 4 term,'
  write(*,*)'-i.4pi<q|V|4><4|V|q>/<4|V|4>/(qdenom*C/4)=',&
       -sqm1g*4*3.1415926*Ftmda(3,1)*Ftmda(1,3)/(qresdenom*Cfactor/4)*f4norm**2&
       /Ftotalg
  write(*,*)
  write(*,*)'    Cancellation:        <q|V|q>          <q|Vw|q>   &
       &       <q|V-Vw|q>'
  write(*,'(a,7f9.4)')'Complex qq internal:',Fintqq*2,Fintqw*2 &
       &,Fintqq*2-Fintqw*2
  write(*,'(a,2F9.4,a,2F9.4)')'<q|V-Vw|q> analytic, external',Fextqqwanal*2&
       ,'    total',(Fextqqwanal+Fintqq-Fintqw)*2
  write(*,'(a,2F9.4,a,2F9.4)')'<q|V-Vw|q> analytic, external',Fextqqwanal*2&
       ,'    Ftmda',(Ftmda(3,3))
  write(*,'(a,7f9.4)')'Complex <q|V-Vw|q>/<q|V|4>=' ,2.*&
       (Fintqq-Fintqw+Fextqqwanal)/(Ftmda(3,1))
  write(*,'(a,f7.4,a,2f8.5,a,f7.4,a,f8.5)')&
       ' kg=',kg,' omega=(',omegag,') psi=',-psig,' kpar=',kpar

  zm=12
  Wg=0.01
  call pfset(3)
  call makezg(isigma)
! |4>=-TS^4=+dphi/dz with psi=+1.
! Code says phigprimeofz=-psig*sinh(zval/4.)/cosh(zval/4.)**5 (psig -ve)
! so phigprime is -psi*dphi/dz and |4> is -phigprime/f4norm.     
  call multiframe(2,1,0)
  call dcharsize(.018,.018)
  call autoplot(zg,-phigprime/f4norm,2*ngz)
  call axis2
  call axlabels('','Perturbation Modes')
  ind2=int(ngz/9.)
  call jdrwstr(wx2nx(zg(ind2)),wy2ny(-phigprime(ind2)/f4norm),'4',0.)
  call dashset(5)
  call polyline(zg,-phigprime/f4norm+real(w2byw4*pmds(:,2)+wqbyw4*pmds(:,3)),2*ngz)
  call legendline(.7,.3,0,' real(total)')
  call dashset(0)
  call color(2)
  call polyline(zg,real(w2byw4*pmds(:,2)),2*ngz)
  call color(1)
  call polyline(zg,real(wqbyw4*pmds(:,3)),2*ngz)
  call dashset(2)
  call color(2)
  call polyline(zg,imag(w2byw4*pmds(:,2)),2*ngz)
  call jdrwstr(wx2nx(zg(ind2)),wy2ny(imag(w2byw4*pmds(ind2,2))),'2',0.)
  call color(1)
  indq=int(ngz*.9)
  call polyline(zg,imag(wqbyw4*pmds(:,3)),2*ngz)
  call jdrwstr(wx2nx(zg(indq)),wy2ny(real(wqbyw4*pmds(indq,3))),'q',0.)
  call color(15)
  call legendline(.7,.1,0,' imag')
  call dashset(0)
  call legendline(.7,.2,0,' real')
  call color(4)
  call dashset(2)
  call dashset(0)
  call color(15)
  xi=.02
  call autoplot(zg,-phig+xi*(-phigprime/f4norm+real(wqbyw4*pmds(:,3))&
       +real(w2byw4*pmds(:,2))),2*ngz)
  call axis2
  call fwrite(xi,iwidth,2,string)
!     write(*,*)trim(ffwrite(xi,iwidth,2)),'%'
  call axlabels('z','!AR!@(!Af!@) for amplitude w!d4!d='//string(1:iwidth))
  call legendline(.7,.8,0,'!Af!@!d0!d+!p!u~!u!q!af!@')
  call dashset(4)
  call polyline(zg,-phig,2*ngz)
  call legendline(.7,.7,0,'!Af!@!d0!d')
  call pltend     
end subroutine plotdent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=20) function ffwrite(x,iwidth,ipoint)
  call fwrite(x,iwidth,ipoint,ffwrite)
end function ffwrite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use shiftgen
integer :: nvs=1,isw=3
real :: ormax=0.,oi=0.
call tsparse(ormax,oi,nvs,isw,vs,ps)
iswin=isw
write(*,*)'isw=',isw
omegag=complex(ormax,oi)
  if(isw-2*(isw/2).eq.1) call testFrepel
  isw=isw/2 ! 2
  if(isw-2*(isw/2).eq.1) call testAttract(iswin)
  isw=isw/2 ! 4
  if(isw-2*(isw/2).eq.1) call Frepelofomega
  isw=isw/2 ! 8
  if(isw-2*(isw/2).eq.1) call plotionforce(.01,1.,0.)
  isw=isw/2 ! 16
  if(isw-2*(isw/2).eq.1) call testdenem
  isw=isw/2 ! 32
  if(isw-2*(isw/2).eq.1)  call plotmodes
  isw=isw/2 ! 64
  if(isw-2*(isw/2).eq.1) call testSumHarm
  isw=isw/2 ! 128
  if(isw-2*(isw/2).eq.1) call fvinfplot
  isw=isw/2 ! 256
  if(isw-2*(isw/2).eq.1) call Frepelofimagomega(1836.)
  isw=isw/2 ! 512
  if(isw-2*(isw/2).eq.1) call plotdent
end program
