! Continuous and discrete perturbation modes for sech^4x equilibrium
program continuous
  integer, parameter :: nx=100,npts=2*nx+1
  real, dimension(-nx:nx) :: x,Pk,Qk,T,S
  real, dimension(-nx:nx,5) :: dmodes
  real :: xmax=5
  real :: k=.02
  character*50 string
  
  do i=-nx,nx
     x(i)=xmax*i/nx
  enddo
  S=1./cosh(x)
  T=tanh(x)

  call pfset(3)
  call multiframe(2,1,1)
  call dcharsize(.018,.018)
  nk=2
  do j=1,nk
     Pk= -15*(k**4 + (28*S**2 - 15)*k**2 + 63*S**4 - 56*S**2 + 8)
     Qk= k**4 + (105*S**2 - 85)*k*2 + 945*S**4 - 1155*S**2 + 274
! Continuum normalization factor to make a delta function.
     fnorm=1.
     do l=1,5
        fnorm=fnorm*(k**2+l**2)
     enddo
     fnorm=sqrt(fnorm*(2.*3.1415926))
     if(j.eq.1)then
        call autoinit(x,(T*Pk*cos(k*x)-k*Qk*sin(k*x))/fnorm,npts)
        call axis
     endif
     call color(j)
     call fwrite(k,iwidth,2,string)
! antisymmetric  
     call iwrite(j,iwidth,string)
     call polyline(x,(T*Pk*cos(k*x)-k*Qk*sin(k*x))/fnorm,npts)
     call legendline(.6,1.-.1*j,258,'Continuum p='//string(1:iwidth))
! symmetric
     call dashset(2)
     call polyline(x,(T*Pk*sin(k*x)+k*Qk*cos(k*x))/fnorm,npts)
     call dashset(0)
     k=k*10.
  enddo
  call color(15)
  call polyline([-xmax,xmax],[0.,0.],2)

! Discrete  
!  dmodes(:,1)=45*S*(21*S**4-28*S**2+9)
!  dmodes(:,3)=-105*S**3*(9*S**2-8)
!  dmodes(:,5)=945*S**5
!  dmodes(:,2)=315*T*S**2*(3*S**2-2)
!  dmodes(:,4)=-945*S**4*T
! Normalized:
  dmodes(:,1)=S*(21*S**4-28*S**2+8)*sqrt(30.)/16
  dmodes(:,3)=-S**3*(9*S**2-8)*sqrt(105.)/16
  dmodes(:,5)=S**5*3*sqrt(35.)/16
  dmodes(:,2)=T*S**2*(3*S**2-2)*sqrt(105.)/4
  dmodes(:,4)=-S**4*T*3*sqrt(70.)/8

  call autoinit(x,dmodes(:,4),npts)
  call axis
  call axlabels('x','')
  call labeline(x,dmodes(:,4),npts,'4',1)
  call labeline(x,dmodes(:,2),npts,'2',1)
  call legendline(.7,.2,258,'Discrete')
  call legendline(.05,.2,0,' antisymmetric')
  call dashset(2)
  il=nint(nx*.25)
  write(*,*)'il=',il,x(il)
  call polyline(x,dmodes(:,5),npts)
  call jdrwstr(wx2nx(x(il)),wy2ny(dmodes(il,5)),'5',0.)
  call polyline(x,dmodes(:,3),npts)
  call jdrwstr(wx2nx(x(il)),wy2ny(dmodes(il,3)),'3',0.)
  call polyline(x,dmodes(:,1),npts)
  call jdrwstr(wx2nx(x(il)),wy2ny(dmodes(il,1)),'1',-.5)
  call legendline(.05,.1,0,' symmetric')
  call pltend
  
end program continuous




