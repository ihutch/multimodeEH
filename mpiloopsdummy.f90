module mpiloops
! Dummy routines to disable parallelization via mpiloops. 
! It should be used in place of module mpiloops to produce a totally 
! serial executable without reference to any mpi routines.
contains
  subroutine mpilbarrier(ierr)
    ierr=0
  end subroutine mpilbarrier
  subroutine mpilgetmyid(myid,nprcsses,ierr)
    myid=0;nprcsses=1;ierr=0
  end subroutine mpilgetmyid
  subroutine mpilprep(id,nproc)
    id=0;nproc=1
  end subroutine mpilprep
  subroutine mpilcommsreal(buf,iactiv,nvals,itag)
    real buf(*),b
    b=buf(1);i=iactiv;i=nvals;i=itag  ! Silence warnings.
  end subroutine mpilcommsreal
  subroutine mpilcommscomplex(buf,iactiv,nvals,itag)
    complex buf(*),b
    b=buf(1);i=iactiv;i=nvals;i=itag
  end subroutine mpilcommscomplex
  subroutine mpilcommsinteger(buf,iactiv,nvals,itag)
    integer buf(*),b
    b=buf(1);i=iactiv;i=nvals;i=itag
  end subroutine mpilcommsinteger
  subroutine mpilstopslaves
  end subroutine mpilstopslaves
  subroutine mpilfreeslaves
  end subroutine mpilfreeslaves
end module mpiloops
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This is a template for using mpiloops to parallelize a loop.
! Parallelizing nested loops benefits from serializing impi and
! parallelizing the nest if nproc much exceeds the length of the inner
! loop(s). 
subroutine mpiltemplate
  use mpiloops
  integer, parameter :: ni=5,nj=5,length=ni*nj
  complex, dimension(length) :: sbuf  ! Result and communication buffer.
  nvals=1
  impi=0
  do j=1,nj
  do i=1,ni
     impi=impi+1
     call mpilprep(id,nproc) ! Return id my rank, nproc the total rank.
     iactiv=mod(impi,nproc)  ! Decide the active rank process iactiv
     if(iactiv.eq.id)then    ! If I am active I do the work needed ...
        sbuf(impi)=complex(impi,1./impi)  ! Trivial example.
        write(*,*)'Process',id,' Set the value',sbuf(impi)
     endif
     call mpilcommscomplex(sbuf(impi),iactiv,nvals,impi)
!     call mpilbarrier(ierr)   ! Not necessary because comms should block.
  enddo
  enddo
  call mpilstopslaves        ! Prevent slaves from continuing, usually desired.
  write(*,*)'Finished:'
  write(*,*)sbuf
end subroutine mpiltemplate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Optional main
!call mpiltemplate
!end
