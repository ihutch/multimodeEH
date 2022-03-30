module mpiloops
! Routines to parallelize loops
  integer, private :: nprcsses,myid,ierr
  logical, private :: lmpiflag=.false.  
contains
!********************************************************************
! All routines are abstracted to isolate mpi calls to this module.
  subroutine mpilbarrier(ierr)
    use mpi
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine mpilbarrier
!********************************************************************
  subroutine mpilgetmyid(myid,nprcsses,ierr)
! If necessary initialize the MPI system.
! Get my MPI id, and the number of processors.
    use mpi
    call MPI_INITIALIZED(lmpiflag,ierr)
    if(.not.lmpiflag) call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcsses, ierr )
  end subroutine mpilgetmyid
!*******************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilprep(theid,thenproc)
! Can be called earlier in the host program if preventing slaves from
! action is needed via a test. Multiple calls are safe.
    integer :: theid,thenproc
    call mpilgetmyid(myid,nprcsses,ierr)
    theid=myid
    thenproc=nprcsses
  end subroutine mpilprep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilcommsreal(buf,iactiv,nvals,tag)
    use mpi
    real buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    datatype=MPI_REAL
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommsreal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilcommsinteger(buf,iactiv,nvals,tag)
    use mpi
    integer buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    datatype=MPI_INTEGER
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommsinteger
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilcommscomplex(buf,iactiv,nvals,tag)
    use mpi
    complex buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    datatype=MPI_COMPLEX
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommscomplex
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilstopslaves
    use mpi
! Stops the slaves by calling a barrier for them but not master.
! If called, nothing further is done by slaves; then they crash.
! Or if mpilfreeslaves is called, they continue (little point in that).
    if(myid.ne.0)call MPI_BARRIER(MPI_COMM_WORLD, IERR)
  end subroutine mpilstopslaves
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilfreeslaves
    use mpi
! Releases the slaves if they have been stopped by mpilstopslaves.
    if(myid.eq.0)call MPI_BARRIER(MPI_COMM_WORLD, IERR)
  end subroutine mpilfreeslaves
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
! Optional main for testing etc.
! call mpiltemplate
! end
