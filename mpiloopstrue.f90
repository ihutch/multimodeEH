module mpiloops
! Routines to parallelize loops. Essentially self-explanatory.
! All routines are abstracted to isolate mpi calls to this module.
! It can be replaced with mpiloopsdummy.f90 if mpi linking is unavailable.
! Alternatively a call to mpilserial forces serial operation but should
! not be used in a multiprocessor environment (i.e. under mpiexec.)
  integer, private :: nprcsses,myid,ierr
  logical, private :: lmpinitd=.false.,lserial=.false.
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilprep(theid,thenproc)
! Normali mpi intialization call Can be called earlier in the host
! program if preventing slaves from action is needed via a
! test. Multiple calls are safe.
    integer :: theid,thenproc
    if(lserial)return
    call mpilgetmyid(myid,nprcsses,ierr)
    theid=myid
    thenproc=nprcsses
  end subroutine mpilprep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilgetmyid(myid,nprcsses,ierr)
! If necessary initialize the MPI system.
! Get my MPI id, and the number of processors.
    use mpi
    if(lserial)return
    call MPI_INITIALIZED(lmpinitd,ierr)
    if(.not.lmpinitd) call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcsses, ierr )
  end subroutine mpilgetmyid
!********************************************************************
  subroutine mpilbarrier(ierr)
    use mpi
    if(lserial)return
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine mpilbarrier
!*******************************************************************
  subroutine mpilcommsreal(buf,iactiv,nvals,tag)
    use mpi
    real buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    if(lserial)return
    datatype=MPI_REAL
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from active slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommsreal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilcommsinteger(buf,iactiv,nvals,tag)
    use mpi
    integer buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    if(lserial)return
    datatype=MPI_INTEGER
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from active slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommsinteger
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilcommscomplex(buf,iactiv,nvals,tag)
    use mpi
    complex buf(*)
    integer iactiv,nvals,datatype,dest,tag,comm,ierr,stat(MPI_STATUS_SIZE)
    if(lserial)return
    datatype=MPI_COMPLEX
    dest=0
    comm=MPI_COMM_WORLD
! Communicate the results to master if necessary.
    if(.not.iactiv.eq.0)then ! Only send/recv from active slaves.
       if(myid.eq.iactiv)call MPI_SEND(buf,nvals,datatype,dest,tag,comm,ierr) 
       if(myid.eq.0)call MPI_RECV(buf,nvals,datatype,iactiv,tag,comm,stat,ierr)
    endif
  end subroutine mpilcommscomplex
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilkillslaves
    use mpi
! Kills the slaves by calling mpifinalize and exiting them.
! Then no more parallelism works. Avoids hang with mpistopslaves
    if(lserial)return
    call MPI_FINALIZE(ierr)
    if(myid.ne.0) call exit
  end subroutine mpilkillslaves
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilstopslaves
    use mpi
! Stops the slaves by calling a barrier for them but not master.
! If called, nothing further is done by slaves; then they crash.
! Or if mpilfreeslaves is called, they continue (little point in that).
    if(lserial)return
    if(myid.ne.0)call MPI_BARRIER(MPI_COMM_WORLD, IERR)
  end subroutine mpilstopslaves
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine mpilfreeslaves
    use mpi
! Releases the slaves if they have been stopped by mpilstopslaves.
    if(lserial)return
    if(myid.eq.0)call MPI_BARRIER(MPI_COMM_WORLD, IERR)
  end subroutine mpilfreeslaves
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine mpilserial
! Force serial operation to speed up initialization not when using mpiexec.
   write(*,'(a,i4,i3)')'Forcing Serial Running'
   lserial=.true.
   nprcsses=1
   myid=0
 end subroutine mpilserial
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
  nproc=0 ! Sometimes needed to silence warnings.
  impi=0
  call mpilprep(id,nproc) ! Return id: my rank, nproc: the total rank.
  do j=1,nj
   do i=1,ni
     iactiv=mod(impi,nproc)  ! Decide the active rank process iactiv
     if(iactiv.eq.id)then    ! If I am active I do the work needed ...
        sbuf(impi)=complex(impi,1./impi)  ! Trivial example.
        write(*,*)'Process',id,' Set the value',sbuf(impi)
     endif
     call mpilcommscomplex(sbuf(impi),iactiv,nvals,impi)
!     call mpilbarrier(ierr)   ! Not necessary because comms should block.
     impi=impi+1
   enddo
  enddo
  call mpilstopslaves        ! Prevent slaves from continuing, usually desired.
! call mpilkillslaves        ! Alternative that kills the slave processes.
  write(*,*)'Finished:'
end subroutine mpiltemplate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Optional main for testing etc.
! call mpiltemplate
! end
