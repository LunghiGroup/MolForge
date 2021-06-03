        module mpi_utils
        use mpi
        implicit none

         integer          :: mpi_id,mpi_nproc,err
         integer          :: mpi_blacs_id,mpi_blacs_nproc
         integer          :: mpi_blacs_world
         integer          :: mpi_phonons_id,mpi_phonons_nproc
         integer          :: mpi_phonons_world

        contains
               
 
         subroutine mpi_dist_nprocess(ntot,nloc,nstart,proc_grid,mpi_comm_loc)
         implicit none
         integer                :: i,k,j,rest,ntot,nloc,nstart
         integer                :: mpi_comm_loc,mpi_id_loc,mpi_nproc_loc
         integer, allocatable   :: proc_grid(:)


          call MPI_COMM_RANK(MPI_COMM_LOC,mpi_id_loc,err)
          call MPI_COMM_SIZE(MPI_COMM_LOC,mpi_nproc_loc,err)
 
        ! create groups of k point to be assigned to each mpi process
         
          nloc=ntot/mpi_nproc_loc
          rest=ntot-nloc*mpi_nproc_loc

          if(mpi_id_loc.lt.rest)then
           nloc=nloc+1
          endif


         ! create an array to map k to mpi_id

          if(allocated(proc_grid)) deallocate(proc_grid)
          allocate(proc_grid(ntot))

          do i=0,mpi_nproc_loc-1
           if(i.lt.rest)then
            k=1
           else
            k=0
           endif
           j=i*int(ntot/mpi_nproc_loc)+min(rest,i)+1
           k=int(ntot/mpi_nproc_loc)+j+k-1
           proc_grid(j:k)=i
          enddo

         ! index where to start for each process

          nstart=mpi_id_loc*int(ntot/mpi_nproc_loc)+min(rest,mpi_id_loc)+1

         return
         end subroutine mpi_dist_nprocess


        end module mpi_utils

