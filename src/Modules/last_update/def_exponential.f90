        module exponential_class
        implicit none


        contains


        subroutine exp_taylor(size_block,A,expA,step_min,max_iter)
        use mpi_utils
        use blacs_utils
        use units_parms
        implicit none
        type(dist_cmplx_mat)    :: A,expA 
        type(dist_cmplx_mat)    :: AA,BB
        integer                 :: size_block,i,ii,jj,j,l
        integer                 :: max_iter,indxl2g
        complex(8)              :: coeff,step_min
        logical, allocatable    :: check(:)


         if(allocated(expA%mat)) call expA%dealloc()
         call expA%set(size_block,size_block,NB,MB)
         allocate(check(mpi_blacs_nproc))

         call AA%set(size_block,size_block,NB,MB)
         call BB%set(size_block,size_block,NB,MB)

         do ii=1,size(expA%mat,1)
          do jj=1,size(expA%mat,2)
           i=indxl2g(ii,NB,myrow,0,nprow)
           j=indxl2g(jj,MB,mycol,0,npcol)
           if(i.ne.j)then
            expA%mat(ii,jj)=(0.0d0,0.0d0)
            AA%mat(ii,jj)=(0.0d0,0.0d0)
           else
            expA%mat(ii,jj)=(1.0d0,0.0d0)
            AA%mat(ii,jj)=(1.0d0,0.0d0)
           endif
          enddo
         enddo
                

        !  Taylor expansion for the starting Propagator               

         i=0
         check=.false.
         do while(.not. all(check) .and. i.lt.max_iter )
           
          i=i+1
          coeff=(1.0d0/dble(i))*(-1.0d0*step_min*cmplx(0.0d0,1.0d0,8))

          call pzgemm('N','N',size_block,size_block,size_block,&
                       coeff,AA%mat,1,1,AA%desc,A%mat,1,1,A%desc,&
                       (0.0d0,0.0d0),BB%mat,1,1,BB%desc) 

          do ii=1,size(expA%mat,1)
           do jj=1,size(expA%mat,2)
            expA%mat(ii,jj)=expA%mat(ii,jj)+BB%mat(ii,jj)
           enddo
          enddo

          i=i+1
          coeff=(1.0d0/dble(i))*(-1.0d0*step_min*cmplx(0.0d0,1.0d0,8))

          call pzgemm('N','N',size_block,size_block,size_block,&
                       coeff,BB%mat,1,1,BB%desc,A%mat,1,1,A%desc,&
                       (0.0d0,0.0d0),AA%mat,1,1,AA%desc) 
        

          do ii=1,size(expA%mat,1)
           do jj=1,size(expA%mat,2)
            expA%mat(ii,jj)=expA%mat(ii,jj)+AA%mat(ii,jj)
           enddo
          enddo

          if( maxval(abs(dble(AA%mat))).lt.1.0d-18 .and. &
              maxval(abs(aimag(AA%mat))).lt.1.0d-18 ) then
           check(mpi_blacs_id+1)=.true.
          else
           check(mpi_blacs_id+1)=.false.
          endif

          do l=0,mpi_blacs_nproc-1
           call mpi_bcast(check(l+1),1,mpi_logical,l,mpi_blacs_world,err)
          enddo

         enddo

         call AA%dealloc()
         call BB%dealloc()

         if(mpi_id.eq.0) &
         write(*,*) 'Taylor expansion converged in ',i,'steps'

        end subroutine exp_taylor


        subroutine mult_exp(size_block,expA,mult_fact)
        use mpi_utils
        use blacs_utils
        implicit none
        type(dist_cmplx_mat)    :: expA 
        type(dist_cmplx_mat)    :: AA,BB
        integer                 :: size_block,i
        integer                 :: mult_fact

         call AA%set(size_block,size_block,NB,MB)
         call BB%set(size_block,size_block,NB,MB)

         BB%mat=expA%mat

         do i=1,mult_fact-1

          call pzgemm('N','N',size_block,size_block,size_block,&
                       (1.0d0,0.0d0),BB%mat,1,1,BB%desc,expA%mat,1,1,expA%desc,&
                       (0.0d0,0.0d0),AA%mat,1,1,AA%desc) 

          BB%mat=AA%mat
            
         enddo

         expA%mat=BB%mat

         call AA%dealloc()
         call BB%dealloc()


        return
        end subroutine mult_exp


        end module exponential_class
