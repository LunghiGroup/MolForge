        module matrix_class
        use mpi
        use mpi_utils
        use blacs_utils
        use sparse_class
        implicit none

         type :: gen_cmplx_matrix
          type(dist_cmplx_mat) :: di_mat
          type(csr_mat_cmplx)  :: sp_mat
         end type gen_cmplx_matrix

         type :: gen_dbl_matrix
          type(dist_dbl_mat) :: di_mat
          type(csr_mat_dbl)  :: sp_mat
         end type gen_dbl_matrix

        end module matrix_class

