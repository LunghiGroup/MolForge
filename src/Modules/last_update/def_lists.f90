        module lists_class
        implicit none
        


        type :: list_node
         class(*), pointer      :: key => null()
         class(list_node), pointer  :: next => null()
         class(list_node), pointer  :: prev => null()
        end type list_node 

        type list
         class(list_node), pointer :: head => null()
         class(list_node), pointer :: tail => null()
         class(list_node), pointer :: node => null()
         integer                   :: nelem         
         contains                 
         procedure                 :: init => init_list
         procedure                 :: skip => skip_node
         procedure                 :: rew  => rewind_node
         procedure                 :: rm => remove_node
         procedure                 :: reboot => reboot_list
         procedure                 :: delete => destroy_list
         procedure                 :: add_node
         procedure                 :: rd_dbl_node
         procedure                 :: rd_cmplx_node
         procedure                 :: rd_int_node
         generic                   :: rd_val => rd_int_node,rd_dbl_node,rd_cmplx_node
        end type list

        
        contains
      

        
!!!!!   lists general functions

        subroutine skip_node(this_list)
        implicit none      
        class(list)  :: this_list
         if( associated(this_list%node%next) ) then
          this_list%node=>this_list%node%next
         else
          this_list%node=>this_list%tail
         endif
        return
        end subroutine skip_node

        subroutine destroy_list(this_list)
        implicit none
        class(list) :: this_list
         do while (this_list%nelem.gt.0)
          this_list%node=>this_list%head
          call this_list%rm()
         enddo
        return
        end subroutine destroy_list
       
        subroutine remove_node(this_list)
        implicit none      
        class(list)  :: this_list
        class(list_node), pointer :: tmp_node

         if (this_list%nelem .eq. 0) return

         if (this_list%nelem .eq. 1)then

          if(associated(this_list%node%key)) then
             deallocate(this_list%node%key)
          endif
          if(associated(this_list%node)) then
             deallocate(this_list%node)
          endif
          this_list%node=>null()
          this_list%head=>null()
          this_list%tail=>null()

          this_list%nelem=this_list%nelem-1

         else

         if( .not. associated(this_list%node%next) .and. &
                   associated(this_list%node%prev) ) then
          this_list%tail=>this_list%node%prev
          tmp_node=>this_list%node
          this_list%node=>this_list%tail
          this_list%tail%next=>null()     
          if(associated(tmp_node%key)) then
            deallocate(tmp_node%key)
            tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
            deallocate(tmp_node)
            tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
         endif


         if( .not. associated(this_list%node%prev) .and. &
                   associated(this_list%node%next) ) then
          this_list%head=>this_list%node%next
          tmp_node=>this_list%node
          this_list%node=>this_list%head     
          this_list%head%prev=>null()     
          if(associated(tmp_node%key)) then
            deallocate(tmp_node%key)
            tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
            deallocate(tmp_node)
            tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
         endif

         if( associated(this_list%node%next) .and. &
             associated(this_list%node%prev) ) then
          this_list%node%prev%next=>this_list%node%next         
          this_list%node%next%prev=>this_list%node%prev         
          tmp_node=>this_list%node
          tmp_node%key=>this_list%node%key
          this_list%node=>this_list%node%next
          if(associated(tmp_node%key)) then
            deallocate(tmp_node%key)
            tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
            deallocate(tmp_node)
            tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
         endif

        endif

        return
        end subroutine remove_node


        subroutine rewind_node(this_list)
        implicit none      
        class(list)  :: this_list                
         if( associated(this_list%node%prev) ) then
          this_list%node=>this_list%node%prev
         else
          this_list%node=>this_list%head           
         endif 
        return
        end subroutine rewind_node


        subroutine reboot_list(this_list)
        implicit none      
        class(list)   :: this_list        
         if ( associated(this_list%node) &
              .and. associated(this_list%head )  ) then                 
          this_list%node=>this_list%head
         endif
        return
        end subroutine reboot_list


        subroutine last_node_list(this_list)
        implicit none      
        class(list)   :: this_list
         this_list%node=>this_list%tail
        return
        end subroutine last_node_list
      

        subroutine init_list(this_list)
        implicit none
        class(list)    :: this_list        
         this_list%nelem=0                        
        return
        end subroutine init_list


        subroutine add_node(this_list,val)
        implicit none
        class(list)                    :: this_list
        class(*),pointer               :: arrow
        class(*),optional              :: val
        class(list_node),pointer       :: tmp_node

         if(this_list%nelem.eq.0)then       
          allocate(this_list%head)
          allocate(this_list%node)
          if(present(val))then
           select type (val)
           type is (integer)
            allocate(integer::this_list%head%key)
            allocate(integer::this_list%node%key)
           type is (double precision)
            allocate(double precision::this_list%head%key)
            allocate(double precision::this_list%node%key)
           type is (complex(8))
            allocate(complex(8)::this_list%head%key)
            allocate(complex(8)::this_list%node%key)
           end select
          endif
          this_list%node=>this_list%head
          this_list%tail=>this_list%head
          this_list%nelem=this_list%nelem+1
          tmp_node=>this_list%head
         else
          allocate(tmp_node)
          if(present(val))then
           select type (val)
           type is (integer)
            allocate(integer::tmp_node%key)
           type is (double precision)
            allocate(double precision::tmp_node%key)
           type is (complex(8))
            allocate(complex(8)::tmp_node%key)
           end select
          endif
          this_list%tail%next=>tmp_node
          tmp_node%prev=>this_list%tail 
          this_list%tail=>tmp_node
          this_list%nelem=this_list%nelem+1
         endif


        
         if ( present(val) ) then
          select type (val)
                  
          type is (integer)
           select type (arrow=>tmp_node%key)
            type is (integer)
            arrow=val
           end select
          type is (double precision)
           select type (arrow=>tmp_node%key)
            type is (double precision)
            arrow=val
           end select
          type is (complex(8))
           select type (arrow=>tmp_node%key)
            type is (complex(8))
            arrow=val
           end select

          end select
         endif

         tmp_node=>null()

        return
        end subroutine add_node


        subroutine rd_int_node(this,val)
        implicit none
        class(list)           :: this
        integer               :: val
        class(*),pointer      :: bho

         select type (bho=>this%node%key)
          type is (integer)
           val=bho
         end select

        return
        end subroutine rd_int_node

        subroutine rd_dbl_node(this,val)
        implicit none
        class(list)           :: this
        double precision      :: val
        class(*),pointer      :: bho

         select type (bho=>this%node%key)
         type is (double precision)
           val=bho
         end select

        return
        end subroutine rd_dbl_node

        subroutine rd_cmplx_node(this,val)
        implicit none
        class(list)           :: this
        complex(8)            :: val
        class(*),pointer      :: bho

         select type (bho=>this%node%key)
          type is (complex(8))
           val=bho
         end select

        return
        end subroutine rd_cmplx_node


        end module lists_class




