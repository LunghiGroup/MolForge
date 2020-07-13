        module parser_class
        implicit none

        contains

        
        subroutine get_line(l,cleaned_line,eof)
        implicit none
        integer                          :: l,size
        character(len = 1)               :: buffer
        character(len = :), allocatable  :: line,cleaned_line
        logical                          :: set,eof


30      continue

        set=.true.
        eof=.false.

        do 
        read(l,"(A)",ADVANCE='NO',SIZE=size,EOR=10,END=20) buffer

         if(buffer=='#')then

          if(set)then
            read(l,*)
            cycle
 
          else
 
           if(len_trim(line).eq.0)then
            read(l,*)
            set=.true.
            cycle        
           endif
           if(len_trim(line).ne.0)then
            read(l,*)
            goto 10
           endif
 
          endif
 
         endif

         if(buffer=='/')then
          read(l,*)
          cycle        
         endif

         if(set)then
          line=buffer(1:size)
          set=.false.
         else
          line = line // buffer
         endif
 
        enddo

10      continue      

        if (set) goto 30
        if (len_trim(line).eq.0) then
         read(l,*)
         goto 30
        endif


        eof=.false.      

        cleaned_line=trim(adjustl(line))

        return

20      continue


        write(*,*) 'End of input file'
        eof=.true.

        return
        end subroutine get_line


                
        subroutine get_word(line,cleaned_word,n)
        implicit none
        character(len = :), allocatable  :: line
        character(len = :), allocatable  :: word
        character(len = :), allocatable  :: cleaned_word
        character(len = 1)               :: buffer
        integer                          :: n,i,nword
        logical                          :: set,skip


         if (allocated(word)) deallocate(word)
         allocate(character(len=1)::word)
         word=''        

         skip=.false.

         nword=0

         do i=1,len_trim(line)

          buffer=line(i:i)

          if(buffer.ne.' ')then
                
            word=word//buffer
            skip=.true.

          else
           if(skip) nword=nword+1
           if(n.eq.nword)then
            goto 10
           else
             if (allocated(word)) deallocate(word)
             allocate(character(len=1)::word)
             word=''
             skip=.false.
           endif

          endif

         enddo

         if (len_trim(word).ne.0) nword=nword+1

         if(n.ne.nword)then
          write(*,*) 'Error parsing words, end of line reached'        
          stop
         endif

10       continue
         
         cleaned_word=trim(adjustl(word))

        return         
        end subroutine get_word



        subroutine to_upper(str)
        character(*), intent(in out) :: str
        integer :: i
 
         do i = 1, len(str)
          select case(str(i:i))
          case("a":"z")
           str(i:i) = achar(iachar(str(i:i))-32)
          end select
         end do 

        return
        end subroutine to_upper
 
 
 
        subroutine to_lower(str)
        character(*), intent(in out) :: str
        integer :: i
 
         do i = 1, len(str)
          select case(str(i:i))
           case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
          end select
         end do  

        return
        end subroutine to_lower




        end module parser_class
