      module system_utils
        use iso_c_binding
        implicit none
      
      contains
      
        function get_pid() result(pid)
          integer :: pid
      
          interface
            function c_getpid() bind(C, name="getpid") result(cpid)
              import :: c_int
              integer(c_int) :: cpid
            end function
          end interface
      
          pid = c_getpid()
        end function get_pid


      
      end module system_utils
