subroutine Sort_Mod_Unique_Int(n, val)


end

program unique_sort
    implicit none
    integer :: i = 0, min_val, max_val
    integer, dimension(10) :: val, unique
    integer, dimension(:), allocatable :: final

    val = [ 3,2,5,7,3,1,4,7,3,3 ]
    min_val = minval(val)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        i = i+1
        min_val = minval(val, mask=val>min_val)
        unique(i) = min_val
    enddo
    allocate(final(i), source=unique(1:i))   !<-- Or, just use unique(1:i) 
    print "(10i5:)", final
end program unique_sort
