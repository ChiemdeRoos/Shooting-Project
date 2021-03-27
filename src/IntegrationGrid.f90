module GridSetup

implicit none

private 
public IntegrationGrid


contains


Subroutine IntegrationGrid(N,h,Grid)
    integer, intent(out)                            :: N
    integer                                         :: i
    real*8                                          :: a,b
    real*8, intent(out)                             :: h
    real*8, dimension(:), allocatable, intent(out)  :: Grid
    
    ! Asks the user for the domain and the number of grid points
    print*, "Please enter the start of the domain."
    read*,a
    print*, "Please enter the end of the domain."
    read*,b
    print*, "Please enter the number of grid points."
    read*,N

    ! calculates interval
    h = (b - a) / (N - 1) ! I'm pretty sure you require n-1 here although this shouldnt matter much if you have many grid points.
    print*, "the interval between grid points is:"
    print"(F6.2)",h

    ! allocating the grid with user defined bounderies
    Allocate(grid(N)) 
    grid = 0.0

    ! filling the grid with x-values using the interval input.
    do i = 1, N
        grid(i) = (i - 1) * h
    enddo

    ! print*,"Grid:"
    ! print"(10F6.2)", grid
end Subroutine IntegrationGrid


end module GridSetup