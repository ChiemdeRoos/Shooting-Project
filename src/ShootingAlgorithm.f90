module ShootingAlgo

use integration_module    

implicit none
        
private 
public ShootingMethod


contains


subroutine ShootingMethod(Matrix_V,EigVal,EigVec,N,h) ! change intents
    real*8, dimension(:,:), intent(in)  :: Matrix_V,EigVec
    real*8, dimension(:)                :: EigVal
    real*8, intent(in)                  :: h
    real*8, dimension(2)                :: dx 
    real*8                              :: inward_integral,outward_integral,EigVal2,error
    integer, intent(in)                 :: N
    integer                             :: xm,i,p,l
    real*8, allocatable, dimension(:)   :: Inward,Outward

    ! allocates space for the vectors that we are going to compare later on and that we fill with the outcome of EQ 10 & 11
    ! we use + 1 here to be able to calculate the derivative in xm
    ! we need another + 1 if the number of grid points is odd
    ! also use O which is the actual point xm.
    if (mod(N,2) == 0) then
        xm = N / 2
    else
        xm = N / 2 + 1
    endif
    allocate(Inward(xm+1))
    allocate(Outward(xm+1))
    Inward = 0.0
    Outward = 0.0

    ! Asks the user for the number of solutions and the desired error. Selects default if entered 0.
    print*,"The first how many solutions do you want?"
    read*, p
    print*,"Desired error? press 0 for default (default 10e-5)"
    read*, error
    
    if (abs(error) < 10e-9) then
        error = 10e-5
    endif
    

    ! loop to obtain all desired solutions
    do i = 1, p
        
        ! initializes the first 2 values of the inward and outward approximation using the first and last two EigVec values
        Inward(1:2) = EigVec(i,1:2)
        Outward(1) = EigVec(i,N)
        Outward(2) = EigVec(i,N-1)
        ! initializes the number of approximations back to the initial value of 1. 
        l = 1
        ! loop that keeps looping until the approximation is sufficient
        do 

            ! uses EQ 10 and 11 to fill the array for comparison
            call InOutFill(Inward,Outward,xm,h,Matrix_V,EigVal(i),N)
            ! calculates derivatives
            call Derivative(Inward,Outward,xm,h,dx)
            ! calculates integrals
            call Newton_cotes(Inward(1:xm)**2,h,1,xm,inward_integral)
            call Newton_cotes(Outward(1:xm)**2,h,1,xm,outward_integral)
            ! evaluates whether the approximation is succesful also calls on a function that calculates delta eigenvalue
            EigVal2 = EigVal(i) - DeltaEigVal(dx,Inward(xm),Outward(xm),inward_integral,outward_integral)!d_EigVal
            if (abs(EigVal(i) - EigVal2) < error) then 
                print*, "Approximation succesful!"
                print*, "Number of approximations required:", l
                exit
            else
                l = l + 1
                EigVal(i) = EigVal2
            endif

        enddo

        ! Creates a solution array and prints the solution.
        call NewSolution(N,xm,Outward,Inward,i)

    enddo

    ! deallocates allocated arrays
    deallocate(inward)
    deallocate(outward)
end subroutine ShootingMethod


! Calculates the derivative over 3 consecutive points
subroutine Derivative(Inward,Outward,xm,h,dx)
    real*8, dimension(:), intent(in)    :: Inward,Outward
    real*8, INTENT(IN)                  :: h
    real*8, dimension(2), intent(out)   :: dx
    integer, intent(in)                 :: xm

    dx(1) = (Inward(xm-1)-Inward(xm+1)) / (2 * h)
    dx(2) = (Outward(xm-1)-Outward(xm+1)) / (2 * h)
end subroutine Derivative 


! Normalizes vectors
subroutine normalize(input) ! change input to the array you use
    real*8, dimension(:), intent(inout) :: input
    real*8                              :: sum 
    integer                             :: i,n

    sum = 0.0
    n = size(input)
    do i = 1, n
        sum = sum + input(i)**2
    enddo

    input = input / sqrt(sum)
end subroutine normalize


subroutine NewSolution(N,xm,Outward,Inward,i)
    integer, intent(in)                 :: N,i,xm
    real*8, intent(in), dimension(:)    :: Inward,Outward
    real*8, allocatable, dimension(:)   :: Solution
    integer                             :: j

    ! allocates space for the solution array
    allocate(Solution(N))

    ! here we can normalize the solutions
    ! creation of the solution is dependent on whether the grid is even or uneven
    Solution(1:xm) = Inward(1:xm)
    if (mod(N,2) == 0) then
        do j = 1, xm
            Solution(xm+j) = Outward(size(Outward)-j)
        enddo
    else
        do j = 2, xm
            Solution(xm+j-1) = Outward(size(Outward)-j)
        enddo
    endif
    
    ! normalizes the solution, prints the solution and deallocates the array
    call normalize(Solution)
    print*, "Solution",i,":"
    print"(F20.8)", Solution
    deallocate(Solution)
end subroutine NewSolution


subroutine InOutFill(Inward,Outward,xm,h,Matrix_V,EigVal,N)
    real*8, dimension(:), intent(inout) :: Inward,Outward
    real*8, dimension(:,:), intent(in)  :: Matrix_V
    integer, intent(in)                 :: xm,N
    real*8, intent(in)                  :: EigVal,h
    integer                             :: i

    ! uses EQ 10 and 11 to fill the array for comparison and normilizes them.
    do i = 2, xm
        Inward(i+1) = -Inward(i-1) + 2 * h**2 *(Matrix_V(i,i) - EigVal + (1 / h**2)) * Inward(i)
        Outward(i+1) = -Outward(i-1) + 2 * h**2 *(Matrix_V(N-i,N-i) - EigVal + (1 / h**2)) * Outward(i)
    enddo
    call normalize(Inward)
    call normalize(Outward)
end subroutine InOutFill


! calculates delta eigenvalue
Function DeltaEigVal(dx,Inward,Outward,inward_integral,outward_integral) result(d_EigVal)
    real*8                              :: d_EigVal
    real*8, dimension(:), intent(in)    :: dx
    real*8, intent(in)                  :: Inward,Outward,inward_integral,outward_integral


    d_EigVal = (1.0 / 2.0) * ((dx(1) / Inward) - (dx(2) / Outward)) * &
    ((1/Outward**2) * outward_integral + &
    (1/Inward**2) * inward_integral)**(-1)
end function DeltaEigVal


end module ShootingAlgo