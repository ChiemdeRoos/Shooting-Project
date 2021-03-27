module Scheme

use Diagonalization

implicit none
    
private 
public ThreePointScheme
    
    
contains


Subroutine ThreePointScheme(N,h,Matrix_V,EigVal,EigVec,Grid)
    Integer, intent(in)                                 :: N
    Integer                                             :: problem
    real*8, intent(in)                                  :: h
    real*8, allocatable, dimension(:,:)                 :: Matrix_S,Matrix_L
    real*8, allocatable, dimension(:,:), intent(out)    :: EigVec,Matrix_V
    real*8, allocatable, dimension(:), intent(out)      :: EigVal
    real*8, dimension(:), intent(in)                    :: Grid

    ! calls the creation of the matrices for S and V, then calculates L from S and V and lastly diagonalizes L to obtain the eigenvalues and eigenvectors.
    call FillS(N,Matrix_S)
    ! print*, "Diagonlized matrix of S:"
    ! print"(10F7.2)", Matrix_S

    print*,"Which problem do you want to solve?"
    print*, "press 1: particle in a box"
    print*, "press 2: Gaussian potential"
    print*, "press 3: morse potential"
    read*,problem
    if (problem < 1 .OR. problem > 3) then 
        print*,"Invalid input, please try again."
        read*,problem
    endif

    ! selects what to solve based on user input. and fills the matrix with the potential.
    select case(problem)
        case(1)
            print*,"You selected Particle in a box."
            call ParticleBox(N,Matrix_V)
        case(2)
            print*,"You selected Gaussian potential."
            call GaussianPotential(N,Matrix_V,Grid)
        case(3)
            print*,"You selected Morse potential."
            call MorsePotential(N,Matrix_V,Grid)
    end select    

    ! fills the L and S matrix
    call Calculate_L(Matrix_L,Matrix_S,Matrix_V,N,h) 
    call Calculate_EigVal_EigVec(Matrix_L,EigVal,EigVec)
    !transposes the EigVec matrix
    EigVec = transpose(EigVec)

    deallocate(Matrix_L)
    deallocate(Matrix_S)
end subroutine ThreePointScheme


Subroutine FillS(N,Matrix_S)
    integer, intent(in)                                 :: N
    integer                                             :: i
    real*8, allocatable, dimension(:,:), intent(out)    :: Matrix_S

    ! Allocates space to matrix_S, initializes it to 0 and locates the number 1 for the first and last row
    allocate(Matrix_S(N,N))
    Matrix_S = 0.0
    Matrix_S(2,1) = 1
    Matrix_S(N - 1,N) = 1

    ! Places -2 on the diaganols
    do i = 1, N
        Matrix_S(i,i) = -2
    enddo

    ! Places all other 1 values
    do i = 2, N - 1
        Matrix_S(i - 1,i) = 1
        Matrix_S(i + 1,i) = 1
    enddo
end subroutine FillS


! If there is time i should expand this subroutine for other solutions than the particle in a box. If i do this i also need to integration grid in this subroutine
! Allocates the potential for a particle in a box
subroutine ParticleBox(N,Matrix_V)
    integer, intent(in)                                 :: N
    integer                                             :: i
    real*8, allocatable, dimension(:,:), intent(out)    :: Matrix_V

    ! allocates the matrix of the potential V to the grid size N, Then fills it with 0.0
    allocate(Matrix_V(N,N))
    Matrix_V = 0.0

    ! Fills the diagonals of the matrix of the potential V with 0 --> for a different example (not particle in a box) the value on the diagnoals may be different than 0.
    do i = 1,N
        Matrix_V(i,i) = 0.0
    enddo
end subroutine ParticleBox


subroutine GaussianPotential(N,Matrix_V,Grid)
    integer, intent(in)                                 :: N
    real*8, dimension(:), intent(in)                    :: Grid
    integer                                             :: i
    real*8, allocatable, dimension(:,:), intent(out)    :: Matrix_V

    ! allocates the matrix of the potential V to the grid size N, Then fills it with 0.0
    allocate(Matrix_V(N,N))
    Matrix_V = 0.0

    ! Fills the diagonals of the matrix of the potential V with 0 --> for a different example (not particle in a box) the value on the diagnoals may be different than 0.
    ! for this example i put V0 = 1 and alpha = 1.02
    do i = 1,N
        Matrix_V(i,i) = -3*exp(-0.1*Grid(i)**2)
    enddo
end subroutine GaussianPotential


subroutine MorsePotential(N,Matrix_V,Grid)
    integer, intent(in)                                 :: N
    real*8, dimension(:), intent(in)                    :: Grid
    integer                                             :: i
    real*8, allocatable, dimension(:,:), intent(out)    :: Matrix_V

    ! allocates the matrix of the potential V to the grid size N, Then fills it with 0.0
    allocate(Matrix_V(N,N))
    Matrix_V = 0.0

    ! Fills the diagonals of the matrix of the potential V with 0 --> for a different example (not particle in a box) the value on the diagnoals may be different than 0.
    ! for this example i used some data i could find for the H2 molecule. Not sure what the units are. D = 0.176, alpha = 1.02, xe = 0.80
    do i = 1,N
        Matrix_V(i,i) = 0.176 * (1 - exp(-1.02 * (Grid(i) - 0.80)**2))
    enddo
end subroutine MorsePotential


! Calculates the L value from formula 7
subroutine Calculate_L(Matrix_L,Matrix_S,Matrix_V,N,h) 
    real*8, allocatable, dimension(:,:), intent(out)    :: Matrix_L
    real*8, allocatable, dimension(:,:), intent(in)     :: Matrix_S, Matrix_V
    integer, intent(in)                                 :: N
    real*8, intent(in)                                  :: h

    ! allocates space for L
    allocate(Matrix_L(N,N))

    ! calculates all values in the matrix L
    Matrix_L = ((- 1) / (2 * h ** 2)) * Matrix_S + Matrix_V
end subroutine calculate_L


subroutine Calculate_EigVal_EigVec(Matrix_L,EigVal,EigVec)
    real*8, allocatable, dimension(:,:), intent(in)     :: Matrix_L
    real*8, allocatable, dimension(:,:), intent(out)    :: EigVec
    real*8, allocatable, dimension(:), intent(out)      :: EigVal

    ! allocates arrays for the eigenvectos and the eigenvalues
    allocate(EigVec(size(Matrix_L,1),size(Matrix_L,2)))
    allocate(EigVal(size(Matrix_L,1)))

    ! diagonlizes L to obtain eigenvectors and eigenvalues
    call diagonalize(Matrix_L,EigVec,EigVal)
end subroutine Calculate_EigVal_EigVec


end module Scheme