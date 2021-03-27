program Shooting

use GridSetup
use Scheme
use ShootingAlgo

implicit none

integer                             :: N
real*8                              :: h
real*8, allocatable, dimension(:,:) :: EigVec, Matrix_V
real*8, allocatable, dimension(:)   :: EigVal, Grid

! This program solves the schrodinger equation numerically using a shooting algorithm.
! This program can solve the following potentials: paritcle in a box, Gaussian, and Morse potential
! Since I was unable to find values that make sense for both the gaussian and morse potential I am currently unsure whether the implementation is succeful. 
! The particle in a box potential should work fine.

! Calls a subroutine that defines the integration grid size and intervals
Call IntegrationGrid(N,h,Grid)
! Executes the 3-point scheme and returns an array L with ?eigenvalues?
Call ThreePointScheme(N,h,Matrix_V,EigVal,EigVec,Grid) 
! Executes the shooting algorithm
Call ShootingMethod(Matrix_V,EigVal,EigVec,N,h)
   
end program Shooting