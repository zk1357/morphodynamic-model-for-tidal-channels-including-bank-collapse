!*******************************************************************************
!   Copyright(C) 2004-2011 Intel Corporation. All Rights Reserved.
!   
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!   
!   *Third Party trademarks are the property of their respective owners.
!   
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!   Content : MKL PARDISO Fortran-90 example
!
!*******************************************************************************
!----------------------------------------------------------------------
! Example program to show the use of the "PARDISO" routine
! for symmetric linear systems
!---------------------------------------------------------------------
INCLUDE 'mkl_pardiso.f90'

SUBROUTINE pardiso_sym_f90( N , NNZ , A, JA , IA , B , X )
USE mkl_pardiso
IMPLICIT NONE
INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER, ALLOCATABLE :: ia( : )
INTEGER, ALLOCATABLE :: ja( : )
INTEGER, ALLOCATABLE :: perm( : )
REAL(KIND=DP), ALLOCATABLE :: a( : )
REAL(KIND=DP), ALLOCATABLE :: b( : )
REAL(KIND=DP), ALLOCATABLE :: x( : )
INTEGER i, idum(1)
REAL(KIND=DP) ddum(1)
!.. Fill all arrays containing matrix data.

nrhs = 1 
maxfct = 1 
mnum = 1
!ALLOCATE( ia ( n + 1 ) )
!ia = (/ 1, 5, 8, 10, 12, 15, 17, 18, 19 /)
!ALLOCATE( ja ( nnz ) )
!ja = (/ 1,    3,       6, 7,    &
!           2, 3,    5,          &
!              3,             8, &
!                 4,       7,    &
!                    5, 6, 7,    &
!                       6,    8, &
!                          7,    &
!                             8 /)
!ALLOCATE( a ( nnz ) )
!a = (/ 7.d0,        1.d0,             2.d0, 7.d0,        &
!             -4.d0, 8.d0,       2.d0,                    &
!                    1.d0,                         5.d0,  &
!                          7.d0,             9.d0,        &
!                                5.d0, 1.d0, 5.d0,        &
!                                     -1.d0,       5.d0,  &
!                                           11.d0,        &
!                                                  5.d0 /)
!ALLOCATE( b ( n ) )
!ALLOCATE( x ( n ) )
!..
!.. Set up PARDISO control parameter
!..
ALLOCATE( iparm ( 64 ) )
!ALLOCATE( PERM(N) )
do i = 1, 64
   iparm(i) = 0
end do 

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x
iparm(8) = 9 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturbe the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations

error  = 0 ! initialize error flag
msglvl = 0 ! no print statistical information
mtype  = -2 ! symmetric, indefinite

!!.. Initiliaze the internal solver memory pointer. This is only
!! necessary for the FIRST call of the PARDISO solver.
!
ALLOCATE ( pt ( 64 ) )
do i = 1, 64
   pt( i )%DUMMY =  0 
end do

!!.. Reordering and Symbolic Factorization, This step also allocates
!! all memory that is necessary for the factorization
!
!phase = 11 ! only reordering and symbolic factorization
!
!CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
!idum, nrhs, iparm, msglvl, ddum, ddum, error)
!    
!WRITE(*,*) 'Reordering completed ... '
!IF (error /= 0) THEN
!   WRITE(*,*) 'The following ERROR was detected: ', error
!   GOTO 1000
!END IF
!WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
!
!!.. Factorization.
!phase = 22 ! only factorization
!CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
!idum, nrhs, iparm, msglvl, ddum, ddum, error)
!WRITE(*,*) 'Factorization completed ... '
!IF (error /= 0) THEN
!   WRITE(*,*) 'The following ERROR was detected: ', error
!   GOTO 1000
!ENDIF
!
!!.. Back substitution and iterative refinement
!iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 13 ! only factorization
!do i = 1, n
!   b(i) = 1.d0
!end do
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
perm, nrhs, iparm, msglvl, b, x, error)
!WRITE(*,*) 'Solve completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   PAUSE
!   GOTO 1000
ENDIF
!WRITE(*,*) 'The solution of the system is '
!DO i = 1, n
!   WRITE(*,*) ' x(',i,') = ', x(i)
!END DO
! pause     
!1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF ( ALLOCATED( ia ) )      DEALLOCATE( ia )
IF ( ALLOCATED( ja ) )      DEALLOCATE( ja )
IF ( ALLOCATED( a ) )       DEALLOCATE( a )
!IF ( ALLOCATED( b ) )       DEALLOCATE( b )
!IF ( ALLOCATED( x ) )       DEALLOCATE( x )
IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )
                            deALLOCATE( pt )
!                            DEALLOCATE(PERM)
!IF (error1 /= 0) THEN
!   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
!   STOP 1
!ENDIF
!
!IF ( error /= 0 ) STOP 1
!!STOP 0

!write(*,*)j

!end do
END SUBROUTINE pardiso_sym_f90
