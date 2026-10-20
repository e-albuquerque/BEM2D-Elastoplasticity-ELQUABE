!-----------------------------------------------------------------------
! MODULE: solver_mod
!> @brief System of linear equations solver module using LAPACK (DGESV).
!-----------------------------------------------------------------------
module solver_mod
    implicit none

    ! Explicit interface for LAPACK's double-precision real linear solver
    interface
        subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            implicit none
            integer, intent(in)    :: N, NRHS, LDA, LDB
            integer, intent(out)   :: IPIV(N), INFO
            real,    intent(inout) :: A(LDA, *), B(LDB, *)
        end subroutine DGESV
    end interface

contains

    !-------------------------------------------------------------------
    !> @brief Solves the linear system A * X = B using LAPACK's DGESV
    !>        (LU decomposition with partial pivoting and row interchanges).
    !-------------------------------------------------------------------
    subroutine SLNPD(A, B, D, N, NX, IPR)
        implicit none
        integer, intent(in) :: N, NX, IPR
        real, intent(inout) :: A(NX, NX), B(NX)
        real, intent(out)   :: D

        integer, allocatable :: IPIV(:)
        integer :: INFO, I

        allocate(IPIV(N))

        ! Solve A * X = B in-place via LAPACK DGESV
        call DGESV(N, 1, A, NX, IPIV, B, NX, INFO)

        if (INFO /= 0) then
            write(IPR, 2) INFO
            D = 0.0
            deallocate(IPIV)
            return
        end if

        ! Compute determinant from U diagonal factors and pivot parity
        D = 1.0
        do I = 1, N
            if (IPIV(I) /= I) D = -D
            D = D * A(I, I)
        end do

        deallocate(IPIV)
        return
2       format(' **** SINGULARITY OR LAPACK DGESV ERROR (INFO =', I5, ')')
    end subroutine SLNPD

end module solver_mod
