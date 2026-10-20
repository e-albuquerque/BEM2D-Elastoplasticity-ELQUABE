!-----------------------------------------------------------------------
! MODULE: solver_mod
!> @brief System of linear equations solver module.
!-----------------------------------------------------------------------
module solver_mod
    implicit none

contains

    !-------------------------------------------------------------------
    !> @brief Solves the linear system A * X = B using Gaussian
    !>        elimination with partial row pivoting (SLNPD).
    !-------------------------------------------------------------------
    subroutine SLNPD(A, B, D, N, NX, IPR)
        implicit none
        integer, intent(in) :: N, NX, IPR
        real, intent(inout) :: A(NX, NX), B(NX)
        real, intent(out)   :: D

        real :: TOL, C
        integer :: I, J, K, K1, L, N1
        logical :: pivoted

        TOL = 1.0e-6
        N1  = N - 1

        do K = 1, N1
            K1 = K + 1
            C  = A(K, K)
            if (abs(C) <= TOL) then
                pivoted = .false.
                do J = K1, N
                    ! Try to interchange rows to get non-zero diagonal coefficient
                    if (abs(A(J, K)) > TOL) then
                        do L = K, N
                            C       = A(K, L)
                            A(K, L) = A(J, L)
                            A(J, L) = C
                        end do
                        C    = B(K)
                        B(K) = B(J)
                        B(J) = C
                        C    = A(K, K)
                        pivoted = .true.
                        exit
                    end if
                end do

                if (.not. pivoted) then
                    write(IPR, 2) K
                    D = 0.0
                    return
                end if
            end if

            ! Divide row by diagonal coefficient
            C = A(K, K)
            do J = K1, N
                A(K, J) = A(K, J) / C
            end do
            B(K) = B(K) / C

            ! Eliminate unknown X(K) from row I
            do I = K1, N
                C = A(I, K)
                do J = K1, N
                    A(I, J) = A(I, J) - C * A(K, J)
                end do
                B(I) = B(I) - C * B(K)
            end do
        end do

        ! Compute last unknown
        if (abs(A(N, N)) <= TOL) then
            write(IPR, 2) N
            D = 0.0
            return
        end if
        B(N) = B(N) / A(N, N)

        ! Apply backsubstitution process to compute remaining unknowns
        do L = 1, N1
            K  = N - L
            K1 = K + 1
            do J = K1, N
                B(K) = B(K) - A(K, J) * B(J)
            end do
        end do

        ! Compute value of determinant
        D = 1.0
        do I = 1, N
            D = D * A(I, I)
        end do

        return
2       format(' **** SINGULARITY IN ROW', I5)
    end subroutine SLNPD

end module solver_mod
