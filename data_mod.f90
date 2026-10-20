!-----------------------------------------------------------------------
! MODULE: data_mod
!> @brief Global variables and dynamic memory allocation for ELQUABE.
!-----------------------------------------------------------------------
module data_mod
    implicit none

    ! Problem dimensions and unit numbers
    integer :: N, NE, L, NN, NX, NX1
    integer :: INP, IPR, IP1, IP2, IP3, IP4

    ! Material properties (Shear modulus GE and Poisson's ratio XNU)
    real :: GE, XNU

    ! Dynamically allocated BEM matrices and nodal vectors
    real, allocatable :: G(:,:), H(:,:)
    real, allocatable :: X(:), Y(:)
    real, allocatable :: FI(:), DFI(:)
    integer, allocatable :: KODE(:)
    real, allocatable :: CX(:), CY(:), SSOL(:), DSOL(:)

contains

    !-------------------------------------------------------------------
    !> @brief Allocates dynamic arrays based on NE (boundary elements)
    !>        and L (internal points).
    !-------------------------------------------------------------------
    subroutine allocate_arrays()
        implicit none

        N   = 2 * NE
        NN  = 2 * N
        NX  = 2 * N
        NX1 = 3 * N

        allocate(G(NX, NX1))
        allocate(H(NX, NX))
        allocate(X(N + 1), Y(N + 1))
        allocate(FI(NX), DFI(NX1))
        allocate(KODE(NX1))

        if (L > 0) then
            allocate(CX(L), CY(L))
            allocate(DSOL(2 * L))
            allocate(SSOL(3 * L))
        else
            allocate(CX(1), CY(1))
            allocate(DSOL(1))
            allocate(SSOL(1))
        end if

        G    = 0.0
        H    = 0.0
        X    = 0.0
        Y    = 0.0
        FI   = 0.0
        DFI  = 0.0
        KODE = 0
        CX   = 0.0
        CY   = 0.0
        DSOL = 0.0
        SSOL = 0.0
    end subroutine allocate_arrays

    !-------------------------------------------------------------------
    !> @brief Deallocates all dynamic arrays.
    !-------------------------------------------------------------------
    subroutine deallocate_arrays()
        implicit none

        if (allocated(G))    deallocate(G)
        if (allocated(H))    deallocate(H)
        if (allocated(X))    deallocate(X)
        if (allocated(Y))    deallocate(Y)
        if (allocated(FI))   deallocate(FI)
        if (allocated(DFI))  deallocate(DFI)
        if (allocated(KODE)) deallocate(KODE)
        if (allocated(CX))   deallocate(CX)
        if (allocated(CY))   deallocate(CY)
        if (allocated(DSOL)) deallocate(DSOL)
        if (allocated(SSOL)) deallocate(SSOL)
    end subroutine deallocate_arrays

end module data_mod
