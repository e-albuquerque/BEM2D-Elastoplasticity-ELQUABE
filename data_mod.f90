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

    ! Dynamically allocated BEM matrices, nodal vectors, and connectivity
    real, allocatable :: G(:,:), H(:,:)
    real, allocatable :: X(:), Y(:)
    integer, allocatable :: CONEC(:,:)
    real, allocatable :: FI(:), DFI(:)
    integer, allocatable :: KODE(:)
    real, allocatable :: CX(:), CY(:), SSOL(:), DSOL(:)

    ! Boundary nodal stresses (N, 4): [sigma_x, tau_xy, sigma_y, sigma_vm]
    ! and internal point von Mises stresses (L)
    real, allocatable :: SBOU(:,:), SVM_INT(:)

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
        allocate(X(N), Y(N))
        allocate(CONEC(NE, 3))
        allocate(FI(NX), DFI(NX1))
        allocate(KODE(NX1))
        allocate(SBOU(N, 4))

        if (L > 0) then
            allocate(CX(L), CY(L))
            allocate(DSOL(2 * L))
            allocate(SSOL(3 * L))
            allocate(SVM_INT(L))
        else
            allocate(CX(1), CY(1))
            allocate(DSOL(1))
            allocate(SSOL(1))
            allocate(SVM_INT(1))
        end if

        G       = 0.0
        H       = 0.0
        X       = 0.0
        Y       = 0.0
        CONEC   = 0
        FI      = 0.0
        DFI     = 0.0
        KODE    = 0
        SBOU    = 0.0
        CX      = 0.0
        CY      = 0.0
        DSOL    = 0.0
        SSOL    = 0.0
        SVM_INT = 0.0
    end subroutine allocate_arrays

    !-------------------------------------------------------------------
    !> @brief Deallocates all dynamic arrays.
    !-------------------------------------------------------------------
    subroutine deallocate_arrays()
        implicit none

        if (allocated(G))       deallocate(G)
        if (allocated(H))       deallocate(H)
        if (allocated(X))       deallocate(X)
        if (allocated(Y))       deallocate(Y)
        if (allocated(CONEC))   deallocate(CONEC)
        if (allocated(FI))      deallocate(FI)
        if (allocated(DFI))     deallocate(DFI)
        if (allocated(KODE))    deallocate(KODE)
        if (allocated(SBOU))    deallocate(SBOU)
        if (allocated(CX))      deallocate(CX)
        if (allocated(CY))      deallocate(CY)
        if (allocated(DSOL))    deallocate(DSOL)
        if (allocated(SSOL))    deallocate(SSOL)
        if (allocated(SVM_INT)) deallocate(SVM_INT)
    end subroutine deallocate_arrays

end module data_mod
