!-----------------------------------------------------------------------
! MODULE: data_mod
!> @brief Global variables and dynamic memory allocation for the 2D
!>        Elastoplastic Boundary Element Method (ELQUABE + BEMECH 2D).
!-----------------------------------------------------------------------
module data_mod
    implicit none

    ! Problem dimensions and unit numbers
    integer :: N, NE, L, NN, NX, NX1
    integer :: INP, IPR, IP1, IP2, IP3, IP4

    ! Material properties (Shear modulus GE and Poisson's ratio XNU)
    real :: GE, XNU

    ! Dynamically allocated elastic BEM matrices, nodal vectors, and connectivity
    real, allocatable :: G(:,:), H(:,:)
    real, allocatable :: X(:), Y(:)
    integer, allocatable :: CONEC(:,:)
    real, allocatable :: FI(:), DFI(:)
    integer, allocatable :: KODE(:)
    integer, allocatable :: FIRST_G_DOF(:)
    real, allocatable :: CX(:), CY(:), SSOL(:), DSOL(:)

    ! Boundary nodal stresses (N, 4): [sigma_x, tau_xy, sigma_y, sigma_vm]
    ! and internal point von Mises stresses (L)
    real, allocatable :: SBOU(:,:), SVM_INT(:)

    !-------------------------------------------------------------------
    ! 2D Elastoplasticity parameters and arrays (from BEMECH 2D)
    !-------------------------------------------------------------------
    integer :: NTP, NCELL, NODC, NSIG, NTSF, NAUTO
    integer :: NCRIT, NINCS, MITER
    real    :: TOLGP, UNIAX, FRICT, HARDI, HARDK, HARDS, EM, FINCR, TOLER
    real    :: PR, PR1, PR2, PR3, PR4, PR5, CON, GAM, ROOT3, CLAME, CJUP
    real    :: ALFA, CEQS, CEQP

    ! Constant 2D indexing and multiplier arrays
    integer, parameter :: INX(4) = [1, 3, 2, 4]
    integer, parameter :: INY(4) = [1, 3, 3, 2]
    real,    parameter :: RE(4)  = [1.0, 1.0, 2.0, 1.0]
    real,    parameter :: DLT(2, 2) = reshape([1.0, 0.0, 0.0, 1.0], [2, 2])

    ! Unified nodal coordinates CD(2, NTP) and internal cell connectivity LNDC(NODC, NCELL)
    real,    allocatable :: CD(:,:)
    integer, allocatable :: LNDC(:,:)
    integer, allocatable :: NNOD(:)

    ! Elastoplastic stress, back-stress, plastic strain, and flow derivatives
    real, allocatable :: ESTRS(:,:), BACKS(:,:), DFDS(:,:), EPSTN(:)
    real, allocatable :: X_B(:)

    ! Influence matrices stored in RAM for 2D plasticity
    real, allocatable :: AB_MAT(:,:), AU_MAT(:,:), ASIG_MAT(:,:)
    real, allocatable :: EB_MAT(:,:), EU_MAT(:,:), ESIG_MAT(:,:)
    real, allocatable :: AC_MAT(:,:), EUI_MAT(:,:), E_MAT(:,:)

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
        NTP = N + L

        allocate(G(NX, NX1))
        allocate(H(NX, NX))
        allocate(X(N), Y(N))
        allocate(CONEC(NE, 3))
        allocate(FI(NX), DFI(NX1))
        allocate(KODE(NX1))
        allocate(FIRST_G_DOF(NN))
        allocate(SBOU(N, 4))
        allocate(CD(2, NTP))
        allocate(X_B(NN))

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

        G           = 0.0
        H           = 0.0
        X           = 0.0
        Y           = 0.0
        CONEC       = 0
        FI          = 0.0
        DFI         = 0.0
        KODE        = 0
        FIRST_G_DOF = 0
        SBOU        = 0.0
        CD          = 0.0
        X_B         = 0.0
        CX          = 0.0
        CY          = 0.0
        DSOL        = 0.0
        SSOL        = 0.0
        SVM_INT     = 0.0
    end subroutine allocate_arrays

    !-------------------------------------------------------------------
    !> @brief Allocates dynamic arrays for 2D elastoplasticity based on
    !>        NCELL, NODC, NSIG, and NTP.
    !-------------------------------------------------------------------
    subroutine allocate_plastic_arrays()
        implicit none

        NTSF = NSIG * NTP

        allocate(NNOD(NTP))
        allocate(EPSTN(NTP))
        allocate(ESTRS(NSIG, NTP))
        allocate(BACKS(NSIG, NTP))
        allocate(DFDS(NSIG, NTP))

        NNOD  = 0
        EPSTN = 0.0
        ESTRS = 0.0
        BACKS = 0.0
        DFDS  = 0.0

        if (NCELL > 0) then
            allocate(LNDC(NODC, NCELL))
            allocate(AB_MAT(NN, NN))
            allocate(ASIG_MAT(NTSF, NN))
            allocate(EB_MAT(NN, NTSF))
            allocate(ESIG_MAT(NTSF, NTSF))
            allocate(AC_MAT(NN, NTSF))
            allocate(E_MAT(NTSF, NTSF))

            LNDC     = 0
            AB_MAT   = 0.0
            ASIG_MAT = 0.0
            EB_MAT   = 0.0
            ESIG_MAT = 0.0
            AC_MAT   = 0.0
            E_MAT    = 0.0

            if (L > 0) then
                allocate(AU_MAT(2 * L, NN))
                allocate(EU_MAT(2 * L, NTSF))
                allocate(EUI_MAT(2 * L, NTSF))
                AU_MAT  = 0.0
                EU_MAT  = 0.0
                EUI_MAT = 0.0
            else
                allocate(AU_MAT(1, NN))
                allocate(EU_MAT(1, NTSF))
                allocate(EUI_MAT(1, NTSF))
                AU_MAT  = 0.0
                EU_MAT  = 0.0
                EUI_MAT = 0.0
            end if
        end if
    end subroutine allocate_plastic_arrays

    !-------------------------------------------------------------------
    !> @brief Deallocates all dynamic arrays.
    !-------------------------------------------------------------------
    subroutine deallocate_arrays()
        implicit none

        if (allocated(G))           deallocate(G)
        if (allocated(H))           deallocate(H)
        if (allocated(X))           deallocate(X)
        if (allocated(Y))           deallocate(Y)
        if (allocated(CONEC))       deallocate(CONEC)
        if (allocated(FI))          deallocate(FI)
        if (allocated(DFI))         deallocate(DFI)
        if (allocated(KODE))        deallocate(KODE)
        if (allocated(FIRST_G_DOF)) deallocate(FIRST_G_DOF)
        if (allocated(SBOU))        deallocate(SBOU)
        if (allocated(CD))          deallocate(CD)
        if (allocated(X_B))         deallocate(X_B)
        if (allocated(CX))          deallocate(CX)
        if (allocated(CY))          deallocate(CY)
        if (allocated(DSOL))        deallocate(DSOL)
        if (allocated(SSOL))        deallocate(SSOL)
        if (allocated(SVM_INT))     deallocate(SVM_INT)

        if (allocated(NNOD))        deallocate(NNOD)
        if (allocated(EPSTN))       deallocate(EPSTN)
        if (allocated(ESTRS))       deallocate(ESTRS)
        if (allocated(BACKS))       deallocate(BACKS)
        if (allocated(DFDS))        deallocate(DFDS)
        if (allocated(LNDC))        deallocate(LNDC)
        if (allocated(AB_MAT))      deallocate(AB_MAT)
        if (allocated(AU_MAT))      deallocate(AU_MAT)
        if (allocated(ASIG_MAT))    deallocate(ASIG_MAT)
        if (allocated(EB_MAT))      deallocate(EB_MAT)
        if (allocated(EU_MAT))      deallocate(EU_MAT)
        if (allocated(ESIG_MAT))    deallocate(ESIG_MAT)
        if (allocated(AC_MAT))      deallocate(AC_MAT)
        if (allocated(EUI_MAT))     deallocate(EUI_MAT)
        if (allocated(E_MAT))       deallocate(E_MAT)
    end subroutine deallocate_arrays

end module data_mod
