!-----------------------------------------------------------------------
! MODULE: assembly_mod
!> @brief 2D Elastostatic and Elastoplastic Boundary Element integration,
!>        matrix assembly, boundary/internal stress recovery, 2D internal
!>        cell integration (6-node triangular and 8-node quadrilateral
!>        cells), and incremental-iterative plasticity solver.
!-----------------------------------------------------------------------
module assembly_mod
    use data_mod
    use solver_mod
    implicit none

    integer, parameter :: NPG = 20
    real, parameter :: PI = 3.141592653589793

    ! 20-point Gauss-Legendre quadrature coordinates and weights (for EXTINEQ and SIGMAEQ)
    real, parameter :: GI(NPG) = [ &
        -0.9931285991850950, -0.9639719272779138, &
        -0.9122344282513259, -0.8391169718222188, &
        -0.7463319064601508, -0.6360536807265150, &
        -0.5108670019508271, -0.3737060887154195, &
        -0.2277858511416451, -0.0765265211334973, &
         0.0765265211334973,  0.2277858511416451, &
         0.3737060887154195,  0.5108670019508271, &
         0.6360536807265150,  0.7463319064601508, &
         0.8391169718222188,  0.9122344282513259, &
         0.9639719272779138,  0.9931285991850950 ]
    real, parameter :: OME(NPG) = [ &
         0.0176140071391509,  0.0406014298003864, &
         0.0626720483341088,  0.0832767415767047, &
         0.1019301198172407,  0.1181945319615186, &
         0.1316886384491769,  0.1420961093183824, &
         0.1491729864726042,  0.1527533871307263, &
         0.1527533871307263,  0.1491729864726042, &
         0.1420961093183824,  0.1316886384491769, &
         0.1181945319615186,  0.1019301198172407, &
         0.0832767415767047,  0.0626720483341088, &
         0.0406014298003864,  0.0176140071391509 ]

    ! 10-point standard Gauss quadrature (paired with 10-point logarithmic Gauss in LOCINEQ)
    real, parameter :: GI10(10) = [ &
        0.9739065285171717, -0.9739065285171717,  0.8650633666889845, -0.8650633666889845, &
        0.6794095682990244, -0.6794095682990244,  0.4333953941292472, -0.4333953941292472, &
        0.1488743389816312, -0.1488743389816312 ]
    real, parameter :: OME10(10) = [ &
        0.0666713443086881,  0.0666713443086881,  0.1494513491505806,  0.1494513491505806, &
        0.2190863625159820,  0.2190863625159820,  0.2692667193099963,  0.2692667193099963, &
        0.2955242247147529,  0.2955242247147529 ]

    ! 10-point logarithmic Gauss quadrature coordinates and weights
    real, parameter :: GIL(10) = [ &
        0.0090426309, 0.0539712662, 0.1353118246, 0.2470524162, &
        0.3802125396, 0.5237923179, 0.6657752055, 0.7941904160, &
        0.8981610912, 0.9688479887 ]
    real, parameter :: OMEL(10) = [ &
        0.1209551319, 0.1863635425, 0.1956608732, 0.1735771421, &
        0.1356956729, 0.0936467585, 0.0557877273, 0.0271598109, &
        0.0095151826, 0.0016381576 ]

    ! Intrinsic coordinates and side node definitions for 8-node quadrilateral cells
    real, parameter :: CORDL_QUAD(2, 8) = reshape([ &
        -1.0, -1.0,   1.0, -1.0,   1.0,  1.0,  -1.0,  1.0, &
         0.0, -1.0,   1.0,  0.0,   0.0,  1.0,  -1.0,  0.0 ], [2, 8])
    integer, parameter :: NODEF_QUAD(3, 4) = reshape([ &
        1, 2, 5,   2, 3, 6,   3, 4, 7,   4, 1, 8 ], [3, 4])

    ! Intrinsic area coordinates (L1, L2) and side node definitions for 6-node triangular cells
    real, parameter :: CORDL_TRI(2, 6) = reshape([ &
        1.0, 0.0,   0.0, 1.0,   0.0, 0.0, &
        0.5, 0.5,   0.0, 0.5,   0.5, 0.0 ], [2, 6])
    integer, parameter :: NODEF_TRI(3, 3) = reshape([ &
        1, 2, 4,   2, 3, 5,   3, 1, 6 ], [3, 3])

contains

    !-------------------------------------------------------------------
    !> @brief Computes the G and H matrices in parallel using OpenMP
    !>        and element connectivity CONEC, then forms the linear
    !>        system A * X = F (where A is stored in H and F in FI).
    !-------------------------------------------------------------------
    subroutine GHMATEQ()
        implicit none
        real :: HW(2, 6), GW(2, 6), CH
        real :: SF_G(NPG, 3), SF_G10(10, 3), SF_GL(10, 6)
        real :: XJA(NPG), ETA1(NPG), ETA2(NPG)
        real :: X1, Y1, X2, Y2, X3, Y3, A, B, C, D
        integer :: EL, I, J, K, LL, NODO
        integer :: NO1, NO2, NO3, NODE_LOC, DIR, NO_GLOB, COL_H, COL_G

        H  = 0.0
        G  = 0.0
        FI = 0.0

        ! Precompute shape functions at 20 standard Gauss points (for EXTINEQ)
        do I = 1, NPG
            SF_G(I, 1) = GI(I) * (GI(I) - 1.0) * 0.5
            SF_G(I, 2) = 1.0 - GI(I)**2
            SF_G(I, 3) = GI(I) * (GI(I) + 1.0) * 0.5
        end do

        ! Precompute shape functions at 10 Gauss and logarithmic Gauss points (for LOCINEQ)
        do I = 1, 10
            SF_G10(I, 1) = GI10(I) * (GI10(I) - 1.0) * 0.5
            SF_G10(I, 2) = 1.0 - GI10(I)**2
            SF_G10(I, 3) = GI10(I) * (GI10(I) + 1.0) * 0.5

            ! Columns 1..3: extreme node singular transformation (NODO = 1 or 3)
            SF_GL(I, 1) = (GIL(I) - 1.0) * (2.0 * GIL(I) - 1.0)
            SF_GL(I, 2) = 4.0 * GIL(I) * (1.0 - GIL(I))
            SF_GL(I, 3) = GIL(I) * (2.0 * GIL(I) - 1.0)
            ! Columns 4..6: mid-node singular transformation (NODO = 2)
            SF_GL(I, 4) = 0.5 * GIL(I) * (GIL(I) - 1.0)
            SF_GL(I, 5) = 1.0 - GIL(I)**2
            SF_GL(I, 6) = 0.5 * GIL(I) * (GIL(I) + 1.0)
        end do

        !$omp parallel do default(shared) &
        !$omp private(EL, LL, I, NO1, NO2, NO3, X1, Y1, X2, Y2, X3, Y3, A, B, C, D, &
        !$omp         XJA, ETA1, ETA2, NODO, HW, GW, K, NODE_LOC, NO_GLOB, DIR, J, COL_H, COL_G) &
        !$omp schedule(dynamic)
        do EL = 1, NE
            NO1 = CONEC(EL, 1)
            NO2 = CONEC(EL, 2)
            NO3 = CONEC(EL, 3)

            X1 = X(NO1); Y1 = Y(NO1)
            X2 = X(NO2); Y2 = Y(NO2)
            X3 = X(NO3); Y3 = Y(NO3)

            ! Compute Jacobians and outward unit normal vectors once for element EL
            A = X3 - 2.0 * X2 + X1
            B = (X3 - X1) / 2.0
            C = Y3 - 2.0 * Y2 + Y1
            D = (Y3 - Y1) / 2.0

            do I = 1, NPG
                XJA(I)  = sqrt((GI(I) * A + B)**2 + (GI(I) * C + D)**2)
                ETA1(I) =  (GI(I) * C + D) / XJA(I)
                ETA2(I) = -(GI(I) * A + B) / XJA(I)
            end do

            ! Loop over collocation (source) points
            do LL = 1, N
                if (LL == NO1) then
                    NODO = 1
                else if (LL == NO2) then
                    NODO = 2
                else if (LL == NO3) then
                    NODO = 3
                else
                    NODO = 0
                end if

                if (NODO > 0) then
                    call EXTINEQ(X(LL), Y(LL), X1, Y1, X2, Y2, X3, Y3, &
                                 SF_G, XJA, ETA1, ETA2, NODO, HW, GW)
                    call LOCINEQ(X1, Y1, X2, Y2, X3, Y3, SF_G10, SF_GL, GW, NODO)
                else
                    call EXTINEQ(X(LL), Y(LL), X1, Y1, X2, Y2, X3, Y3, &
                                 SF_G, XJA, ETA1, ETA2, NODO, HW, GW)
                end if

                do K = 1, 2
                    do NODE_LOC = 1, 3
                        NO_GLOB = CONEC(EL, NODE_LOC)
                        do DIR = 1, 2
                            J     = 2 * NODE_LOC - 2 + DIR
                            COL_G = 6 * EL - 6 + J
                            COL_H = 2 * NO_GLOB - 2 + DIR

                            G(2*LL-2+K, COL_G) = G(2*LL-2+K, COL_G) + GW(K, J)
                            !$omp atomic
                            H(2*LL-2+K, COL_H) = H(2*LL-2+K, COL_H) + HW(K, J)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        ! Compute the diagonal coefficients of the H matrix (rigid body motion)
        do I = 1, N
            H(2*I-1, 2*I-1) = 0.0
            H(2*I,   2*I-1) = 0.0
            H(2*I-1, 2*I)   = 0.0
            H(2*I,   2*I)   = 0.0
            do J = 1, N
                if (I /= J) then
                    H(2*I-1, 2*I-1) = H(2*I-1, 2*I-1) - H(2*I-1, 2*J-1)
                    H(2*I,   2*I-1) = H(2*I,   2*I-1) - H(2*I,   2*J-1)
                    H(2*I-1, 2*I)   = H(2*I-1, 2*I)   - H(2*I-1, 2*J)
                    H(2*I,   2*I)   = H(2*I,   2*I)   - H(2*I,   2*J)
                end if
            end do

            if (H(2*I-1, 2*I-1) < 0.0) then
                H(2*I-1, 2*I-1) = 1.0 + H(2*I-1, 2*I-1)
                H(2*I,   2*I)   = 1.0 + H(2*I,   2*I)
            end if
        end do

        ! Reorder the columns of the system of equations in accordance
        ! with the boundary conditions and connectivity CONEC to form
        ! system matrix A stored in H
        FIRST_G_DOF = 0
        do EL = 1, NE
            do NODE_LOC = 1, 3
                NO_GLOB = CONEC(EL, NODE_LOC)
                do DIR = 1, 2
                    J     = 2 * NODE_LOC - 2 + DIR
                    COL_G = 6 * EL - 6 + J
                    COL_H = 2 * NO_GLOB - 2 + DIR

                    if (KODE(COL_G) <= 0) then
                        if (FIRST_G_DOF(COL_H) == 0) then
                            FIRST_G_DOF(COL_H) = COL_G
                            do K = 1, NN
                                CH = H(K, COL_H)
                                H(K, COL_H) = -G(K, COL_G) * GE
                                G(K, COL_G) = -CH
                            end do
                        else
                            do K = 1, NN
                                H(K, COL_H) = H(K, COL_H) - G(K, COL_G) * GE
                                G(K, COL_G) = 0.0
                            end do
                        end if
                    end if
                end do
            end do
        end do

        ! Save system matrix AB_MAT for plasticity before H is factored by SLNPD
        if (NCELL > 0) then
            AB_MAT = H
        end if

        ! Form the right-hand side vector F stored in FI
        do I = 1, NN
            FI(I) = 0.0
            do J = 1, 6 * NE
                FI(I) = FI(I) + G(I, J) * DFI(J)
            end do
        end do

    end subroutine GHMATEQ

    !-------------------------------------------------------------------
    !> @brief Evaluates the hyperbolic sine (sinh) coordinate
    !>        transformation (Johnston & Elliott, 2005) to treat
    !>        near-singular integrals on boundary elements.
    !-------------------------------------------------------------------
    subroutine SINHTRANS(U, A, B_IN, XI, JT)
        implicit none
        real, intent(in)  :: U, A, B_IN
        real, intent(out) :: XI, JT
        real :: B, MI, ETA

        B = B_IN
        if (B < 1.0e-8) B = 1.0e-8

        MI  = 0.5 * (asinh((1.0 + A) / B) + asinh((1.0 - A) / B))
        ETA = 0.5 * (asinh((1.0 + A) / B) - asinh((1.0 - A) / B))
        XI  = A + B * sinh(MI * U - ETA)
        JT  = B * MI * cosh(MI * U - ETA)
    end subroutine SINHTRANS

    !-------------------------------------------------------------------
    !> @brief Computes the HW and GW matrices that relate a source
    !>        node (XP, YP) with a boundary element using NPG-point
    !>        Gauss quadrature.
    !-------------------------------------------------------------------
    subroutine EXTINEQ(XP, YP, X1, Y1, X2, Y2, X3, Y3, SF_G, XJA, ETA1, ETA2, NODO, HW, GW)
        implicit none
        real, intent(in)    :: XP, YP, X1, Y1, X2, Y2, X3, Y3
        real, intent(in)    :: SF_G(NPG, 3), XJA(NPG), ETA1(NPG), ETA2(NPG)
        integer, intent(in) :: NODO
        real, intent(out)   :: HW(2, 6), GW(2, 6)

        real :: A, B, C, D, DE
        real :: DX_EL, DY_EL, NORM_DELEM2, NORM_DELEM
        real :: EET, N1G, N2G, N3G, PSX, PSY, B_DIST, XI, JT, W_EFF
        real :: F1, F2, F3, XCO, YCO, XJ, E1, E2, RA, RD1, RD2, RDN
        integer :: I

        HW = 0.0
        GW = 0.0
        DE = 4.0 * PI * PR1

        if (NODO == 0) then
            A = X3 - 2.0 * X2 + X1
            B = (X3 - X1) / 2.0
            C = Y3 - 2.0 * Y2 + Y1
            D = (Y3 - Y1) / 2.0

            DX_EL       = X3 - X1
            DY_EL       = Y3 - Y1
            NORM_DELEM2 = DX_EL**2 + DY_EL**2
            NORM_DELEM  = sqrt(NORM_DELEM2)

            EET = 2.0 * (DX_EL * (XP - X1) + DY_EL * (YP - Y1)) / NORM_DELEM2 - 1.0

            N1G = 0.5 * EET * (EET - 1.0)
            N2G = 1.0 - EET**2
            N3G = 0.5 * EET * (1.0 + EET)

            PSX    = N1G * X1 + N2G * X2 + N3G * X3
            PSY    = N1G * Y1 + N2G * Y2 + N3G * Y3
            B_DIST = 2.0 * sqrt((PSX - XP)**2 + (PSY - YP)**2) / NORM_DELEM
        end if

        do I = 1, NPG
            if (NODO == 0) then
                call SINHTRANS(GI(I), EET, B_DIST, XI, JT)
                F1    = 0.5 * XI * (XI - 1.0)
                F2    = 1.0 - XI**2
                F3    = 0.5 * XI * (XI + 1.0)
                XJ    = sqrt((XI * A + B)**2 + (XI * C + D)**2)
                E1    =  (XI * C + D) / XJ
                E2    = -(XI * A + B) / XJ
                W_EFF = OME(I) * XJ * JT
            else
                F1    = SF_G(I, 1)
                F2    = SF_G(I, 2)
                F3    = SF_G(I, 3)
                XJ    = XJA(I)
                E1    = ETA1(I)
                E2    = ETA2(I)
                W_EFF = OME(I) * XJ
            end if

            XCO = X1 * F1 + X2 * F2 + X3 * F3
            YCO = Y1 * F1 + Y2 * F2 + Y3 * F3
            RA  = sqrt((XP - XCO)**2 + (YP - YCO)**2)
            RD1 = (XCO - XP) / RA
            RD2 = (YCO - YP) / RA
            RDN = RD1 * E1 + RD2 * E2

            GW(1, 1) = GW(1, 1) + (PR5 * log(1.0 / RA) + RD1**2) * W_EFF * F1 / (2.0 * DE * GE)
            GW(1, 2) = GW(1, 2) + RD1 * RD2 * W_EFF * F1 / (2.0 * DE * GE)
            GW(2, 1) = GW(1, 2)
            GW(2, 2) = GW(2, 2) + (PR5 * log(1.0 / RA) + RD2**2) * W_EFF * F1 / (2.0 * DE * GE)
            HW(1, 1) = HW(1, 1) - RDN * (PR2 + 2.0 * RD1**2) / (RA * DE) * W_EFF * F1
            HW(1, 2) = HW(1, 2) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E1 * RD2 - E2 * RD1)) * &
                                  W_EFF * F1 / (RA * DE)
            HW(2, 1) = HW(2, 1) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E2 * RD1 - E1 * RD2)) * &
                                  W_EFF * F1 / (RA * DE)
            HW(2, 2) = HW(2, 2) - RDN * (PR2 + 2.0 * RD2**2) * W_EFF * F1 / (RA * DE)

            GW(1, 3) = GW(1, 3) + (PR5 * log(1.0 / RA) + RD1**2) * W_EFF * F2 / (2.0 * DE * GE)
            GW(1, 4) = GW(1, 4) + RD1 * RD2 * W_EFF * F2 / (2.0 * DE * GE)
            GW(2, 3) = GW(1, 4)
            GW(2, 4) = GW(2, 4) + (PR5 * log(1.0 / RA) + RD2**2) * W_EFF * F2 / (2.0 * DE * GE)
            HW(1, 3) = HW(1, 3) - RDN * (PR2 + 2.0 * RD1**2) / (RA * DE) * W_EFF * F2
            HW(1, 4) = HW(1, 4) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E1 * RD2 - E2 * RD1)) * &
                                  W_EFF * F2 / (RA * DE)
            HW(2, 3) = HW(2, 3) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E2 * RD1 - E1 * RD2)) * &
                                  W_EFF * F2 / (RA * DE)
            HW(2, 4) = HW(2, 4) - RDN * (PR2 + 2.0 * RD2**2) * W_EFF * F2 / (RA * DE)

            GW(1, 5) = GW(1, 5) + (PR5 * log(1.0 / RA) + RD1**2) * W_EFF * F3 / (2.0 * DE * GE)
            GW(1, 6) = GW(1, 6) + RD1 * RD2 * W_EFF * F3 / (2.0 * DE * GE)
            GW(2, 5) = GW(1, 6)
            GW(2, 6) = GW(2, 6) + (PR5 * log(1.0 / RA) + RD2**2) * W_EFF * F3 / (2.0 * DE * GE)
            HW(1, 5) = HW(1, 5) - RDN * (PR2 + 2.0 * RD1**2) / (RA * DE) * W_EFF * F3
            HW(1, 6) = HW(1, 6) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E1 * RD2 - E2 * RD1)) * &
                                  W_EFF * F3 / (RA * DE)
            HW(2, 5) = HW(2, 5) - (RDN * 2.0 * RD1 * RD2 + PR2 * (E2 * RD1 - E1 * RD2)) * &
                                  W_EFF * F3 / (RA * DE)
            HW(2, 6) = HW(2, 6) - RDN * (PR2 + 2.0 * RD2**2) * W_EFF * F3 / (RA * DE)
        end do

    end subroutine EXTINEQ

    !-------------------------------------------------------------------
    !> @brief Computes the GW matrix when the collocation point is one
    !>        of the nodes of the integration element (singular/logarithmic
    !>        quadrature).
    !-------------------------------------------------------------------
    subroutine LOCINEQ(XG1, YG1, XG2, YG2, XG3, YG3, SF_G, SF_GL, GW, NODO)
        implicit none
        real, intent(in)    :: XG1, YG1, XG2, YG2, XG3, YG3
        real, intent(in)    :: SF_G(10, 3), SF_GL(10, 6)
        integer, intent(in) :: NODO
        real, intent(out)   :: GW(2, 6)

        real :: X1, Y1, X2, Y2, X3, Y3
        real :: A1, B1, A2, B2, A, B, C, CONT1, CONT2
        real :: T1, T2, T3, F1, F2, F3, FL1, FL2, FL3, FLN1, FLN2, FLN3
        real :: XJA1, XJA11, XJA2, XLO, S1, S2, S3
        integer :: I

        X1 = 0.0; Y1 = 0.0; X2 = 0.0; Y2 = 0.0; X3 = 0.0; Y3 = 0.0
        A1 = 0.0; B1 = 0.0; A2 = 0.0; B2 = 0.0
        XJA1 = 0.0; XJA11 = 0.0; XJA2 = 0.0; XLO = 0.0
        S1 = 0.0; S2 = 0.0; S3 = 0.0

        select case (NODO)
        case (1)
            X3 = XG3 - XG1; Y3 = YG3 - YG1
            X2 = XG2 - XG1; Y2 = YG2 - YG1
            A1 = (X3 - 2.0 * X2) * 0.5; B1 = X2
            A2 = (Y3 - 2.0 * Y2) * 0.5; B2 = Y2
        case (2)
            X3 = XG3 - XG2; Y3 = YG3 - YG2
            X1 = XG1 - XG2; Y1 = YG1 - YG2
            A1 = X1 + X3; B1 = X3 - X1
            A2 = Y1 + Y3; B2 = Y3 - Y1
        case (3)
            X2 = XG2 - XG3; Y2 = YG2 - YG3
            X1 = XG1 - XG3; Y1 = YG1 - YG3
            A1 = (X1 - 2.0 * X2) * 0.5; B1 = -X2
            A2 = (Y1 - 2.0 * Y2) * 0.5; B2 = -Y2
        end select

        GW    = 0.0
        A     = A1**2 + A2**2
        B     = 2.0 * (A1 * B1 + A2 * B2)
        C     = B1**2 + B2**2
        CONT1 = PR5 / (8.0 * PI * GE * PR1)
        CONT2 = CONT1 / PR5

        do I = 1, 10
            T1 = ((A1 * GI10(I) + B1)**2) / ((A2 * GI10(I) + B2)**2 + (A1 * GI10(I) + B1)**2)
            T2 = ((A2 * GI10(I) + B2)**2) / ((A2 * GI10(I) + B2)**2 + (A1 * GI10(I) + B1)**2)
            T3 = (A1 * GI10(I) + B1) * (A2 * GI10(I) + B2) / ((A1 * GI10(I) + B1)**2 + (A2 * GI10(I) + B2)**2)

            F1   = SF_G(I, 1); F2   = SF_G(I, 2); F3   = SF_G(I, 3)
            FL1  = SF_GL(I, 1); FL2  = SF_GL(I, 2); FL3  = SF_GL(I, 3)
            FLN1 = SF_GL(I, 4); FLN2 = SF_GL(I, 5); FLN3 = SF_GL(I, 6)

            select case (NODO)
            case (1)
                XJA1 = sqrt((4.0 * A1 * GIL(I) - 2.0 * A1 + 0.5 * X3)**2 + &
                            (4.0 * A2 * GIL(I) - 2.0 * A2 + 0.5 * Y3)**2) * 2.0
                XJA2 = sqrt((A1 * GI10(I) * 2.0 + 0.5 * X3)**2 + (A2 * GI10(I) * 2.0 + 0.5 * Y3)**2)
                XLO  = -log(2.0 * sqrt((GI10(I) * A1 + B1)**2 + (GI10(I) * A2 + B2)**2))
                S3   = CONT1 * (FL3 * XJA1 * OMEL(I) + F3 * XJA2 * XLO * OME10(I))
                S2   = CONT1 * (FL2 * XJA1 * OMEL(I) + F2 * XJA2 * XLO * OME10(I))
                S1   = CONT1 * (FL1 * XJA1 * OMEL(I) + F1 * XJA2 * XLO * OME10(I))
            case (2)
                XJA1  = sqrt((0.5 * B1 - A1 * GIL(I))**2 + (0.5 * B2 - A2 * GIL(I))**2)
                XJA11 = sqrt((0.5 * B1 + A1 * GIL(I))**2 + (0.5 * B2 + A2 * GIL(I))**2)
                XJA2  = sqrt((0.5 * B1 + A1 * GI10(I))**2 + (0.5 * B2 + A2 * GI10(I))**2)
                XLO   = -0.5 * log((GI10(I) * A1 * 0.5 + B1 * 0.5)**2 + (GI10(I) * A2 * 0.5 + B2 * 0.5)**2)
                S3    = CONT1 * ((FLN1 * XJA1 + FLN3 * XJA11) * OMEL(I) + F3 * XJA2 * XLO * OME10(I))
                S2    = CONT1 * (FLN2 * (XJA1 + XJA11) * OMEL(I) + F2 * XJA2 * XLO * OME10(I))
                S1    = CONT1 * ((FLN3 * XJA1 + FLN1 * XJA11) * OMEL(I) + F1 * XJA2 * XLO * OME10(I))
            case (3)
                XJA1 = sqrt((2.0 * A1 - 4.0 * A1 * GIL(I) - 0.5 * X1)**2 + &
                            (2.0 * A2 - 4.0 * A2 * GIL(I) - 0.5 * Y1)**2) * 2.0
                XJA2 = sqrt((2.0 * A1 * GI10(I) - 0.5 * X1)**2 + (2.0 * A2 * GI10(I) - 0.5 * Y1)**2)
                XLO  = -log(2.0 * sqrt((A1 * GI10(I) + B1)**2 + (A2 * GI10(I) + B2)**2))
                S3   = CONT1 * (FL1 * XJA1 * OMEL(I) + F3 * XJA2 * XLO * OME10(I))
                S2   = CONT1 * (FL2 * XJA1 * OMEL(I) + F2 * XJA2 * XLO * OME10(I))
                S1   = CONT1 * (FL3 * XJA1 * OMEL(I) + F1 * XJA2 * XLO * OME10(I))
            end select

            GW(1, 5) = GW(1, 5) + S3 + CONT2 * F3 * T1 * XJA2 * OME10(I)
            GW(1, 3) = GW(1, 3) + S2 + CONT2 * F2 * T1 * XJA2 * OME10(I)
            GW(1, 1) = GW(1, 1) + S1 + CONT2 * F1 * T1 * XJA2 * OME10(I)
            GW(1, 6) = GW(1, 6) + CONT2 * T3 * XJA2 * F3 * OME10(I)
            GW(1, 4) = GW(1, 4) + CONT2 * T3 * XJA2 * F2 * OME10(I)
            GW(1, 2) = GW(1, 2) + CONT2 * T3 * XJA2 * F1 * OME10(I)
            GW(2, 5) = GW(1, 6)
            GW(2, 3) = GW(1, 4)
            GW(2, 1) = GW(1, 2)
            GW(2, 6) = GW(2, 6) + S3 + CONT2 * T2 * XJA2 * F3 * OME10(I)
            GW(2, 4) = GW(2, 4) + S2 + CONT2 * T2 * XJA2 * F2 * OME10(I)
            GW(2, 2) = GW(2, 2) + S1 + CONT2 * T2 * XJA2 * F1 * OME10(I)
        end do

    end subroutine LOCINEQ

    !-------------------------------------------------------------------
    !> @brief Rearranges FI and DFI arrays using CONEC to store boundary
    !>        displacements and tractions, computes boundary stresses via
    !>        shape function derivatives, evaluates displacements and
    !>        stresses at internal points, and builds AU_MAT and ASIG_MAT
    !>        when NCELL > 0.
    !-------------------------------------------------------------------
    subroutine INTEREQ()
        implicit none
        real :: HW(2, 6), GW(2, 6)
        real :: D11(6), D12(6), D22(6), S11(6), S12(6), S22(6)
        real :: SF_G(NPG, 3)
        real :: XJA_EL(NE, NPG), ETA1_EL(NE, NPG), ETA2_EL(NE, NPG)
        real :: A, B, C, D
        real :: CH, SX, TXY, SY, SZ
        integer :: EL, I, J, K, ISG0
        integer :: NO1, NO2, NO3, NODE_LOC, DIR, NO_GLOB, COL_H, COL_G
        integer :: FIRST_DOF_LOC(NN)

        ! Save boundary unknown solution vector X_B before rearranging FI and DFI
        X_B = FI

        ! Rearrange FI and DFI arrays using element connectivity CONEC
        ! so FI stores all boundary displacements and DFI stores all tractions
        FIRST_DOF_LOC = 0
        do EL = 1, NE
            do NODE_LOC = 1, 3
                NO_GLOB = CONEC(EL, NODE_LOC)
                do DIR = 1, 2
                    J     = 2 * NODE_LOC - 2 + DIR
                    COL_G = 6 * EL - 6 + J
                    COL_H = 2 * NO_GLOB - 2 + DIR

                    if (KODE(COL_G) <= 0) then
                        if (FIRST_DOF_LOC(COL_H) == 0) then
                            FIRST_DOF_LOC(COL_H) = COL_G
                            CH         = FI(COL_H) * GE
                            FI(COL_H)  = DFI(COL_G)
                            DFI(COL_G) = CH
                        else
                            DFI(COL_G) = DFI(FIRST_DOF_LOC(COL_H))
                        end if
                    end if
                end do
            end do
        end do

        ! Compute boundary nodal stresses and boundary rows of ASIG_MAT & ESIG_MAT
        call BOUND_STRESS()

        ! Store normalized elastic boundary stresses in ESTRS(:, 1:N)
        do I = 1, N
            ESTRS(1, I) = SBOU(I, 1) / (2.0 * GE)
            ESTRS(2, I) = SBOU(I, 3) / (2.0 * GE)
            ESTRS(3, I) = SBOU(I, 2) / (2.0 * GE)
            if (NSIG == 4) then
                ESTRS(4, I) = PR * (ESTRS(1, I) + ESTRS(2, I))
            end if
        end do

        ! Compute the values of stresses and displacements at internal points
        if (L == 0) return

        do I = 1, NPG
            SF_G(I, 1) = GI(I) * (GI(I) - 1.0) * 0.5
            SF_G(I, 2) = 1.0 - GI(I)**2
            SF_G(I, 3) = GI(I) * (GI(I) + 1.0) * 0.5
        end do

        do EL = 1, NE
            NO1 = CONEC(EL, 1)
            NO2 = CONEC(EL, 2)
            NO3 = CONEC(EL, 3)
            A = X(NO3) - 2.0 * X(NO2) + X(NO1)
            B = (X(NO3) - X(NO1)) / 2.0
            C = Y(NO3) - 2.0 * Y(NO2) + Y(NO1)
            D = (Y(NO3) - Y(NO1)) / 2.0
            do I = 1, NPG
                XJA_EL(EL, I)  = sqrt((GI(I) * A + B)**2 + (GI(I) * C + D)**2)
                ETA1_EL(EL, I) =  (GI(I) * C + D) / XJA_EL(EL, I)
                ETA2_EL(EL, I) = -(GI(I) * A + B) / XJA_EL(EL, I)
            end do
        end do

        !$omp parallel do default(shared) &
        !$omp private(K, EL, NO1, NO2, NO3, NODE_LOC, NO_GLOB, DIR, J, COL_G, COL_H, &
        !$omp         HW, GW, D11, D12, D22, S11, S12, S22, SX, TXY, SY, SZ, ISG0) &
        !$omp schedule(dynamic)
        do K = 1, L
            SSOL(3*K-2) = 0.0
            SSOL(3*K-1) = 0.0
            SSOL(3*K)   = 0.0
            DSOL(2*K-1) = 0.0
            DSOL(2*K)   = 0.0
            ISG0        = NSIG * (N + K - 1)

            do EL = 1, NE
                NO1 = CONEC(EL, 1)
                NO2 = CONEC(EL, 2)
                NO3 = CONEC(EL, 3)

                call EXTINEQ(CX(K), CY(K), X(NO1), Y(NO1), X(NO2), Y(NO2), &
                             X(NO3), Y(NO3), SF_G, XJA_EL(EL, :), ETA1_EL(EL, :), &
                             ETA2_EL(EL, :), 0, HW, GW)
                call SIGMAEQ(CX(K), CY(K), X(NO1), Y(NO1), X(NO2), Y(NO2), &
                             X(NO3), Y(NO3), D11, D12, D22, S11, S12, S22)

                do NODE_LOC = 1, 3
                    NO_GLOB = CONEC(EL, NODE_LOC)
                    do DIR = 1, 2
                        J     = 2 * NODE_LOC - 2 + DIR
                        COL_G = 6 * EL - 6 + J
                        COL_H = 2 * NO_GLOB - 2 + DIR

                        SSOL(3*K-2) = SSOL(3*K-2) + D11(J) * DFI(COL_G) - S11(J) * FI(COL_H)
                        SSOL(3*K-1) = SSOL(3*K-1) + D12(J) * DFI(COL_G) - S12(J) * FI(COL_H)
                        SSOL(3*K)   = SSOL(3*K)   + D22(J) * DFI(COL_G) - S22(J) * FI(COL_H)
                        DSOL(2*K-1) = DSOL(2*K-1) + GW(1, J) * DFI(COL_G) - HW(1, J) * FI(COL_H)
                        DSOL(2*K)   = DSOL(2*K)   + GW(2, J) * DFI(COL_G) - HW(2, J) * FI(COL_H)

                        ! Assemble AU_MAT and ASIG_MAT for plasticity (INNERPS)
                        if (NCELL > 0) then
                            if (FIRST_G_DOF(COL_H) == 0) then
                                ! COL_H unknown is displacement u(COL_H)
                                AU_MAT(2*K-1, COL_H)     = AU_MAT(2*K-1, COL_H)     - HW(1, J)
                                AU_MAT(2*K,   COL_H)     = AU_MAT(2*K,   COL_H)     - HW(2, J)
                                ASIG_MAT(ISG0+1, COL_H)  = ASIG_MAT(ISG0+1, COL_H)  - S11(J) / (2.0 * GE)
                                ASIG_MAT(ISG0+2, COL_H)  = ASIG_MAT(ISG0+2, COL_H)  - S22(J) / (2.0 * GE)
                                ASIG_MAT(ISG0+3, COL_H)  = ASIG_MAT(ISG0+3, COL_H)  - S12(J) / (2.0 * GE)
                            else if (KODE(COL_G) <= 0) then
                                ! COL_H unknown is traction t(COL_G) / GE
                                AU_MAT(2*K-1, COL_H)     = AU_MAT(2*K-1, COL_H)     + GW(1, J) * GE
                                AU_MAT(2*K,   COL_H)     = AU_MAT(2*K,   COL_H)     + GW(2, J) * GE
                                ASIG_MAT(ISG0+1, COL_H)  = ASIG_MAT(ISG0+1, COL_H)  + 0.5 * D11(J)
                                ASIG_MAT(ISG0+2, COL_H)  = ASIG_MAT(ISG0+2, COL_H)  + 0.5 * D22(J)
                                ASIG_MAT(ISG0+3, COL_H)  = ASIG_MAT(ISG0+3, COL_H)  + 0.5 * D12(J)
                            end if
                        end if
                    end do
                end do
            end do

            if (NCELL > 0 .and. NSIG == 4) then
                ASIG_MAT(ISG0+4, 1:NN) = PR * (ASIG_MAT(ISG0+1, 1:NN) + ASIG_MAT(ISG0+2, 1:NN))
            end if

            SX  = SSOL(3*K-2)
            TXY = SSOL(3*K-1)
            SY  = SSOL(3*K)
            ESTRS(1, N + K) = SX  / (2.0 * GE)
            ESTRS(2, N + K) = SY  / (2.0 * GE)
            ESTRS(3, N + K) = TXY / (2.0 * GE)
            if (NSIG == 4) then
                SZ = PR * (SX + SY)
                ESTRS(4, N + K) = SZ / (2.0 * GE)
                SVM_INT(K) = sqrt(0.5 * ((SX - SY)**2 + (SY - SZ)**2 + (SZ - SX)**2) + 3.0 * TXY**2)
            else
                SVM_INT(K) = sqrt(SX**2 - SX * SY + SY**2 + 3.0 * TXY**2)
            end if
        end do
        !$omp end parallel do

    end subroutine INTEREQ

    !-------------------------------------------------------------------
    !> @brief Computes the derivatives of the quadratic continuous shape
    !>        functions with respect to the intrinsic coordinate QSI.
    !-------------------------------------------------------------------
    subroutine CALC_DFFORMA(QSI, DN1DQSI, DN2DQSI, DN3DQSI)
        implicit none
        real, intent(in)  :: QSI
        real, intent(out) :: DN1DQSI, DN2DQSI, DN3DQSI

        DN1DQSI =  QSI - 0.5
        DN2DQSI = -2.0 * QSI
        DN3DQSI =  QSI + 0.5
    end subroutine CALC_DFFORMA

    !-------------------------------------------------------------------
    !> @brief Computes the boundary stress tensor (sigma_x, tau_xy,
    !>        sigma_y), von Mises equivalent stress, and when NCELL > 0,
    !>        assembles the boundary rows of ASIG_MAT (BDSTRS/BSCOEF) and
    !>        ESIG_MAT (PBSCOEF).
    !-------------------------------------------------------------------
    subroutine BOUND_STRESS()
        implicit none
        integer :: EL, NO, NO1, NO2, NO3, NO_GLOBAL, I
        integer :: ID, M_GLOB, J, COL_H, COL_G, MP, M, N_IDX, K_IDX, L_IDX, NTM3, ROW_S, COL_S, ISG0
        integer :: NODE_COUNT(N)
        real :: X1, Y1, X2, Y2, X3, Y3
        real :: QSI, DN1DQSI, DN2DQSI, DN3DQSI, DN_LOC(3)
        real :: DXDQSI, DYDQSI, DGAMADQSI
        real :: SX_TAN, SY_TAN, NX_VEC, NY_VEC
        real :: COSL(2, 2), LIJ(2, 2), LJI(2, 2), S_LOC(2, 2), SIGMA_GLOB(2, 2)
        real :: TX, TY, S11, S12, S22, E22, DVDQSI
        real :: UNO1(2), UNO2(2), UNO3(2), V1(2), V2(2), V3(2)
        real :: SX, TXY, SY, SZ, TM, TM4, TRM

        SBOU       = 0.0
        NODE_COUNT = 0

        if (NCELL > 0) then
            ASIG_MAT(1:NSIG*N, 1:NN)   = 0.0
            ESIG_MAT(1:NSIG*N, 1:NTSF) = 0.0
        end if

        do EL = 1, NE
            NO1 = CONEC(EL, 1)
            NO2 = CONEC(EL, 2)
            NO3 = CONEC(EL, 3)

            X1 = X(NO1); Y1 = Y(NO1)
            X2 = X(NO2); Y2 = Y(NO2)
            X3 = X(NO3); Y3 = Y(NO3)

            UNO1(1) = FI(2*NO1 - 1); UNO1(2) = FI(2*NO1)
            UNO2(1) = FI(2*NO2 - 1); UNO2(2) = FI(2*NO2)
            UNO3(1) = FI(2*NO3 - 1); UNO3(2) = FI(2*NO3)

            do NO = 1, 3
                QSI = real(NO - 2)
                call CALC_DFFORMA(QSI, DN1DQSI, DN2DQSI, DN3DQSI)

                DXDQSI    = DN1DQSI * X1 + DN2DQSI * X2 + DN3DQSI * X3
                DYDQSI    = DN1DQSI * Y1 + DN2DQSI * Y2 + DN3DQSI * Y3
                DGAMADQSI = sqrt(DXDQSI**2 + DYDQSI**2)

                ! Unit tangent vector (SX_TAN, SY_TAN) and outward normal (NX_VEC, NY_VEC)
                SX_TAN = DXDQSI / DGAMADQSI
                SY_TAN = DYDQSI / DGAMADQSI
                NX_VEC =  SY_TAN
                NY_VEC = -SX_TAN

                ! COSL(1, :) = tangent L1i, COSL(2, :) = normal L2i (matching BEMECH AXES_COS)
                COSL(1, 1) = SX_TAN; COSL(1, 2) = SY_TAN
                COSL(2, 1) = NX_VEC; COSL(2, 2) = NY_VEC

                DN_LOC(1) = DN1DQSI / DGAMADQSI
                DN_LOC(2) = DN2DQSI / DGAMADQSI
                DN_LOC(3) = DN3DQSI / DGAMADQSI

                ! Transformation matrix (global -> local normal-tangential)
                LIJ(1, 1) =  NX_VEC; LIJ(1, 2) = NY_VEC
                LIJ(2, 1) = -NY_VEC; LIJ(2, 2) = NX_VEC

                NO_GLOBAL = CONEC(EL, NO)
                ISG0      = NSIG * (NO_GLOBAL - 1)

                TX = DFI(6*EL - 6 + 2*NO - 1)
                TY = DFI(6*EL - 6 + 2*NO)

                S11 = LIJ(1, 1) * TX + LIJ(1, 2) * TY
                S12 = LIJ(2, 1) * TX + LIJ(2, 2) * TY

                V1 = matmul(LIJ, UNO1)
                V2 = matmul(LIJ, UNO2)
                V3 = matmul(LIJ, UNO3)

                DVDQSI = DN1DQSI * V1(2) + DN2DQSI * V2(2) + DN3DQSI * V3(2)
                E22    = DVDQSI / DGAMADQSI
                S22    = (2.0 * GE / PR1) * E22 + (PR / PR1) * S11

                LJI(1, 1) = NX_VEC; LJI(1, 2) = -NY_VEC
                LJI(2, 1) = NY_VEC; LJI(2, 2) =  NX_VEC

                S_LOC(1, 1) = S11; S_LOC(1, 2) = S12
                S_LOC(2, 1) = S12; S_LOC(2, 2) = S22

                SIGMA_GLOB = matmul(LJI, matmul(S_LOC, transpose(LJI)))

                SBOU(NO_GLOBAL, 1) = SBOU(NO_GLOBAL, 1) + SIGMA_GLOB(1, 1)
                SBOU(NO_GLOBAL, 2) = SBOU(NO_GLOBAL, 2) + SIGMA_GLOB(1, 2)
                SBOU(NO_GLOBAL, 3) = SBOU(NO_GLOBAL, 3) + SIGMA_GLOB(2, 2)
                NODE_COUNT(NO_GLOBAL) = NODE_COUNT(NO_GLOBAL) + 1

                ! Assemble boundary rows of ASIG_MAT and ESIG_MAT (BDSTRS, BSCOEF, PBSCOEF)
                if (NCELL > 0) then
                    ! 1. Displacement sensitivity coefficients (Eqs. 4.86 & 4.89)
                    do ID = 1, 3
                        M_GLOB = CONEC(EL, ID)
                        do J = 1, 2
                            COL_H = 2 * M_GLOB - 2 + J
                            if (FIRST_G_DOF(COL_H) == 0) then
                                MP = 0
                                do M = 1, 2
                                    do N_IDX = M, 2
                                        MP    = MP + 1
                                        ROW_S = INX(MP)
                                        TM    = (1.0 / PR1) * COSL(1, M) * COSL(1, N_IDX) * &
                                                COSL(1, J) * DN_LOC(ID)
                                        ASIG_MAT(ISG0 + ROW_S, COL_H) = ASIG_MAT(ISG0 + ROW_S, COL_H) + TM
                                    end do
                                end do
                                if (NSIG == 4) then
                                    TM4 = (PR / PR1) * COSL(1, J) * DN_LOC(ID)
                                    ASIG_MAT(ISG0 + 4, COL_H) = ASIG_MAT(ISG0 + 4, COL_H) + TM4
                                end if
                            end if
                        end do
                    end do

                    ! 2. Traction sensitivity coefficients (Eqs. 4.87 & 4.90)
                    do J = 1, 2
                        COL_G = 6 * EL - 6 + 2 * NO - 2 + J
                        COL_H = 2 * NO_GLOBAL - 2 + J
                        if (KODE(COL_G) <= 0) then
                            MP = 0
                            do M = 1, 2
                                do N_IDX = M, 2
                                    MP    = MP + 1
                                    ROW_S = INX(MP)
                                    TM    = 0.5 * ((COSL(2, M) * COSL(1, N_IDX) + &
                                                    COSL(1, M) * COSL(2, N_IDX)) * COSL(1, J) + &
                                                   (DLT(M, N_IDX) - (PR2 / PR1) * &
                                                    COSL(1, M) * COSL(1, N_IDX)) * COSL(2, J))
                                    ASIG_MAT(ISG0 + ROW_S, COL_H) = ASIG_MAT(ISG0 + ROW_S, COL_H) + TM
                                end do
                            end do
                            if (NSIG == 4) then
                                TM4 = 0.5 * (PR / PR1) * COSL(2, J)
                                ASIG_MAT(ISG0 + 4, COL_H) = ASIG_MAT(ISG0 + 4, COL_H) + TM4
                            end if
                        end if
                    end do

                    ! 3. Initial plastic stress coefficients on boundary (PBSCOEF, Eqs. 9.19 & 9.20)
                    MP = 0
                    do M = 1, 2
                        do N_IDX = M, 2
                            MP    = MP + 1
                            ROW_S = INX(MP)
                            NTM3  = 0
                            do K_IDX = 1, 2
                                do L_IDX = K_IDX, 2
                                    NTM3  = NTM3 + 1
                                    COL_S = INX(NTM3)
                                    TRM   = COSL(1, M) * COSL(1, N_IDX) * &
                                            (DLT(K_IDX, L_IDX) - (1.0 / PR1) * COSL(2, K_IDX) * COSL(2, L_IDX))
                                    ESIG_MAT(ISG0 + ROW_S, ISG0 + COL_S) = &
                                        ESIG_MAT(ISG0 + ROW_S, ISG0 + COL_S) - TRM * (2.0 - DLT(K_IDX, L_IDX))
                                end do
                            end do
                        end do
                    end do
                    if (NSIG == 4) then
                        NTM3 = 0
                        do K_IDX = 1, 2
                            do L_IDX = K_IDX, 2
                                NTM3  = NTM3 + 1
                                COL_S = INX(NTM3)
                                ESIG_MAT(ISG0 + 4, ISG0 + COL_S) = &
                                    ESIG_MAT(ISG0 + 4, ISG0 + COL_S) + &
                                    (PR / PR1) * COSL(2, K_IDX) * COSL(2, L_IDX) * (2.0 - DLT(K_IDX, L_IDX))
                            end do
                        end do
                        ESIG_MAT(ISG0 + 4, ISG0 + 4) = ESIG_MAT(ISG0 + 4, ISG0 + 4) - 1.0
                    end if
                end if
            end do
        end do

        ! Average at shared extreme nodes and compute von Mises equivalent stress
        do I = 1, N
            if (NODE_COUNT(I) > 0) then
                SBOU(I, 1:3) = SBOU(I, 1:3) / real(NODE_COUNT(I))
                if (NCELL > 0) then
                    ISG0 = NSIG * (I - 1)
                    ASIG_MAT(ISG0+1:ISG0+NSIG, 1:NN)   = ASIG_MAT(ISG0+1:ISG0+NSIG, 1:NN)   / real(NODE_COUNT(I))
                    ESIG_MAT(ISG0+1:ISG0+NSIG, 1:NTSF) = ESIG_MAT(ISG0+1:ISG0+NSIG, 1:NTSF) / real(NODE_COUNT(I))
                end if
            end if
            SX  = SBOU(I, 1)
            TXY = SBOU(I, 2)
            SY  = SBOU(I, 3)
            if (NSIG == 4) then
                SZ = PR * (SX + SY)
                SBOU(I, 4) = sqrt(0.5 * ((SX - SY)**2 + (SY - SZ)**2 + (SZ - SX)**2) + 3.0 * TXY**2)
            else
                SBOU(I, 4) = sqrt(SX**2 - SX * SY + SY**2 + 3.0 * TXY**2)
            end if
        end do

    end subroutine BOUND_STRESS

    !-------------------------------------------------------------------
    !> @brief Computes the S and D stress kernel matrices using Gauss
    !>        quadrature with the sinh transformation (SINHTRANS) to
    !>        evaluate stress at an internal point (XP, YP).
    !-------------------------------------------------------------------
    subroutine SIGMAEQ(XP, YP, X1, Y1, X2, Y2, X3, Y3, D11, D12, D22, S11, S12, S22)
        implicit none
        real, intent(in)  :: XP, YP, X1, Y1, X2, Y2, X3, Y3
        real, intent(out) :: D11(6), D12(6), D22(6), S11(6), S12(6), S22(6)

        real :: FA, AL, A, B, C, D, DE
        real :: DX_EL, DY_EL, NORM_DELEM2, NORM_DELEM
        real :: EET, N1G, N2G, N3G, PSX, PSY, B_DIST, XI, JT, W_EFF
        real :: F1, F2, F3, XCO, YCO, XJ, E1, E2, RA, RD1, RD2, RDN
        integer :: I

        D11 = 0.0; D12 = 0.0; D22 = 0.0
        S11 = 0.0; S12 = 0.0; S22 = 0.0

        FA = PR4
        AL = PR2
        A  = X3 - 2.0 * X2 + X1
        B  = (X3 - X1) / 2.0
        C  = Y3 - 2.0 * Y2 + Y1
        D  = (Y3 - Y1) / 2.0
        DE = 4.0 * PI * PR1

        DX_EL       = X3 - X1
        DY_EL       = Y3 - Y1
        NORM_DELEM2 = DX_EL**2 + DY_EL**2
        NORM_DELEM  = sqrt(NORM_DELEM2)

        EET = 2.0 * (DX_EL * (XP - X1) + DY_EL * (YP - Y1)) / NORM_DELEM2 - 1.0

        N1G = 0.5 * EET * (EET - 1.0)
        N2G = 1.0 - EET**2
        N3G = 0.5 * EET * (1.0 + EET)

        PSX    = N1G * X1 + N2G * X2 + N3G * X3
        PSY    = N1G * Y1 + N2G * Y2 + N3G * Y3
        B_DIST = 2.0 * sqrt((PSX - XP)**2 + (PSY - YP)**2) / NORM_DELEM

        do I = 1, NPG
            call SINHTRANS(GI(I), EET, B_DIST, XI, JT)
            F1    = 0.5 * XI * (XI - 1.0)
            F2    = 1.0 - XI**2
            F3    = 0.5 * XI * (XI + 1.0)
            XJ    = sqrt((XI * A + B)**2 + (XI * C + D)**2)
            E1    =  (XI * C + D) / XJ
            E2    = -(XI * A + B) / XJ
            W_EFF = OME(I) * XJ * JT

            XCO = X1 * F1 + X2 * F2 + X3 * F3
            YCO = Y1 * F1 + Y2 * F2 + Y3 * F3
            RA  = sqrt((XP - XCO)**2 + (YP - YCO)**2)
            RD1 = (XCO - XP) / RA
            RD2 = (YCO - YP) / RA
            RDN = RD1 * E1 + RD2 * E2

            D11(1) = D11(1) + (AL * RD1 + 2.0 * RD1**3) * W_EFF * F1 / (DE * RA)
            D11(2) = D11(2) + (2.0 * RD1**2 * RD2 - AL * RD2) * W_EFF * F1 / (DE * RA)
            D11(3) = D11(3) + (AL * RD1 + 2.0 * RD1**3) * W_EFF * F2 / (DE * RA)
            D11(4) = D11(4) + (2.0 * RD1**2 * RD2 - AL * RD2) * W_EFF * F2 / (DE * RA)
            D11(5) = D11(5) + (AL * RD1 + 2.0 * RD1**3) * W_EFF * F3 / (DE * RA)
            D11(6) = D11(6) + (2.0 * RD1**2 * RD2 - AL * RD2) * W_EFF * F3 / (DE * RA)

            D12(1) = D12(1) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F1 / (DE * RA) * W_EFF
            D12(2) = D12(2) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F1 / (DE * RA) * W_EFF
            D12(3) = D12(3) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F2 / (DE * RA) * W_EFF
            D12(4) = D12(4) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F2 / (DE * RA) * W_EFF
            D12(5) = D12(5) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F3 / (DE * RA) * W_EFF
            D12(6) = D12(6) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F3 / (DE * RA) * W_EFF

            D22(1) = D22(1) + (2.0 * RD1 * RD2**2 - AL * RD1) * F1 / (DE * RA) * W_EFF
            D22(2) = D22(2) + (AL * RD2 + 2.0 * RD2**3) * F1 / (DE * RA) * W_EFF
            D22(3) = D22(3) + (2.0 * RD1 * RD2**2 - AL * RD1) * F2 / (DE * RA) * W_EFF
            D22(4) = D22(4) + (AL * RD2 + 2.0 * RD2**3) * F2 / (DE * RA) * W_EFF
            D22(5) = D22(5) + (2.0 * RD1 * RD2**2 - AL * RD1) * F3 / (DE * RA) * W_EFF
            D22(6) = D22(6) + (AL * RD2 + 2.0 * RD2**3) * F3 / (DE * RA) * W_EFF

            S11(1) = S11(1) + (2.0 * RDN * (AL * RD1 + PR * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * PR * E1 * &
                     RD1**2 + AL * (2.0 * E1 * RD1**2 + 2.0 * E1) - FA * E1) * 2.0 * GE * F1 / (DE * RA**2) * &
                     W_EFF
            S11(2) = S11(2) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * PR * E1 * RD1 * RD2 + &
                     AL * 2.0 * E2 * RD1**2 - FA * E2) * 2.0 * GE * F1 / (DE * RA**2) * W_EFF
            S11(3) = S11(3) + (2.0 * RDN * (AL * RD1 + PR * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * PR * E1 * &
                     RD1**2 + AL * (2.0 * E1 * RD1**2 + 2.0 * E1) - FA * E1) * 2.0 * GE * F2 / (DE * RA**2) * &
                     W_EFF
            S11(4) = S11(4) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * PR * E1 * RD1 * RD2 + &
                     AL * 2.0 * E2 * RD1**2 - FA * E2) * 2.0 * GE * F2 / (DE * RA**2) * W_EFF
            S11(5) = S11(5) + (2.0 * RDN * (AL * RD1 + PR * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * PR * E1 * &
                     RD1**2 + AL * (2.0 * E1 * RD1**2 + 2.0 * E1) - FA * E1) * 2.0 * GE * F3 / (DE * RA**2) * &
                     W_EFF
            S11(6) = S11(6) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * PR * E1 * RD1 * RD2 + &
                     AL * 2.0 * E2 * RD1**2 - FA * E2) * 2.0 * GE * F3 / (DE * RA**2) * W_EFF

            S12(1) = S12(1) + (2.0 * RDN * (PR * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * PR * (E1 * RD2 * &
                     RD1 + E2 * RD1**2) + AL * (2.0 * E1 * RD1 * RD2 + E2)) * 2.0 * GE * F1 / &
                     (DE * RA**2) * W_EFF
            S12(2) = S12(2) + (2.0 * RDN * (PR * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * PR * (E1 * RD2**2 + &
                     E2 * RD1 * RD2) + AL * (2.0 * E2 * RD1 * RD2 + E1)) * 2.0 * GE * F1 / &
                     (DE * RA**2) * W_EFF
            S12(3) = S12(3) + (2.0 * RDN * (PR * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * PR * (E1 * RD2 * &
                     RD1 + E2 * RD1**2) + AL * (2.0 * E1 * RD1 * RD2 + E2)) * 2.0 * GE * F2 / &
                     (DE * RA**2) * W_EFF
            S12(4) = S12(4) + (2.0 * RDN * (PR * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * PR * (E1 * RD2**2 + &
                     E2 * RD1 * RD2) + AL * (2.0 * E2 * RD1 * RD2 + E1)) * 2.0 * GE * F2 / &
                     (DE * RA**2) * W_EFF
            S12(5) = S12(5) + (2.0 * RDN * (PR * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * PR * (E1 * RD2 * &
                     RD1 + E2 * RD1**2) + AL * (2.0 * E1 * RD1 * RD2 + E2)) * 2.0 * GE * F3 / &
                     (DE * RA**2) * W_EFF
            S12(6) = S12(6) + (2.0 * RDN * (PR * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * PR * (E1 * RD2**2 + &
                     E2 * RD1 * RD2) + AL * (2.0 * E2 * RD1 * RD2 + E1)) * 2.0 * GE * F3 / &
                     (DE * RA**2) * W_EFF

            S22(1) = S22(1) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * PR * E2 * RD1 * RD2 + &
                     AL * 2.0 * E1 * RD2**2 - FA * E1) * 2.0 * GE * F1 / (DE * RA**2) * W_EFF
            S22(2) = S22(2) + (2.0 * RDN * (AL * RD2 + 2.0 * PR * RD2 - 4.0 * RD2**3) + 4.0 * PR * E2 * &
                     RD2**2 + AL * (2.0 * E2 * RD2**2 + 2.0 * E2) - FA * E2) * 2.0 * GE * F1 / &
                     (DE * RA**2) * W_EFF
            S22(3) = S22(3) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * PR * E2 * RD1 * RD2 + &
                     AL * 2.0 * E1 * RD2**2 - FA * E1) * 2.0 * GE * F2 / (DE * RA**2) * W_EFF
            S22(4) = S22(4) + (2.0 * RDN * (AL * RD2 + 2.0 * PR * RD2 - 4.0 * RD2**3) + 4.0 * PR * E2 * &
                     RD2**2 + AL * (2.0 * E2 * RD2**2 + 2.0 * E2) - FA * E2) * 2.0 * GE * F2 / &
                     (DE * RA**2) * W_EFF
            S22(5) = S22(5) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * PR * E2 * RD1 * RD2 + &
                     AL * 2.0 * E1 * RD2**2 - FA * E1) * 2.0 * GE * F3 / (DE * RA**2) * W_EFF
            S22(6) = S22(6) + (2.0 * RDN * (AL * RD2 + 2.0 * PR * RD2 - 4.0 * RD2**3) + 4.0 * PR * E2 * &
                     RD2**2 + AL * (2.0 * E2 * RD2**2 + 2.0 * E2) - FA * E2) * 2.0 * GE * F3 / &
                     (DE * RA**2) * W_EFF
        end do

    end subroutine SIGMAEQ

    !===================================================================
    ! 2D ELASTOPLASTIC DOMAIN INTEGRATION & SOLVER ROUTINES (BEMECH 2D)
    !===================================================================

    !-------------------------------------------------------------------
    !> @brief Computes 2D initial stress influence matrices EB_MAT,
    !>        EU_MAT, and ESIG_MAT over internal cells (NODC = 6 or 8).
    !-------------------------------------------------------------------
    subroutine NL_COEFS()
        implicit none
        write(*, '(/10X, "Forming 2D PLASTIC Coefficient Matrices (NCELL=", I4, ", NODC=", I2, ")...")') &
            NCELL, NODC

        ! 1. Initial stress coefficients E_ijk for boundary & interior displacements (NTU = 4)
        call INT_CELL_2D(4, 2, 1, NTP)

        ! 2. Initial stress coefficients E_ijkl for interior stresses (NTU = 5)
        if (L > 0) then
            call INT_CELL_2D(5, NSIG, N + 1, NTP)
        end if
    end subroutine NL_COEFS

    !-------------------------------------------------------------------
    !> @brief Integrates E_ijk (NTU=4) or E_ijkl (NTU=5) over all 2D
    !>        internal cells (NODC = 6 triangles or NODC = 8 quads)
    !>        in parallel using OpenMP.
    !-------------------------------------------------------------------
    subroutine INT_CELL_2D(NTU, NROW, NBEG, NEND)
        implicit none
        integer, intent(in) :: NTU, NROW, NBEG, NEND

        integer :: IPOIN, ICELL, ID, LND, NCS, I, J, K, L_IDX, NP, JP, MP, IP_COL, NGS_SING
        real    :: RLEVL, CP(2), CK(2, NODC), PKERN(6, NSIG * NODC), ROW_COEF(NROW, NTSF), TERM

        RLEVL    = real(2 - 5 + NTU)  ! RLEVL = 1 for E_ijk (NTU=4), 2 for E_ijkl (NTU=5)
        NGS_SING = 6
        if (NAUTO /= 0) NGS_SING = iabs(NAUTO)

        !$omp parallel do default(shared) &
        !$omp private(IPOIN, ICELL, ID, LND, NCS, I, J, K, L_IDX, NP, JP, MP, IP_COL, &
        !$omp         CP, CK, PKERN, ROW_COEF, TERM) &
        !$omp schedule(dynamic)
        do IPOIN = NBEG, NEND
            CP(1:2)               = CD(1:2, IPOIN)
            ROW_COEF(1:NROW, 1:NTSF) = 0.0

            do ICELL = 1, NCELL
                PKERN(1:NROW, 1:NSIG*NODC) = 0.0
                NCS = 0
                do ID = 1, NODC
                    LND = LNDC(ID, ICELL)
                    if (LND == IPOIN) NCS = ID
                    CK(1:2, ID) = CD(1:2, LND)
                end do

                if (NCS == 0) then
                    call ADAPTINT_CELL_2D(NODC, CK, CP, NTU, RLEVL, NGS_SING, PKERN)
                else
                    call SIN_CELL_2D(NODC, NCS, NTU, CK, CP, RLEVL, NGS_SING, PKERN)
                end if

                ! Assemble cell contributions into ROW_COEF
                do ID = 1, NODC
                    NP = NSIG * (LNDC(ID, ICELL) - 1)
                    JP = NSIG * (ID - 1)
                    do I = 1, NROW
                        MP = I
                        if (NTU == 5) MP = INX(I)
                        do J = 1, NSIG
                            ROW_COEF(I, NP + J) = ROW_COEF(I, NP + J) + RE(J) * PKERN(MP, JP + INX(J))
                        end do
                    end do
                end do
            end do

            ! Free terms F_ijkl for interior stresses (NTU == 5, Eq. 9.3)
            if (NTU == 5) then
                IP_COL = NSIG * (IPOIN - 1)
                NP     = 0
                do K = 1, 2
                    do L_IDX = K, 2
                        NP = NP + 1
                        MP = 0
                        do I = 1, 2
                            do J = I, 2
                                MP   = MP + 1
                                TERM = -1.0 / (8.0 * PR1) * (2.0 - DLT(K, L_IDX)) * &
                                       (0.5 * CJUP * (DLT(I, K) * DLT(J, L_IDX) + DLT(I, L_IDX) * DLT(J, K)) + &
                                        (1.0 - GAM * PR) * DLT(I, J) * DLT(K, L_IDX))
                                ROW_COEF(INX(MP), IP_COL + INX(NP)) = &
                                    ROW_COEF(INX(MP), IP_COL + INX(NP)) + TERM
                            end do
                        end do
                        if (NSIG == 4) then
                            ROW_COEF(4, IP_COL + INX(NP)) = &
                                ROW_COEF(4, IP_COL + INX(NP)) + (PR / (2.0 * PR1)) * DLT(K, L_IDX)
                        end if
                    end do
                end do
                if (NSIG == 4) then
                    ROW_COEF(4, IP_COL + 4) = ROW_COEF(4, IP_COL + 4) - 1.0
                end if
            end if

            ! Store assembled rows in EB_MAT / EU_MAT (when NTU == 4) or ESIG_MAT (when NTU == 5)
            if (NTU == 4) then
                if (IPOIN <= N) then
                    EB_MAT(2*IPOIN-1, 1:NTSF) = ROW_COEF(1, 1:NTSF)
                    EB_MAT(2*IPOIN,   1:NTSF) = ROW_COEF(2, 1:NTSF)
                else
                    EU_MAT(2*(IPOIN-N)-1, 1:NTSF) = ROW_COEF(1, 1:NTSF)
                    EU_MAT(2*(IPOIN-N),   1:NTSF) = ROW_COEF(2, 1:NTSF)
                end if
            else if (NTU == 5) then
                ESIG_MAT(NSIG*(IPOIN-1)+1 : NSIG*IPOIN, 1:NTSF) = ROW_COEF(1:NSIG, 1:NTSF)
            end if
        end do
        !$omp end parallel do

    end subroutine INT_CELL_2D

    !-------------------------------------------------------------------
    !> @brief Evaluates shape functions SHAP(1:NODC) and vector
    !>        RI = X(XI) - CP for a 2D internal cell (NODC = 6 or 8)
    !>        or 4-node linear sub-quadrilateral (NODC = 4).
    !-------------------------------------------------------------------
    subroutine SHAPEF_CELL_2D(NODC_IN, XI, CK, CP, SHAP, RI, R2)
        implicit none
        integer, intent(in) :: NODC_IN
        real, intent(in)    :: XI(2), CK(2, NODC_IN), CP(2)
        real, intent(out)   :: SHAP(NODC_IN), RI(2), R2

        integer :: I, ID
        real    :: L1, L2, L3, C1, C2, WL

        if (NODC_IN == 6) then
            ! 6-node quadratic triangle in area coordinates (L1, L2, L3 = 1 - L1 - L2)
            L1      = XI(1)
            L2      = XI(2)
            L3      = 1.0 - L1 - L2
            SHAP(1) = L1 * (2.0 * L1 - 1.0)
            SHAP(2) = L2 * (2.0 * L2 - 1.0)
            SHAP(3) = L3 * (2.0 * L3 - 1.0)
            SHAP(4) = 4.0 * L1 * L2
            SHAP(5) = 4.0 * L2 * L3
            SHAP(6) = 4.0 * L3 * L1
        else
            ! 4-node or 8-node quadrilateral on [-1, 1]^2
            do I = 1, 4
                C1      = CORDL_QUAD(1, I)
                C2      = CORDL_QUAD(2, I)
                SHAP(I) = 0.25 * (1.0 + C1 * XI(1)) * (1.0 + C2 * XI(2))
            end do
            if (NODC_IN == 8) then
                do I = 1, 4
                    C1        = CORDL_QUAD(1, I)
                    C2        = CORDL_QUAD(2, I)
                    SHAP(I)   = SHAP(I) * (C1 * XI(1) + C2 * XI(2) - 1.0)
                    C1        = CORDL_QUAD(1, I + 4)
                    C2        = CORDL_QUAD(2, I + 4)
                    WL        = C1 * XI(1) + C2 * XI(2)
                    SHAP(I+4) = 0.5 * (WL + 1.0) * (1.0 - (C1 * XI(2))**2 - (C2 * XI(1))**2)
                end do
            end if
        end if

        R2 = 0.0
        do I = 1, 2
            RI(I) = -CP(I)
            do ID = 1, NODC_IN
                RI(I) = RI(I) + SHAP(ID) * CK(I, ID)
            end do
            R2 = R2 + RI(I) * RI(I)
        end do
    end subroutine SHAPEF_CELL_2D

    !-------------------------------------------------------------------
    !> @brief Evaluates shape function derivatives DN(2, NODC_IN),
    !>        Jacobian matrix GD(2, 2), and determinant FJCB for a 2D
    !>        cell (NODC_IN = 6, 8, or 4).
    !-------------------------------------------------------------------
    subroutine DSHAPE_CELL_2D(NODC_IN, XI, CK, FJCB, DN, GD)
        implicit none
        integer, intent(in) :: NODC_IN
        real, intent(in)    :: XI(2), CK(2, NODC_IN)
        real, intent(out)   :: FJCB, DN(2, NODC_IN), GD(2, 2)

        integer :: I, J, ID
        real    :: L1, L2, L3, C1, C2, S, T

        if (NODC_IN == 6) then
            L1       = XI(1)
            L2       = XI(2)
            L3       = 1.0 - L1 - L2
            DN(1, 1) =  4.0 * L1 - 1.0;        DN(2, 1) =  0.0
            DN(1, 2) =  0.0;                   DN(2, 2) =  4.0 * L2 - 1.0
            DN(1, 3) = -(4.0 * L3 - 1.0);      DN(2, 3) = -(4.0 * L3 - 1.0)
            DN(1, 4) =  4.0 * L2;              DN(2, 4) =  4.0 * L1
            DN(1, 5) = -4.0 * L2;              DN(2, 5) =  4.0 * (L3 - L2)
            DN(1, 6) =  4.0 * (L3 - L1);       DN(2, 6) = -4.0 * L1
        else
            do I = 1, 4
                C1       = CORDL_QUAD(1, I)
                C2       = CORDL_QUAD(2, I)
                DN(1, I) = 0.25 * C1 * (1.0 + C2 * XI(2))
                DN(2, I) = 0.25 * C2 * (1.0 + C1 * XI(1))
            end do
            if (NODC_IN == 8) then
                do I = 1, 4
                    C1         = CORDL_QUAD(1, I)
                    C2         = CORDL_QUAD(2, I)
                    S          = C1 * XI(1) + C2 * XI(2) - 1.0
                    DN(1, I)   = DN(1, I) * S + 0.25 * (1.0 + C1 * XI(1)) * (1.0 + C2 * XI(2)) * C1
                    DN(2, I)   = DN(2, I) * S + 0.25 * (1.0 + C1 * XI(1)) * (1.0 + C2 * XI(2)) * C2
                    C1         = CORDL_QUAD(1, I + 4)
                    C2         = CORDL_QUAD(2, I + 4)
                    S          = 1.0 + C1 * XI(1) + C2 * XI(2)
                    T          = 1.0 - (C1 * XI(2))**2 - (C2 * XI(1))**2
                    DN(1, I+4) = 0.5 * C1 * T - C2 * C2 * XI(1) * S
                    DN(2, I+4) = 0.5 * C2 * T - C1 * C1 * XI(2) * S
                end do
            end if
        end if

        do I = 1, 2
            do J = 1, 2
                GD(I, J) = 0.0
                do ID = 1, NODC_IN
                    GD(I, J) = GD(I, J) + DN(J, ID) * CK(I, ID)
                end do
            end do
        end do

        FJCB = abs(GD(1, 1) * GD(2, 2) - GD(1, 2) * GD(2, 1))
    end subroutine DSHAPE_CELL_2D

    !-------------------------------------------------------------------
    !> @brief Evaluates 1D line shape functions and outward normal on a
    !>        3-node quadratic cell side [node1(-1), node2(+1), node3(0)].
    !-------------------------------------------------------------------
    subroutine SHAPEF_SIDE_2D(XI1, CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
        implicit none
        real, intent(in)  :: XI1, CSID(2, 3), CP(2)
        real, intent(out) :: SHAP_S(3), RI(2), R2, COSB(2), FJCB

        real :: DN_S(3), GD_S(2)
        integer :: I, ID

        ! In BEMECH, a 3-node side has Node 1 at xi=-1, Node 2 at xi=+1, Node 3 at xi=0
        SHAP_S(1) = -0.5 * XI1 * (1.0 - XI1)
        SHAP_S(2) =  0.5 * XI1 * (1.0 + XI1)
        SHAP_S(3) =  1.0 - XI1 * XI1

        DN_S(1) = -0.5 * (1.0 - 2.0 * XI1)
        DN_S(2) =  0.5 * (1.0 + 2.0 * XI1)
        DN_S(3) = -2.0 * XI1

        R2 = 0.0
        do I = 1, 2
            RI(I)   = -CP(I)
            GD_S(I) = 0.0
            do ID = 1, 3
                RI(I)   = RI(I)   + SHAP_S(ID) * CSID(I, ID)
                GD_S(I) = GD_S(I) + DN_S(ID)   * CSID(I, ID)
            end do
            R2 = R2 + RI(I) * RI(I)
        end do

        FJCB    = sqrt(GD_S(1)**2 + GD_S(2)**2)
        COSB(1) =  GD_S(2) / FJCB
        COSB(2) = -GD_S(1) / FJCB
    end subroutine SHAPEF_SIDE_2D

    !-------------------------------------------------------------------
    !> @brief Adaptive integration over a non-singular 2D cell
    !>        (NODC = 6 triangle or NODC = 8 quadrilateral).
    !-------------------------------------------------------------------
    subroutine ADAPTINT_CELL_2D(NODC_IN, CK, CP, NTU, RLEVL, NGS_SING, PKERN)
        implicit none
        integer, intent(in)    :: NODC_IN, NTU, NGS_SING
        real,    intent(in)    :: CK(2, NODC_IN), CP(2), RLEVL
        real,    intent(inout) :: PKERN(6, NSIG * NODC_IN)

        integer, parameter :: MGAUS = 10
        integer :: ISID, IG, I, J, ISUB1, ISUB2, NGSP1, NGSP2, NSUB1, NSUB2, ID, N1_S, N2_S, N3_S
        real    :: POSG1(MGAUS), WEIT1(MGAUS), POSG2(MGAUS), WEIT2(MGAUS)
        real    :: CSID(2, 3), SHAP_S(3), RI(2), COSB(2), SHAP(NODC_IN), DN(2, NODC_IN), GD(2, 2)
        real    :: SHAP4(4), DN4(2, 4), CSUB(2, 4), XIC(2), XI(2), ETA(2)
        real    :: AL(2), DGS(2), AVL, DISL, R2, FJCB, FJCBL, WJCB
        real    :: LA(2), LB(2), LC(2), CEN_L(2), H_STEP, RHO, S_PAR, AREA2

        ! Compute side lengths of the cell
        call GAUSSV(NGS_SING, POSG1, WEIT1)

        if (NODC_IN == 6) then
            AVL = 0.0
            do ISID = 1, 3
                N1_S = NODEF_TRI(1, ISID)
                N2_S = NODEF_TRI(2, ISID)
                N3_S = NODEF_TRI(3, ISID)
                CSID(:, 1) = CK(:, N1_S)
                CSID(:, 2) = CK(:, N2_S)
                CSID(:, 3) = CK(:, N3_S)
                AL(1) = 0.0
                do IG = 1, NGS_SING
                    call SHAPEF_SIDE_2D(POSG1(IG), CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
                    AL(1) = AL(1) + FJCB * WEIT1(IG)
                end do
                AVL = AVL + AL(1) / 3.0
            end do
            AL(1) = AVL; AL(2) = AVL

            ! Estimate minimum distance DISL from CP to the 6-node triangle
            call MINDIST_TRI_2D(CK, CP, AVL, XIC, DISL)

            ! Choose Gauss order NGSP1 and number of triangle subdivisions NSUB1
            call CHOSEGP_1D(RLEVL, AVL, DISL, 1, MGAUS, NGSP1, NSUB1)

            if (NSUB1 == 1) then
                call GAUSSV(NGSP1, POSG1, WEIT1)
                LA = [1.0, 0.0]; LB = [0.0, 1.0]; LC = [0.0, 0.0]
                do I = 1, NGSP1
                    RHO = 0.5 * (1.0 + POSG1(I))
                    FJCBL = 0.25 * RHO
                    do J = 1, NGSP1
                        S_PAR = 0.5 * (1.0 + POSG1(J))
                        XI(1) = (1.0 - RHO) * LC(1) + RHO * ((1.0 - S_PAR) * LA(1) + S_PAR * LB(1))
                        XI(2) = (1.0 - RHO) * LC(2) + RHO * ((1.0 - S_PAR) * LA(2) + S_PAR * LB(2))
                        call SHAPEF_CELL_2D(6, XI, CK, CP, SHAP, RI, R2)
                        call DSHAPE_CELL_2D(6, XI, CK, FJCB, DN, GD)
                        WJCB = FJCB * FJCBL * WEIT1(I) * WEIT1(J)
                        call EVAL_KE_2D(RI, R2, WJCB, 6, NTU, 0, SHAP, PKERN)
                    end do
                end do
            else
                ! Subdivide the triangle into NSUB1^2 congruent sub-triangles
                H_STEP = 1.0 / real(NSUB1)
                AREA2  = H_STEP * H_STEP
                do ISUB1 = 0, NSUB1 - 1
                    do ISUB2 = 0, NSUB1 - 1 - ISUB1
                        ! 1. Upward sub-triangle
                        LA = [real(ISUB1 + 1) * H_STEP, real(ISUB2)     * H_STEP]
                        LB = [real(ISUB1)     * H_STEP, real(ISUB2 + 1) * H_STEP]
                        LC = [real(ISUB1)     * H_STEP, real(ISUB2)     * H_STEP]
                        CEN_L = (LA + LB + LC) / 3.0
                        call SHAPEF_CELL_2D(6, CEN_L, CK, CP, SHAP, RI, R2)
                        call CHOSEGP_1D(RLEVL, AVL / real(NSUB1), sqrt(R2), 0, MGAUS, NGSP1, NSUB2)
                        call GAUSSV(NGSP1, POSG1, WEIT1)
                        do I = 1, NGSP1
                            RHO   = 0.5 * (1.0 + POSG1(I))
                            FJCBL = 0.25 * RHO * AREA2
                            do J = 1, NGSP1
                                S_PAR = 0.5 * (1.0 + POSG1(J))
                                XI(1) = (1.0 - RHO) * LC(1) + RHO * ((1.0 - S_PAR) * LA(1) + S_PAR * LB(1))
                                XI(2) = (1.0 - RHO) * LC(2) + RHO * ((1.0 - S_PAR) * LA(2) + S_PAR * LB(2))
                                call SHAPEF_CELL_2D(6, XI, CK, CP, SHAP, RI, R2)
                                call DSHAPE_CELL_2D(6, XI, CK, FJCB, DN, GD)
                                WJCB = FJCB * FJCBL * WEIT1(I) * WEIT1(J)
                                call EVAL_KE_2D(RI, R2, WJCB, 6, NTU, 0, SHAP, PKERN)
                            end do
                        end do

                        ! 2. Downward sub-triangle (if inside reference triangle)
                        if (ISUB1 + ISUB2 < NSUB1 - 1) then
                            LA = [real(ISUB1 + 1) * H_STEP, real(ISUB2)     * H_STEP]
                            LB = [real(ISUB1)     * H_STEP, real(ISUB2 + 1) * H_STEP]
                            LC = [real(ISUB1 + 1) * H_STEP, real(ISUB2 + 1) * H_STEP]
                            CEN_L = (LA + LB + LC) / 3.0
                            call SHAPEF_CELL_2D(6, CEN_L, CK, CP, SHAP, RI, R2)
                            call CHOSEGP_1D(RLEVL, AVL / real(NSUB1), sqrt(R2), 0, MGAUS, NGSP1, NSUB2)
                            call GAUSSV(NGSP1, POSG1, WEIT1)
                            do I = 1, NGSP1
                                RHO   = 0.5 * (1.0 + POSG1(I))
                                FJCBL = 0.25 * RHO * AREA2
                                do J = 1, NGSP1
                                    S_PAR = 0.5 * (1.0 + POSG1(J))
                                    XI(1) = (1.0 - RHO) * LC(1) + RHO * ((1.0 - S_PAR) * LA(1) + S_PAR * LB(1))
                                    XI(2) = (1.0 - RHO) * LC(2) + RHO * ((1.0 - S_PAR) * LA(2) + S_PAR * LB(2))
                                    call SHAPEF_CELL_2D(6, XI, CK, CP, SHAP, RI, R2)
                                    call DSHAPE_CELL_2D(6, XI, CK, FJCB, DN, GD)
                                    WJCB = FJCB * FJCBL * WEIT1(I) * WEIT1(J)
                                    call EVAL_KE_2D(RI, R2, WJCB, 6, NTU, 0, SHAP, PKERN)
                                end do
                            end do
                        end if
                    end do
                end do
            end if

        else
            ! 8-node quadrilateral cell (exact BEMECH ADAPTINT logic)
            AVL = 0.0
            do ISID = 1, 2
                AL(ISID) = 0.0
                XI       = 0.0
                do IG = 1, NGS_SING
                    XI(ISID) = POSG1(IG)
                    call DSHAPE_CELL_2D(8, XI, CK, FJCB, DN, GD)
                    FJCB     = sqrt(GD(1, ISID)**2 + GD(2, ISID)**2)
                    AL(ISID) = AL(ISID) + FJCB * WEIT1(IG)
                end do
                AVL = AVL + 0.5 * AL(ISID)
            end do

            call MINDIST_QUAD_2D(8, CK, CP, AVL, XIC, DISL)
            call CHOSEGP_1D(RLEVL, AL(1), DISL, 1, MGAUS, NGSP1, NSUB1)
            call CHOSEGP_1D(RLEVL, AL(2), DISL, 1, MGAUS, NGSP2, NSUB2)

            if (NSUB1 + NSUB2 == 2) then
                call GAUSSV(NGSP1, POSG1, WEIT1)
                call GAUSSV(NGSP2, POSG2, WEIT2)
                do J = 1, NGSP2
                    do I = 1, NGSP1
                        XI = [POSG1(I), POSG2(J)]
                        call SHAPEF_CELL_2D(8, XI, CK, CP, SHAP, RI, R2)
                        call DSHAPE_CELL_2D(8, XI, CK, FJCB, DN, GD)
                        WJCB = FJCB * WEIT1(I) * WEIT2(J)
                        call EVAL_KE_2D(RI, R2, WJCB, 8, NTU, 0, SHAP, PKERN)
                    end do
                end do
            else
                DGS(1) = 2.0 / real(NSUB1)
                DGS(2) = 2.0 / real(NSUB2)
                AL(1)  = AL(1) / real(NSUB1)
                AL(2)  = AL(2) / real(NSUB2)
                do ID = 1, 4
                    CSUB(2, ID) = -1.0 + (CORDL_QUAD(2, ID) + 1.0) / real(NSUB2)
                end do
                do ISUB2 = 1, NSUB2
                    do ID = 1, 4
                        CSUB(1, ID) = -1.0 + (CORDL_QUAD(1, ID) + 1.0) / real(NSUB1)
                    end do
                    do ISUB1 = 1, NSUB1
                        call MINDIST_QUAD_2D(4, CSUB, XIC, 1.0, XI, DISL)
                        call SHAPEF_CELL_2D(4, XI, CSUB, [0.0, 0.0], SHAP4, RI, R2)
                        call SHAPEF_CELL_2D(8, RI, CK, CP, SHAP, XI, R2)
                        call CHOSEGP_1D(RLEVL, AL(1), sqrt(R2), 0, MGAUS, NGSP1, I)
                        call CHOSEGP_1D(RLEVL, AL(2), sqrt(R2), 0, MGAUS, NGSP2, I)
                        call GAUSSV(NGSP1, POSG1, WEIT1)
                        call GAUSSV(NGSP2, POSG2, WEIT2)
                        do J = 1, NGSP2
                            do I = 1, NGSP1
                                ETA = [POSG1(I), POSG2(J)]
                                call SHAPEF_CELL_2D(4, ETA, CSUB, [0.0, 0.0], SHAP4, XI, R2)
                                call DSHAPE_CELL_2D(4, ETA, CSUB, FJCBL, DN4, GD)
                                call SHAPEF_CELL_2D(8, XI, CK, CP, SHAP, RI, R2)
                                call DSHAPE_CELL_2D(8, XI, CK, FJCB, DN, GD)
                                WJCB = FJCB * FJCBL * WEIT1(I) * WEIT2(J)
                                call EVAL_KE_2D(RI, R2, WJCB, 8, NTU, 0, SHAP, PKERN)
                            end do
                        end do
                        CSUB(1, 1:4) = CSUB(1, 1:4) + DGS(1)
                    end do
                    CSUB(2, 1:4) = CSUB(2, 1:4) + DGS(2)
                end do
            end if
        end if

    end subroutine ADAPTINT_CELL_2D

    !-------------------------------------------------------------------
    !> @brief Singular cell integration (NCS > 0) for 2D cells
    !>        (NODC = 6 triangles or NODC = 8 quadrilaterals) using
    !>        degenerate sub-triangle subdivision and Eq. (9.10).
    !-------------------------------------------------------------------
    subroutine SIN_CELL_2D(NODC_IN, NCS, NTU, CK, CP, RLEVL, NGS_SING, PKERN)
        implicit none
        integer, intent(in)    :: NODC_IN, NCS, NTU, NGS_SING
        real,    intent(in)    :: CK(2, NODC_IN), CP(2), RLEVL
        real,    intent(inout) :: PKERN(6, NSIG * NODC_IN)

        integer :: MSUB, ISN, ISUB, ID, I, J, N1_S, N2_S, N3_S
        integer :: NSB(4)
        real    :: POSG(10), WEIT(10), CSID(2, 3), CSUB(2, 4), XIS(2), LA(2), LB(2), LC(2)
        real    :: ETA(2), XI(2), RI(2), SHAP4(4), DN4(2, 4), SHAP(NODC_IN), DN(2, NODC_IN), GD(2, 2)
        real    :: R2, FJCB, FJCBL, WJCB, RHO, S_PAR, AREA2

        ISN = NCS
        if (NTU == 4) ISN = 0

        call GAUSSV(NGS_SING, POSG, WEIT)

        if (NODC_IN == 6) then
            MSUB = 3
            LC   = CORDL_TRI(1:2, NCS)
            do ISUB = 1, MSUB
                NSB(ISUB) = 1
                do ID = 1, 3
                    if (NODEF_TRI(ID, ISUB) == NCS) NSB(ISUB) = 0
                end do
                if (NSB(ISUB) == 0) cycle

                N1_S = NODEF_TRI(1, ISUB)
                N2_S = NODEF_TRI(2, ISUB)
                N3_S = NODEF_TRI(3, ISUB)

                ! Second integral in Eq. (9.10): line integral over non-adjacent side
                if (NTU == 5) then
                    CSID(1:2, 1) = CK(1:2, N1_S)
                    CSID(1:2, 2) = CK(1:2, N2_S)
                    CSID(1:2, 3) = CK(1:2, N3_S)
                    call CELL_BOUND_2D(NODC_IN, NTU, ISN, CP, CSID, RLEVL, NGS_SING, PKERN)
                end if

                ! First integral in Eq. (9.10): polar/Duffy degenerate sub-triangle collapsing at LC
                LA    = CORDL_TRI(1:2, N1_S)
                LB    = CORDL_TRI(1:2, N2_S)
                AREA2 = abs((LA(1) - LC(1)) * (LB(2) - LC(2)) - (LA(2) - LC(2)) * (LB(1) - LC(1)))

                do I = 1, NGS_SING
                    RHO   = 0.5 * (1.0 + POSG(I))
                    FJCBL = 0.25 * RHO * AREA2
                    do J = 1, NGS_SING
                        S_PAR = 0.5 * (1.0 + POSG(J))
                        XI(1) = (1.0 - RHO) * LC(1) + RHO * ((1.0 - S_PAR) * LA(1) + S_PAR * LB(1))
                        XI(2) = (1.0 - RHO) * LC(2) + RHO * ((1.0 - S_PAR) * LA(2) + S_PAR * LB(2))
                        call SHAPEF_CELL_2D(6, XI, CK, CP, SHAP, RI, R2)
                        call DSHAPE_CELL_2D(6, XI, CK, FJCB, DN, GD)
                        WJCB = FJCB * FJCBL * WEIT(I) * WEIT(J)
                        call EVAL_KE_2D(RI, R2, WJCB, 6, NTU, ISN, SHAP, PKERN)
                    end do
                end do
            end do

        else
            ! 8-node quadrilateral cell (exact BEMECH SIN_CELL + SETDSUB logic)
            MSUB = 4
            XIS  = CORDL_QUAD(1:2, NCS)
            do ISUB = 1, MSUB
                NSB(ISUB) = 1
                do ID = 1, 3
                    if (NODEF_QUAD(ID, ISUB) == NCS) NSB(ISUB) = 0
                end do
                if (NSB(ISUB) == 0) cycle

                N1_S = NODEF_QUAD(1, ISUB)
                N2_S = NODEF_QUAD(2, ISUB)
                N3_S = NODEF_QUAD(3, ISUB)

                if (NTU == 5) then
                    CSID(1:2, 1) = CK(1:2, N1_S)
                    CSID(1:2, 2) = CK(1:2, N2_S)
                    CSID(1:2, 3) = CK(1:2, N3_S)
                    call CELL_BOUND_2D(NODC_IN, NTU, ISN, CP, CSID, RLEVL, NGS_SING, PKERN)
                end if

                ! SETDSUB for 2D quad: merge nodes not on side ISUB into XIS
                do ID = 1, 4
                    if (ID == N1_S .or. ID == N2_S) then
                        CSUB(1:2, ID) = CORDL_QUAD(1:2, ID)
                    else
                        CSUB(1:2, ID) = XIS(1:2)
                    end if
                end do

                do J = 1, NGS_SING
                    do I = 1, NGS_SING
                        ETA = [POSG(I), POSG(J)]
                        call SHAPEF_CELL_2D(4, ETA, CSUB, [0.0, 0.0], SHAP4, XI, R2)
                        call DSHAPE_CELL_2D(4, ETA, CSUB, FJCBL, DN4, GD)
                        call SHAPEF_CELL_2D(8, XI, CK, CP, SHAP, RI, R2)
                        call DSHAPE_CELL_2D(8, XI, CK, FJCB, DN, GD)
                        WJCB = FJCB * FJCBL * WEIT(I) * WEIT(J)
                        call EVAL_KE_2D(RI, R2, WJCB, 8, NTU, ISN, SHAP, PKERN)
                    end do
                end do
            end do
        end if

    end subroutine SIN_CELL_2D

    !-------------------------------------------------------------------
    !> @brief Evaluates the boundary line integral (second integral in
    !>        Eq. 9.10) over a 3-node quadratic side CSID of a singular
    !>        cell when NTU == 5.
    !-------------------------------------------------------------------
    subroutine CELL_BOUND_2D(NODC_IN, NTU, NCS, CP, CSID, RLEVL, NGS_SING, PKERN)
        implicit none
        integer, intent(in)    :: NODC_IN, NTU, NCS, NGS_SING
        real,    intent(in)    :: CP(2), CSID(2, 3), RLEVL
        real,    intent(inout) :: PKERN(6, NSIG * NODC_IN)

        integer, parameter :: MGAUS = 10
        integer :: IG, NGSP1, NSUB1, ISUB1, DUMMY_SUB
        real    :: POSG(MGAUS), WEIT(MGAUS), SHAP_S(3), RI(2), COSB(2), DUMMY_SHAP(1)
        real    :: AL_SIDE, DISL, R2, FJCB, RN, WJCB, XI1, DGS, X_START, X_MID

        ! Compute side length AL_SIDE and minimum distance DISL from CP to CSID
        call GAUSSV(NGS_SING, POSG, WEIT)
        AL_SIDE = 0.0
        DISL    = 1.0e20
        do IG = 1, NGS_SING
            call SHAPEF_SIDE_2D(POSG(IG), CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
            AL_SIDE = AL_SIDE + FJCB * WEIT(IG)
            if (sqrt(R2) < DISL) DISL = sqrt(R2)
        end do

        call CHOSEGP_1D(RLEVL, AL_SIDE, DISL, 1, MGAUS, NGSP1, NSUB1)

        if (NSUB1 == 1) then
            call GAUSSV(NGSP1, POSG, WEIT)
            do IG = 1, NGSP1
                call SHAPEF_SIDE_2D(POSG(IG), CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
                RN   = RI(1) * COSB(1) + RI(2) * COSB(2)
                WJCB = FJCB * WEIT(IG) * RN
                if (NTU == 5) WJCB = WJCB * log(sqrt(R2))
                call EVAL_KE_2D(RI, R2, WJCB, 1, NTU, NCS, DUMMY_SHAP, PKERN)
            end do
        else
            DGS     = 2.0 / real(NSUB1)
            X_START = -1.0
            do ISUB1 = 1, NSUB1
                X_MID = X_START + 0.5 * DGS
                call SHAPEF_SIDE_2D(X_MID, CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
                call CHOSEGP_1D(RLEVL, AL_SIDE / real(NSUB1), sqrt(R2), 0, MGAUS, NGSP1, DUMMY_SUB)
                call GAUSSV(NGSP1, POSG, WEIT)
                do IG = 1, NGSP1
                    XI1 = X_MID + 0.5 * DGS * POSG(IG)
                    call SHAPEF_SIDE_2D(XI1, CSID, CP, SHAP_S, RI, R2, COSB, FJCB)
                    RN   = RI(1) * COSB(1) + RI(2) * COSB(2)
                    WJCB = 0.5 * DGS * FJCB * WEIT(IG) * RN
                    if (NTU == 5) WJCB = WJCB * log(sqrt(R2))
                    call EVAL_KE_2D(RI, R2, WJCB, 1, NTU, NCS, DUMMY_SHAP, PKERN)
                end do
                X_START = X_START + DGS
            end do
        end if

    end subroutine CELL_BOUND_2D

    !-------------------------------------------------------------------
    !> @brief Evaluates the 2D initial stress kernels E_ikl (NTU=4) and
    !>        E_ijkl (NTU=5) and accumulates them into PKERN.
    !-------------------------------------------------------------------
    subroutine EVAL_KE_2D(RI, R2, WJCB, NODD, NTU, NCS, SHAP, PKERN)
        implicit none
        integer, intent(in)    :: NODD, NTU, NCS
        real,    intent(in)    :: RI(2), R2, WJCB, SHAP(*)
        real,    intent(inout) :: PKERN(6, NSIG * NODC)

        integer :: K, L_IDX, I, J, ID, IP, NP, MP
        real    :: R1, CONP, Z(2), VKE(2, 2, 2, 2), TM1, TM2, SHP

        R1   = sqrt(R2)
        CONP = CON * WJCB / R1
        Z    = RI / R1

        do K = 1, 2
            do L_IDX = K, 2
                do I = 1, 2
                    do J = I, 2
                        if (NTU == 4) then
                            ! Initial stress kernel E_ikl for displacements
                            TM1 = PR2 * (DLT(I, K) * Z(L_IDX) + DLT(I, L_IDX) * Z(K)) - DLT(K, L_IDX) * Z(I)
                            VKE(I, I, K, L_IDX) = CONP * (TM1 + 2.0 * Z(I) * Z(K) * Z(L_IDX))
                            exit
                        else
                            ! Initial stress kernel E_ijkl for interior stresses
                            TM1 = PR2 * (DLT(I, K) * DLT(J, L_IDX) + DLT(I, L_IDX) * DLT(J, K) - &
                                         DLT(I, J) * DLT(K, L_IDX) + 2.0 * DLT(I, J) * Z(K) * Z(L_IDX)) + &
                                  2.0 * DLT(K, L_IDX) * Z(I) * Z(J)
                            TM2 = 2.0 * PR * (DLT(I, K) * Z(J) * Z(L_IDX) + DLT(J, K) * Z(I) * Z(L_IDX) + &
                                              DLT(I, L_IDX) * Z(J) * Z(K) + DLT(J, L_IDX) * Z(I) * Z(K)) - &
                                  2.0 * GAM * Z(I) * Z(J) * Z(K) * Z(L_IDX)
                            VKE(I, J, K, L_IDX) = -(CONP / R1) * (TM1 + TM2)
                        end if
                    end do
                end do
            end do
        end do

        do ID = 1, NODD
            if (NODD == 1) then
                IP  = NSIG * (NCS - 1)
                SHP = 1.0
            else
                IP  = NSIG * (ID - 1)
                SHP = SHAP(ID)
                if (ID == NCS) SHP = SHAP(ID) - 1.0
            end if

            NP = 0
            do K = 1, 2
                do L_IDX = K, 2
                    NP = NP + 1
                    MP = 0
                    do I = 1, 2
                        do J = I, 2
                            MP = MP + 1
                            PKERN(MP, IP + NP) = PKERN(MP, IP + NP) + VKE(I, J, K, L_IDX) * SHP
                            if (NTU == 4) exit
                        end do
                    end do
                    if (NSIG == 4 .and. NTU == 5) then
                        PKERN(4, IP + NP) = PKERN(4, IP + NP) + &
                                            PR * (VKE(1, 1, K, L_IDX) + VKE(2, 2, K, L_IDX)) * SHP
                    end if
                end do
            end do
        end do

    end subroutine EVAL_KE_2D

    !-------------------------------------------------------------------
    !> @brief Determines Gauss integration order NGSP and number of
    !>        sub-elements NSUB using Eqs. (4.28)-(4.32) of Gao & Davies.
    !-------------------------------------------------------------------
    subroutine CHOSEGP_1D(RLEVL, AL_IN, DISL_IN, NFLAG, MGAUS, NGSP, NSUB)
        implicit none
        real,    intent(in)  :: RLEVL, AL_IN, DISL_IN
        integer, intent(in)  :: NFLAG, MGAUS
        integer, intent(out) :: NGSP, NSUB

        real :: PB, WFA, ALM, DISL, ALI

        NSUB = 1
        if (NAUTO > 0) then
            NGSP = NAUTO
            return
        end if

        DISL = max(DISL_IN, 1.0e-12)
        PB   = sqrt(RLEVL * 2.0 / 3.0 + 0.4)
        WFA  = PB * log(abs(TOLGP) / 2.0)

        if (TOLGP >= 0.0) then
            ALM = AL_IN
            if (ALM > 3.9 * DISL) ALM = 3.9 * DISL
            NGSP = int(0.5 * WFA / log(ALM / DISL / 4.0))
        else
            NGSP = int(-WFA / 10.0 * ((8.0 / 3.0 * AL_IN / DISL)**(3.0 / 4.0) + 1.0))
        end if

        if (NGSP < 2) NGSP = 2

        if (NFLAG == 1 .and. NAUTO < 0) then
            if (NGSP > MGAUS) then
                NGSP = MGAUS
                ALI  = AL_IN
                if (TOLGP < 0.0) ALI = 3.0 / 8.0 * DISL * (-10.0 * real(NGSP) / WFA - 1.0)**(4.0 / 3.0)
                if (TOLGP > 0.0) ALI = 4.0 * DISL * (TOLGP / 2.0)**(0.5 * PB / real(NGSP))
                NSUB = int(AL_IN / ALI + 0.95)
                if (NSUB < 1) NSUB = 1
                if (NSUB > 12) NSUB = 12
            end if
        else
            if (NGSP > MGAUS) NGSP = MGAUS
        end if
    end subroutine CHOSEGP_1D

    !-------------------------------------------------------------------
    !> @brief Finds the minimum distance DISL from point CP to a 6-node
    !>        quadratic triangle CK(2, 6).
    !-------------------------------------------------------------------
    subroutine MINDIST_TRI_2D(CK, CP, AVL, XIC, DISL)
        implicit none
        real, intent(in)  :: CK(2, 6), CP(2), AVL
        real, intent(out) :: XIC(2), DISL

        integer :: I, J, NGRID
        real    :: L1, L2, SHAP(6), RI(2), R2, RM

        NGRID = 6
        DISL  = AVL * 1000.0
        XIC   = [1.0 / 3.0, 1.0 / 3.0]

        do I = 0, NGRID
            L1 = real(I) / real(NGRID)
            do J = 0, NGRID - I
                L2 = real(J) / real(NGRID)
                call SHAPEF_CELL_2D(6, [L1, L2], CK, CP, SHAP, RI, R2)
                RM = sqrt(R2)
                if (RM < DISL) then
                    DISL = RM
                    XIC  = [L1, L2]
                end if
            end do
        end do
    end subroutine MINDIST_TRI_2D

    !-------------------------------------------------------------------
    !> @brief Finds the minimum distance DISL from point CP to a 4- or
    !>        8-node quadrilateral CF(2, NODE_IN) via Newton iteration.
    !-------------------------------------------------------------------
    subroutine MINDIST_QUAD_2D(NODE_IN, CF, CP, AVL, XIC, DISL)
        implicit none
        integer, intent(in) :: NODE_IN
        real, intent(in)    :: CF(2, NODE_IN), CP(2), AVL
        real, intent(out)   :: XIC(2), DISL

        integer :: ITER, I
        real    :: TOLE, XI(2), RI(2), SHAP(NODE_IN), DN(2, NODE_IN), GD(2, 2), GDT(2, 2)
        real    :: RQ2, RM, FJCB

        TOLE = 1.0e-6
        XI   = 0.0
        DISL = AVL * 1000.0
        XIC  = 0.0

        do ITER = 1, 20
            call SHAPEF_CELL_2D(NODE_IN, XI, CF, CP, SHAP, RI, RQ2)
            RM = sqrt(RQ2)
            if (abs(RM - DISL) / AVL < TOLE) return
            if (RM < DISL) then
                DISL = RM
                XIC  = XI
            else
                return
            end if

            call DSHAPE_CELL_2D(NODE_IN, XI, CF, FJCB, DN, GD)
            FJCB = GD(1, 1) * GD(2, 2) - GD(1, 2) * GD(2, 1)
            if (abs(FJCB) < 1.0e-12) return
            GDT(1, 1) =  GD(2, 2) / FJCB
            GDT(2, 2) =  GD(1, 1) / FJCB
            GDT(1, 2) = -GD(1, 2) / FJCB
            GDT(2, 1) = -GD(2, 1) / FJCB

            do I = 1, 2
                XI(I) = XI(I) - (GDT(I, 1) * RI(1) + GDT(I, 2) * RI(2))
                if (abs(XI(I)) > 1.0) XI(I) = sign(1.0, XI(I))
            end do
        end do
    end subroutine MINDIST_QUAD_2D

    !-------------------------------------------------------------------
    !> @brief Returns 1D Gauss-Legendre points GP and weights GW for
    !>        orders NGAUS = 2 to 10 (from BEMECH GAUSSV).
    !-------------------------------------------------------------------
    subroutine GAUSSV(NGAUS, GP, GW)
        implicit none
        integer, intent(in) :: NGAUS
        real, intent(out)   :: GP(NGAUS), GW(NGAUS)
        integer :: KGAUS, IGASH, JGASH

        select case (NGAUS)
        case (2)
            GP(1) = -0.5773502691896258; GW(1) = 1.0000000000000000
        case (3)
            GP(1) = -0.7745966692414834; GP(2) = 0.0000000000000000
            GW(1) =  0.5555555555555556; GW(2) = 0.8888888888888889
        case (4)
            GP(1) = -0.8611363115940526; GP(2) = -0.3399810435848563
            GW(1) =  0.3478548451374539; GW(2) =  0.6521451548625461
        case (5)
            GP(1) = -0.9061798459386640; GP(2) = -0.5384693101056831; GP(3) = 0.0
            GW(1) =  0.2369268850561891; GW(2) =  0.4786286704993665; GW(3) = 0.5688888888888889
        case (6)
            GP(1) = -0.9324695142031520; GP(2) = -0.6612093864662645; GP(3) = -0.2386191860831969
            GW(1) =  0.1713244923791703; GW(2) =  0.3607615730481386; GW(3) =  0.4679139345726910
        case (7)
            GP(1) = -0.9491079123427585; GP(2) = -0.7415311855993944
            GP(3) = -0.4058451513773972; GP(4) =  0.0000000000000000
            GW(1) =  0.1294849661688697; GW(2) =  0.2797053914892767
            GW(3) =  0.3818300505051189; GW(4) =  0.4179591836734694
        case (8)
            GP(1) = -0.9602898564975362; GP(2) = -0.7966664774136267
            GP(3) = -0.5255324099163290; GP(4) = -0.1834346424956498
            GW(1) =  0.1012285362903763; GW(2) =  0.2223810344533745
            GW(3) =  0.3137066458778873; GW(4) =  0.3626837833783620
        case (9)
            GP(1) = -0.9681602395076261; GP(2) = -0.8360311073266358
            GP(3) = -0.6133714327005904; GP(4) = -0.3242534234038089; GP(5) = 0.0
            GW(1) =  0.0812743883615744; GW(2) =  0.1806481606948574
            GW(3) =  0.2606106964029355; GW(4) =  0.3123470770400028; GW(5) = 0.3302393550012598
        case (10)
            GP(1) = -0.9739065285171717; GP(2) = -0.8650633666889845
            GP(3) = -0.6794095682990244; GP(4) = -0.4333953941292472; GP(5) = -0.1488743389816312
            GW(1) =  0.0666713443086881; GW(2) =  0.1494513491505806
            GW(3) =  0.2190863625159820; GW(4) =  0.2692667193099963; GW(5) =  0.2955242247147529
        end select

        KGAUS = NGAUS / 2
        do IGASH = 1, KGAUS
            JGASH     = NGAUS + 1 - IGASH
            GP(JGASH) = -GP(IGASH)
            GW(JGASH) =  GW(IGASH)
        end do
    end subroutine GAUSSV

    !-------------------------------------------------------------------
    !> @brief Solves the non-linear elastoplastic system of equations
    !>        by incremental plastic multiplier iteration (NL_SOLVE).
    !-------------------------------------------------------------------
    subroutine NL_SOLVE()
        implicit none
        integer :: IP, I, INCRD, EL, NODE_LOC, NO_GLOB, DIR, J, COL_G, COL_H, K
        integer :: FIRST_DOF_LOC(NN)
        real    :: EFSIG, YIELD, DUMMY_DF(NSIG), TOL, FMULT, FSUM, FSCAL, WOVER, CH
        real    :: STRES_TOT(NSIG), SX, SY, TXY, EQSTR
        real, allocatable :: YC_B(:), YC_U(:), YE_0(:,:), DFI_ORIG(:)

        write(*, '(/10X, "Solving 2D Non-Linear Elastoplastic System Equations...")')

        ! Save original prescribed boundary condition vector DFI before plastic updates
        allocate(DFI_ORIG(NX1))
        DFI_ORIG = 0.0
        do EL = 1, NE
            do NODE_LOC = 1, 3
                NO_GLOB = CONEC(EL, NODE_LOC)
                do DIR = 1, 2
                    J     = 2 * NODE_LOC - 2 + DIR
                    COL_G = 6 * EL - 6 + J
                    COL_H = 2 * NO_GLOB - 2 + DIR
                    if (KODE(COL_G) <= 0) then
                        DFI_ORIG(COL_G) = FI(COL_H)
                    else
                        DFI_ORIG(COL_G) = DFI(COL_G)
                    end if
                end do
            end do
        end do

        ! 1. Find maximum equivalent elastic stress across all NTP nodes
        EFSIG = 0.0
        do IP = 1, NTP
            call PL_FLOW(ESTRS(:, IP), YIELD, DUMMY_DF, 0)
            if (EFSIG < YIELD) EFSIG = YIELD
        end do

        write(IPR, 10) 2.0 * GE * UNIAX, 2.0 * GE * EFSIG
        write(*,   10) 2.0 * GE * UNIAX, 2.0 * GE * EFSIG
10      format(/2X, 'STRESS LIMIT == ', 1PE13.6, 5X, 'MAXIMUM EQU. STRESS == ', 1PE13.6)

        if (EFSIG <= UNIAX) then
            write(*,   '(/10X, "All nodes are in elastic state (EFSIG <= UNIAX).")')
            write(IPR, '(/2X, "ALL NODES ARE IN ELASTIC STATE.")')
            deallocate(DFI_ORIG)
            return
        end if

        ! 2. Determine number of load increments
        INCRD = 6 * int(sqrt(40.0 * (EFSIG / UNIAX - 1.0) + 1.0))
        write(IPR, '(/7X, "DEFAULT INCREMENT == ", I5)') INCRD
        if (NINCS == 0) NINCS = INCRD
        write(IPR, '(7X, "ACTUAL INCREMENTS == ", I5)') NINCS

        TOL   = 0.001 * UNIAX
        FMULT = 1.0
        FSUM  = 1.0
        FSCAL = UNIAX / EFSIG

        do I = 1, iabs(NINCS) - 1
            FMULT = FMULT * FINCR
            FSUM  = FSUM + FMULT
        end do
        WOVER = (1.0 - FSCAL) / FSUM

        ! Store first increment elastic increments
        allocate(YC_B(NN), YC_U(max(1, 2*L)), YE_0(NSIG, NTP))
        YC_B = X_B * WOVER
        if (L > 0) YC_U(1:2*L) = DSOL(1:2*L) * WOVER
        YE_0 = ESTRS * WOVER

        ! Scale state to initial yield load factor FSCAL
        X_B = X_B * FSCAL
        if (L > 0) DSOL(1:2*L) = DSOL(1:2*L) * FSCAL
        ESTRS = ESTRS * FSCAL

        ! 3. Form coupled plastic influence matrices AC_MAT, EUI_MAT, E_MAT
        call FORM_PLASTIC_MATRICES()

        ! 4. Incremental-iterative plastic multiplier loop
        call P_M_ITER(TOL, FSCAL, WOVER, YC_B, YC_U, YE_0)

        ! 5. Reconstruct final boundary displacements FI and tractions DFI from X_B
        FI(1:NN)    = X_B(1:NN)
        DFI(1:6*NE) = DFI_ORIG(1:6*NE)
        FIRST_DOF_LOC = 0
        do EL = 1, NE
            do NODE_LOC = 1, 3
                NO_GLOB = CONEC(EL, NODE_LOC)
                do DIR = 1, 2
                    J     = 2 * NODE_LOC - 2 + DIR
                    COL_G = 6 * EL - 6 + J
                    COL_H = 2 * NO_GLOB - 2 + DIR
                    if (KODE(COL_G) <= 0) then
                        if (FIRST_DOF_LOC(COL_H) == 0) then
                            FIRST_DOF_LOC(COL_H) = COL_G
                            CH         = FI(COL_H) * GE
                            FI(COL_H)  = DFI(COL_G)
                            DFI(COL_G) = CH
                        else
                            DFI(COL_G) = DFI(FIRST_DOF_LOC(COL_H))
                        end if
                    end if
                end do
            end do
        end do

        ! 6. Store final physical stresses (ESTRS + BACKS)*2*GE in SBOU and SSOL/SVM_INT
        do I = 1, N
            STRES_TOT  = (ESTRS(:, I) + BACKS(:, I)) * (2.0 * GE)
            SX         = STRES_TOT(1)
            SY         = STRES_TOT(2)
            TXY        = STRES_TOT(3)
            SBOU(I, 1) = SX
            SBOU(I, 2) = TXY
            SBOU(I, 3) = SY
            call PL_FLOW(ESTRS(:, I) + BACKS(:, I), EQSTR, DUMMY_DF, 0)
            SBOU(I, 4) = EQSTR * (2.0 * GE)
        end do

        do K = 1, L
            IP            = N + K
            STRES_TOT     = (ESTRS(:, IP) + BACKS(:, IP)) * (2.0 * GE)
            SX            = STRES_TOT(1)
            SY            = STRES_TOT(2)
            TXY           = STRES_TOT(3)
            SSOL(3*K - 2) = SX
            SSOL(3*K - 1) = TXY
            SSOL(3*K)     = SY
            call PL_FLOW(ESTRS(:, IP) + BACKS(:, IP), EQSTR, DUMMY_DF, 0)
            SVM_INT(K)    = EQSTR * (2.0 * GE)
        end do

        deallocate(YC_B, YC_U, YE_0, DFI_ORIG)
    end subroutine NL_SOLVE

    !-------------------------------------------------------------------
    !> @brief Forms AC_MAT = AB_MAT^-1 * EB_MAT,
    !>        EUI_MAT = EU_MAT + AU_MAT * AC_MAT, and
    !>        E_MAT   = ESIG_MAT + ASIG_MAT * AC_MAT.
    !-------------------------------------------------------------------
    subroutine FORM_PLASTIC_MATRICES()
        implicit none

        ! Solve AB_MAT * AC_MAT = EB_MAT via LAPACK DGESV
        AC_MAT = EB_MAT
        call SLNPD_MULTI(AB_MAT, AC_MAT, NN, NTSF, IPR)

        ! EUI_MAT = EU_MAT + AU_MAT * AC_MAT
        if (L > 0) then
            EUI_MAT = EU_MAT + matmul(AU_MAT, AC_MAT)
        end if

        ! E_MAT = ESIG_MAT + ASIG_MAT * AC_MAT
        E_MAT = ESIG_MAT + matmul(ASIG_MAT, AC_MAT)
    end subroutine FORM_PLASTIC_MATRICES

    !-------------------------------------------------------------------
    !> @brief Incremental-iterative plastic multiplier solver (P_M_ITER).
    !-------------------------------------------------------------------
    subroutine P_M_ITER(TOL, FSCAL, WOVER, YC_B, YC_U, YE_0)
        implicit none
        real, intent(in)    :: TOL, WOVER, YC_B(NN), YC_U(*), YE_0(NSIG, NTP)
        real, intent(inout) :: FSCAL

        integer :: INCRS, IP, JP, I, J, ITER, NLAMD, IL, JL, J0
        real    :: FMULT, YIELD, PRER, PSI, PSIB, HEPSI, CLAMD, RNORM, RERR, D_DET
        real    :: EQUP2, DPLAS, CBACK, SS, DF_IP(NSIG), DF_JP(NSIG)
        real, allocatable :: YE(:,:), DLAMD(:), A_LAMD(:,:)

        allocate(YE(NSIG, NTP), DLAMD(NTP))
        FMULT = 1.0

        do INCRS = 1, iabs(NINCS)
            FSCAL = FSCAL + WOVER * FMULT

            ! Scale trial stress to yield surface and identify yielding nodes
            call SIG_SCALE(TOL, FMULT, YC_B, YC_U, YE_0, YE, NLAMD)

            if (NLAMD > 0) then
                allocate(A_LAMD(NLAMD, NLAMD + 1))
                DLAMD(1:NLAMD) = 0.0

                do IP = 1, NTP
                    call PL_FLOW(ESTRS(:, IP), YIELD, DFDS(:, IP), 1)
                    ESTRS(:, IP) = ESTRS(:, IP) + YE(:, IP)
                end do

                PRER = 1.0e20
                do ITER = 1, iabs(MITER)
                    A_LAMD = 0.0

                    do IP = 1, NTP
                        IL = NNOD(IP)
                        if (IL == 0) cycle

                        call DF_MATRX(DFDS(:, IP), DF_IP, PSI, HEPSI)
                        CLAMD = 1.0 - HARDS * HEPSI / PSI

                        do I = 1, NSIG
                            A_LAMD(IL, NLAMD + 1) = A_LAMD(IL, NLAMD + 1) + &
                                                    RE(I) * DFDS(I, IP) * YE(I, IP) / PSI
                        end do

                        do JP = 1, NTP
                            JL = NNOD(JP)
                            if (JL == 0) cycle
                            call DF_MATRX(DFDS(:, JP), DF_JP, PSIB, HEPSI)
                            J0 = NSIG * (JP - 1)
                            do I = 1, NSIG
                                do J = 1, NSIG
                                    A_LAMD(IL, JL) = A_LAMD(IL, JL) - &
                                        (RE(I) * DFDS(I, IP) / PSI) * E_MAT(NSIG*(IP-1)+I, J0+J) * DF_JP(J)
                                end do
                            end do
                        end do

                        A_LAMD(IL, IL) = 1.0 - CLAMD + A_LAMD(IL, IL)
                    end do

                    ! Check convergence
                    RNORM = 0.0
                    do IP = 1, NLAMD
                        SS = A_LAMD(IP, NLAMD + 1)
                        do JP = 1, NLAMD
                            SS = SS - A_LAMD(IP, JP) * DLAMD(JP)
                        end do
                        A_LAMD(IP, NLAMD + 1) = SS
                        RNORM = RNORM + SS * SS
                    end do

                    RERR = sqrt(RNORM) / UNIAX
                    write(*,   95) RERR, NLAMD, ITER, INCRS
                    write(IPR, 95) RERR, NLAMD, ITER, INCRS
                    if (RERR < TOLER) exit
                    if (RERR > 10.0 * PRER .and. RERR > 1.0) then
                        write(*,   65)
                        write(IPR, 65)
                        exit
                    end if
                    PRER = RERR

                    ! Solve A_LAMD(1:NLAMD, 1:NLAMD) * DELTA_LAMD = A_LAMD(1:NLAMD, NLAMD+1)
                    call SLNPD(A_LAMD(1:NLAMD, 1:NLAMD), A_LAMD(1:NLAMD, NLAMD+1), &
                               D_DET, NLAMD, NLAMD, IPR)

                    ! Update plastic multipliers, boundary unknowns, internal displacements, and stresses
                    DLAMD(1:NLAMD) = DLAMD(1:NLAMD) + A_LAMD(1:NLAMD, NLAMD + 1)
                    call UPDATEV(NLAMD, A_LAMD(1:NLAMD, NLAMD + 1))

                    do IP = 1, NTP
                        call PL_FLOW(ESTRS(:, IP), YIELD, DFDS(:, IP), 1)
                    end do
                end do

                ! Update equivalent plastic strain, back-stress, and stress at converged step
                do IP = 1, NTP
                    if (NNOD(IP) == 0) cycle
                    if (DLAMD(NNOD(IP)) <= 0.0) cycle
                    EQUP2 = 0.0
                    do I = 1, NSIG
                        EQUP2 = EQUP2 + RE(I) * DFDS(I, IP) * DFDS(I, IP)
                    end do
                    DPLAS        = CEQP * sqrt(EQUP2) * DLAMD(NNOD(IP))
                    EPSTN(IP)    = EPSTN(IP) + DPLAS
                    CBACK        = DPLAS * HARDK / (UNIAX + HARDI * EPSTN(IP))
                    BACKS(:, IP) = BACKS(:, IP) + CBACK * ESTRS(:, IP)
                    ESTRS(:, IP) = ESTRS(:, IP) - CBACK * ESTRS(:, IP)
                end do

                deallocate(A_LAMD)
            end if

            FMULT = FMULT * FINCR
        end do

        deallocate(YE, DLAMD)
65      format(//, ' ** WARNING: CURRENT INCREMENT VARIATION IS LARGE **', /)
95      format(5X, 'RERR ==', 1PE14.6, ',  NYIELD ==', I4, ',  ITER ==', I3, ',  INCREM == ', I3)
    end subroutine P_M_ITER

    !-------------------------------------------------------------------
    !> @brief Scales trial elastic stress to yield surface and determines
    !>        excess stress load YE (SIG_SCALE).
    !-------------------------------------------------------------------
    subroutine SIG_SCALE(TOL, FMULT, YC_B, YC_U, YE_0, YE, NLAMD)
        implicit none
        real,    intent(in)  :: TOL, FMULT, YC_B(NN), YC_U(*), YE_0(NSIG, NTP)
        real,    intent(out) :: YE(NSIG, NTP)
        integer, intent(out) :: NLAMD

        integer :: IP
        real    :: ROVER, SIGY, YILD0, YIELD, SIGMT(NSIG), FGRAD(NSIG)

        NLAMD = 0
        X_B   = X_B + YC_B * FMULT
        if (L > 0) DSOL(1:2*L) = DSOL(1:2*L) + YC_U(1:2*L) * FMULT

        do IP = 1, NTP
            NNOD(IP) = 0
            ROVER    = 1.0
            SIGMT    = ESTRS(:, IP) + YE_0(:, IP) * FMULT
            SIGY     = UNIAX + HARDI * EPSTN(IP)

            call PL_FLOW(ESTRS(:, IP), YILD0, FGRAD, 0)
            call PL_FLOW(SIGMT,        YIELD, FGRAD, 0)

            if (YIELD <= SIGY) then
                ROVER = 0.0
            else
                if (YILD0 < SIGY - TOL) then
                    call SIGCROSS(ESTRS(:, IP), SIGMT, SIGY, YILD0, YIELD, ROVER)
                end if
                NLAMD    = NLAMD + 1
                NNOD(IP) = NLAMD
            end if

            ESTRS(:, IP) = ESTRS(:, IP) + YE_0(:, IP) * FMULT * (1.0 - ROVER)
            YE(:, IP)    = YE_0(:, IP)  * FMULT * ROVER
        end do
    end subroutine SIG_SCALE

    !-------------------------------------------------------------------
    !> @brief Finds yield surface crossing factor ALFA and overload
    !>        factor ROVER = 1 - ALFA (SIGCROSS).
    !-------------------------------------------------------------------
    subroutine SIGCROSS(SIGMN, SIGMT, SIGY, YILD0, YIELD_IN, ROVER)
        implicit none
        real, intent(in)  :: SIGMN(NSIG), SIGMT(NSIG), SIGY, YILD0, YIELD_IN
        real, intent(out) :: ROVER

        integer :: ITER, I
        real    :: TOLE, ALFAP, ALFAC, RESF, RESFC, YIELD, EQ
        real    :: SIGMA(NSIG), FGRAD(NSIG)

        TOLE  = 1.0e-5 * SIGY
        YIELD = YIELD_IN
        ALFAP = (SIGY - YILD0) / (YIELD - YILD0)
        ALFAC = ALFAP
        RESF  = abs(YIELD - SIGY)

        do ITER = 1, 10
            SIGMA = SIGMN + ALFAC * (SIGMT - SIGMN)
            call PL_FLOW(SIGMA, YIELD, FGRAD, 1)
            RESFC = YIELD - SIGY
            if (abs(RESFC) < TOLE) then
                ROVER = 1.0 - ALFAC
                return
            end if
            if (abs(RESFC) < RESF) then
                ALFAP = ALFAC
                RESF  = abs(RESFC)
            end if
            EQ = 0.0
            do I = 1, NSIG
                EQ = EQ + RE(I) * FGRAD(I) * (SIGMT(I) - SIGMN(I))
            end do
            if (abs(EQ) > 1.0e-30) ALFAC = ALFAC - RESFC / EQ
        end do

        ALFAC = ALFAP
        ROVER = 1.0 - ALFAC
    end subroutine SIGCROSS

    !-------------------------------------------------------------------
    !> @brief Evaluates equivalent stress YIELD and yield function
    !>        derivatives DFDS for 2D problems (PL_FLOW).
    !-------------------------------------------------------------------
    subroutine PL_FLOW(SIGMA, YIELD, DFDS_OUT, NFL)
        implicit none
        real,    intent(in)  :: SIGMA(NSIG)
        real,    intent(out) :: YIELD, DFDS_OUT(NSIG)
        integer, intent(in)  :: NFL

        integer :: I, J, K, IDX
        real    :: SMEAN, DEVIA(3, 3), VARJ2, VARJ3, SQUJ2, SINT3, THETA, SNPHI

        DEVIA = 0.0
        SMEAN = (SIGMA(1) + SIGMA(2)) / 3.0
        if (NSIG == 4) SMEAN = SMEAN + SIGMA(4) / 3.0

        IDX = 0
        do I = 1, 2
            do J = 1, 2
                IDX = IDX + 1
                DEVIA(I, J) = SIGMA(INY(IDX)) - SMEAN * DLT(I, J)
            end do
        end do

        VARJ2 = 0.0
        VARJ3 = 0.0
        do I = 1, 2
            do J = 1, 2
                VARJ2 = VARJ2 + 0.5 * DEVIA(I, J) * DEVIA(I, J)
                do K = 1, 2
                    VARJ3 = VARJ3 + DEVIA(I, J) * DEVIA(J, K) * DEVIA(K, I) / 3.0
                end do
            end do
        end do

        if (NSIG == 3) DEVIA(3, 3) = -SMEAN
        if (NSIG == 4) DEVIA(3, 3) = SIGMA(4) - SMEAN
        VARJ2 = VARJ2 + 0.5 * DEVIA(3, 3) * DEVIA(3, 3)
        VARJ3 = VARJ3 + (DEVIA(3, 3)**3) / 3.0

        SQUJ2 = sqrt(VARJ2)
        if (SQUJ2 < 1.0e-12 * UNIAX) then
            SINT3 = 0.0
        else
            SINT3 = -3.0 * ROOT3 * VARJ3 / (2.0 * VARJ2 * SQUJ2)
            if (abs(SINT3) > 1.0) SINT3 = sign(1.0, SINT3)
        end if
        THETA = asin(SINT3) / 3.0

        select case (NCRIT)
        case (1) ! Tresca
            YIELD = 2.0 * cos(THETA) * SQUJ2
        case (2) ! von Mises
            YIELD = ROOT3 * SQUJ2
        case (3) ! Mohr-Coulomb
            SNPHI = sin(FRICT)
            YIELD = (SMEAN * SNPHI + SQUJ2 * (cos(THETA) - sin(THETA) * SNPHI / ROOT3)) / CEQS
        case (4) ! Drucker-Prager
            YIELD = (ALFA * 3.0 * SMEAN + SQUJ2) / CEQS
        end select

        if (NFL == 0) return
        call DF_DSIG(SQUJ2, THETA, VARJ2, DEVIA, DFDS_OUT)
    end subroutine PL_FLOW

    !-------------------------------------------------------------------
    !> @brief Evaluates derivatives of the yield function df/d(sigma_ij)
    !>        for 2D problems (DF_DSIG).
    !-------------------------------------------------------------------
    subroutine DF_DSIG(SQUJ2, THETA, VARJ2, DEVIA, DFDS_OUT)
        implicit none
        real, intent(in)  :: SQUJ2, THETA, VARJ2, DEVIA(3, 3)
        real, intent(out) :: DFDS_OUT(NSIG)

        integer :: I, J, K, L_IDX, M
        real    :: CONS1, CONS2, CONS3, TANTH, TANT3, SINTH, COSTH, COST3, ABTHE
        real    :: VECA1(4), VECA2(4), VECA3(4)

        CONS1 = 0.0
        if (NCRIT == 3) CONS1 = sin(FRICT) / 3.0
        if (NCRIT == 4) CONS1 = ALFA

        VECA1(1:2) = 1.0
        VECA1(3)   = 0.0
        if (NSIG == 4) VECA1(4) = 1.0

        if (SQUJ2 < 1.0e-12 * UNIAX) then
            DFDS_OUT(1:NSIG) = (CONS1 / CEQS) * VECA1(1:NSIG)
            return
        end if

        TANTH = tan(THETA)
        TANT3 = tan(3.0 * THETA)
        SINTH = sin(THETA)
        COSTH = cos(THETA)
        COST3 = cos(3.0 * THETA)

        L_IDX = 0
        do I = 1, 2
            do J = I, 2
                L_IDX    = L_IDX + 1
                M        = INX(L_IDX)
                VECA2(M) = DEVIA(I, J) / (2.0 * SQUJ2)
                VECA3(M) = 0.0
                do K = 1, 2
                    VECA3(M) = VECA3(M) + DEVIA(I, K) * DEVIA(K, J) - VARJ2 * DLT(I, J) * (2.0 / 3.0)
                end do
            end do
        end do

        if (NSIG == 4) then
            VECA2(4) = DEVIA(3, 3) / (2.0 * SQUJ2)
            VECA3(4) = DEVIA(3, 3) * DEVIA(3, 3) - VARJ2 * (2.0 / 3.0)
        end if

        CONS2 = 0.0
        CONS3 = 0.0
        ABTHE = abs(THETA * (180.0 / PI))

        select case (NCRIT)
        case (1) ! Tresca
            if (ABTHE < 29.0) then
                CONS2 = 2.0 * (COSTH + SINTH * TANT3)
                CONS3 = ROOT3 * SINTH / (VARJ2 * COST3)
            else
                CONS2 = ROOT3
            end if
        case (2) ! von Mises
            CONS2 = ROOT3
        case (3) ! Mohr-Coulomb
            if (ABTHE < 29.0) then
                CONS2 = COSTH * (1.0 + TANTH * TANT3 + CONS1 * (TANT3 - TANTH) * ROOT3)
                CONS3 = (ROOT3 * SINTH + 3.0 * CONS1 * COSTH) / (2.0 * VARJ2 * COST3)
            else
                CONS2 = 0.5 * ROOT3 * (1.0 - sign(1.0, THETA) * CONS1)
            end if
        case (4) ! Drucker-Prager
            CONS2 = 1.0
        end select

        if (NCRIT > 2) then
            CONS1 = CONS1 / CEQS
            CONS2 = CONS2 / CEQS
            CONS3 = CONS3 / CEQS
        end if

        do I = 1, NSIG
            DFDS_OUT(I) = CONS1 * VECA1(I) + CONS2 * VECA2(I) + CONS3 * VECA3(I)
        end do
    end subroutine DF_DSIG

    !-------------------------------------------------------------------
    !> @brief Computes [D_e]{df/dsigma} stored in DF_OUT, PSI, and HEPSI
    !>        (DF_MATRX, Eqs. 10.6 & 7.68).
    !-------------------------------------------------------------------
    subroutine DF_MATRX(DFDS_IN, DF_OUT, PSI, HEPSI)
        implicit none
        real, intent(in)  :: DFDS_IN(NSIG)
        real, intent(out) :: DF_OUT(NSIG), PSI, HEPSI

        integer :: I, J, L_IDX
        real    :: AVKK, EQUP2

        AVKK = DFDS_IN(1) + DFDS_IN(2)
        if (NSIG == 4) AVKK = AVKK + DFDS_IN(4)

        L_IDX = 0
        do I = 1, 2
            do J = I, 2
                L_IDX = L_IDX + 1
                DF_OUT(INX(L_IDX)) = CLAME * DLT(I, J) * AVKK + DFDS_IN(INX(L_IDX))
            end do
        end do
        if (NSIG == 4) DF_OUT(4) = CLAME * AVKK + DFDS_IN(4)

        PSI   = 0.0
        EQUP2 = 0.0
        do I = 1, NSIG
            EQUP2 = EQUP2 + RE(I) * DFDS_IN(I) * DFDS_IN(I)
            PSI   = PSI   + RE(I) * DFDS_IN(I) * DF_OUT(I)
        end do
        HEPSI = CEQP * sqrt(EQUP2)
        PSI   = PSI + HARDS * HEPSI
    end subroutine DF_MATRX

    !-------------------------------------------------------------------
    !> @brief Updates boundary unknowns X_B, interior displacements DSOL,
    !>        and stresses ESTRS using plastic multiplier increments
    !>        DELTA_LAMD (UPDATEV).
    !-------------------------------------------------------------------
    subroutine UPDATEV(NLAMD, DELTA_LAMD)
        implicit none
        integer, intent(in) :: NLAMD
        real,    intent(in) :: DELTA_LAMD(NLAMD)

        integer :: I, J, IP, JP, J0
        real    :: PSI, HEPSI, DF_JP(NSIG), S_VEC(NSIG), SS

        ! 1. Pre-multiply [df] * DELTA_LAMD for all yielded nodes into a vector of length NTSF
        real :: D_SIGP(NTSF)
        D_SIGP = 0.0
        do JP = 1, NTP
            if (NNOD(JP) == 0) cycle
            call DF_MATRX(DFDS(:, JP), DF_JP, PSI, HEPSI)
            J0 = NSIG * (JP - 1)
            do J = 1, NSIG
                D_SIGP(J0 + J) = DF_JP(J) * DELTA_LAMD(NNOD(JP))
            end do
        end do

        ! 2. Update boundary unknowns X_B = X_B + AC_MAT * D_SIGP
        do I = 1, NN
            SS = X_B(I)
            do J = 1, NTSF
                SS = SS + AC_MAT(I, J) * D_SIGP(J)
            end do
            X_B(I) = SS
        end do

        ! 3. Update interior displacements DSOL = DSOL + EUI_MAT * D_SIGP
        if (L > 0) then
            do I = 1, 2 * L
                SS = DSOL(I)
                do J = 1, NTSF
                    SS = SS + EUI_MAT(I, J) * D_SIGP(J)
                end do
                DSOL(I) = SS
            end do
        end if

        ! 4. Update elastic trial stresses ESTRS = ESTRS + E_MAT * D_SIGP
        do IP = 1, NTP
            S_VEC = ESTRS(:, IP)
            do I = 1, NSIG
                do J = 1, NTSF
                    S_VEC(I) = S_VEC(I) + E_MAT(NSIG * (IP - 1) + I, J) * D_SIGP(J)
                end do
            end do
            ESTRS(:, IP) = S_VEC
        end do
    end subroutine UPDATEV

end module assembly_mod
