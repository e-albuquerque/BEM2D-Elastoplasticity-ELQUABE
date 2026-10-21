!-----------------------------------------------------------------------
! MODULE: assembly_mod
!> @brief Boundary element integration, matrix assembly, and stress
!>        recovery routines (GHMATEQ, EXTINEQ, LOCINEQ, INTEREQ,
!>        BOUND_STRESS, CALC_DFFORMA, SIGMAEQ) using element connectivity
!>        CONEC and parallelized with OpenMP.
!-----------------------------------------------------------------------
module assembly_mod
    use data_mod
    implicit none

contains

    !-------------------------------------------------------------------
    !> @brief Computes the G and H matrices in parallel using OpenMP
    !>        and element connectivity CONEC, then forms the linear
    !>        system A * X = F (where A is stored in H and F in FI).
    !-------------------------------------------------------------------
    subroutine GHMATEQ()
        implicit none
        real :: HW(2, 6), GW(2, 6), CH
        integer :: EL, I, J, K, LL, NODO
        integer :: NO1, NO2, NO3, NODE_LOC, DIR, NO_GLOB, COL_H, COL_G
        integer :: FIRST_G_DOF(NN)

        H  = 0.0
        G  = 0.0
        FI = 0.0

        ! Compute the GW and HW matrices for each collocation point
        ! and each boundary element in parallel across CPU threads
        !$omp parallel do default(shared) private(LL, EL, NO1, NO2, NO3, NODO, HW, GW, K, NODE_LOC, NO_GLOB, DIR, J, COL_H, COL_G) schedule(dynamic)
        do LL = 1, N
            do EL = 1, NE
                NO1 = CONEC(EL, 1)
                NO2 = CONEC(EL, 2)
                NO3 = CONEC(EL, 3)

                ! Check if collocation node LL belongs to element EL
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
                    call EXTINEQ(X(LL), Y(LL), X(NO1), Y(NO1), X(NO2), Y(NO2), X(NO3), Y(NO3), HW, GW)
                    call LOCINEQ(X(NO1), Y(NO1), X(NO2), Y(NO2), X(NO3), Y(NO3), GW, NODO)
                else
                    call EXTINEQ(X(LL), Y(LL), X(NO1), Y(NO1), X(NO2), Y(NO2), X(NO3), Y(NO3), HW, GW)
                end if

                ! Plug the GW and HW matrices into the general G and H matrices
                ! using element connectivity CONEC(EL, 1:3)
                do K = 1, 2
                    do NODE_LOC = 1, 3
                        NO_GLOB = CONEC(EL, NODE_LOC)
                        do DIR = 1, 2
                            J     = 2 * NODE_LOC - 2 + DIR
                            COL_G = 6 * EL - 6 + J
                            COL_H = 2 * NO_GLOB - 2 + DIR

                            G(2*LL-2+K, COL_G) = G(2*LL-2+K, COL_G) + GW(K, J)
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

            ! Add one to the diagonal coefficients for external problems
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
                            ! First element prescribing displacement at this global DOF:
                            ! swap column COL_H of H with column COL_G of -G * GE
                            FIRST_G_DOF(COL_H) = COL_G
                            do K = 1, NN
                                CH = H(K, COL_H)
                                H(K, COL_H) = -G(K, COL_G) * GE
                                G(K, COL_G) = -CH
                            end do
                        else
                            ! Second element also prescribing displacement at this shared node:
                            ! unknown traction is continuous across the node
                            do K = 1, NN
                                H(K, COL_H) = H(K, COL_H) - G(K, COL_G) * GE
                                G(K, COL_G) = 0.0
                            end do
                        end if
                    end if
                end do
            end do
        end do

        ! Form the right-hand side vector F stored in FI
        do I = 1, NN
            FI(I) = 0.0
            do J = 1, 6 * NE
                FI(I) = FI(I) + G(I, J) * DFI(J)
            end do
        end do

    end subroutine GHMATEQ

    !-------------------------------------------------------------------
    !> @brief Computes the HW and GW matrices that relate a collocation
    !>        node (XP, YP) with a boundary element using 10-point
    !>        Gauss quadrature (non-singular integration).
    !-------------------------------------------------------------------
    subroutine EXTINEQ(XP, YP, X1, Y1, X2, Y2, X3, Y3, HW, GW)
        implicit none
        real, intent(in)  :: XP, YP, X1, Y1, X2, Y2, X3, Y3
        real, intent(out) :: HW(2, 6), GW(2, 6)

        real, parameter :: GI(10) = [ &
            0.9739065285, -0.9739065285,  0.8650633666, -0.8650633666, &
            0.6794095683, -0.6794095682,  0.4333953941, -0.4333953941, &
            0.1488743389, -0.1488743389 ]
        real, parameter :: OME(10) = [ &
            0.0666713443, 0.0666713443,  0.1494513491, 0.1494513491, &
            0.2190863625, 0.2190863625,  0.2692667193, 0.2692667193, &
            0.2955242247, 0.2955242247 ]

        real :: A, B, C, D, DE
        real :: F1, F2, F3, XCO, YCO, XJA, ETA1, ETA2, RA, RD1, RD2, RDN
        integer :: I

        HW = 0.0
        GW = 0.0

        A  = X3 - 2.0 * X2 + X1
        B  = (X3 - X1) / 2.0
        C  = Y3 - 2.0 * Y2 + Y1
        D  = (Y3 - Y1) / 2.0
        DE = 4.0 * 3.141592 * (1.0 - XNU)

        do I = 1, 10
            ! Shape functions at the integration point
            F1 = GI(I) * (GI(I) - 1.0) * 0.5
            F2 = 1.0 - GI(I)**2
            F3 = GI(I) * (GI(I) + 1.0) * 0.5

            ! Geometrical properties at the integration point
            XCO  = X1 * F1 + X2 * F2 + X3 * F3
            YCO  = Y1 * F1 + Y2 * F2 + Y3 * F3
            XJA  = sqrt((GI(I) * A + B)**2 + (GI(I) * C + D)**2)
            ETA1 =  (GI(I) * C + D) / XJA
            ETA2 = -(GI(I) * A + B) / XJA
            RA   = sqrt((XP - XCO)**2 + (YP - YCO)**2)
            RD1  = (XCO - XP) / RA
            RD2  = (YCO - YP) / RA
            RDN  = RD1 * ETA1 + RD2 * ETA2

            ! Compute GW and HW matrices
            GW(1, 1) = GW(1, 1) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD1**2) * OME(I) * XJA * F1 / (2.0 * DE * GE)
            GW(1, 2) = GW(1, 2) + RD1 * RD2 * OME(I) * XJA * F1 / (2.0 * DE * GE)
            GW(2, 1) = GW(1, 2)
            GW(2, 2) = GW(2, 2) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD2**2) * OME(I) * XJA * F1 / (2.0 * DE * GE)
            HW(1, 1) = HW(1, 1) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD1**2) / (RA * DE) * OME(I) * XJA * F1
            HW(1, 2) = HW(1, 2) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA1 * RD2 - ETA2 * RD1)) * &
                                  OME(I) * XJA * F1 / (RA * DE)
            HW(2, 1) = HW(2, 1) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA2 * RD1 - ETA1 * RD2)) * &
                                  OME(I) * XJA * F1 / (RA * DE)
            HW(2, 2) = HW(2, 2) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD2**2) * OME(I) * XJA * F1 / (RA * DE)

            GW(1, 3) = GW(1, 3) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD1**2) * OME(I) * XJA * F2 / (2.0 * DE * GE)
            GW(1, 4) = GW(1, 4) + RD1 * RD2 * OME(I) * XJA * F2 / (2.0 * DE * GE)
            GW(2, 3) = GW(1, 4)
            GW(2, 4) = GW(2, 4) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD2**2) * OME(I) * XJA * F2 / (2.0 * DE * GE)
            HW(1, 3) = HW(1, 3) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD1**2) / (RA * DE) * OME(I) * XJA * F2
            HW(1, 4) = HW(1, 4) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA1 * RD2 - ETA2 * RD1)) * &
                                  OME(I) * XJA * F2 / (RA * DE)
            HW(2, 3) = HW(2, 3) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA2 * RD1 - ETA1 * RD2)) * &
                                  OME(I) * XJA * F2 / (RA * DE)
            HW(2, 4) = HW(2, 4) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD2**2) * OME(I) * XJA * F2 / (RA * DE)

            GW(1, 5) = GW(1, 5) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD1**2) * OME(I) * XJA * F3 / (2.0 * DE * GE)
            GW(1, 6) = GW(1, 6) + RD1 * RD2 * OME(I) * XJA * F3 / (2.0 * DE * GE)
            GW(2, 5) = GW(1, 6)
            GW(2, 6) = GW(2, 6) + ((3.0 - 4.0 * XNU) * log(1.0 / RA) + RD2**2) * OME(I) * XJA * F3 / (2.0 * DE * GE)
            HW(1, 5) = HW(1, 5) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD1**2) / (RA * DE) * OME(I) * XJA * F3
            HW(1, 6) = HW(1, 6) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA1 * RD2 - ETA2 * RD1)) * &
                                  OME(I) * XJA * F3 / (RA * DE)
            HW(2, 5) = HW(2, 5) - (RDN * 2.0 * RD1 * RD2 + (1.0 - 2.0 * XNU) * (ETA2 * RD1 - ETA1 * RD2)) * &
                                  OME(I) * XJA * F3 / (RA * DE)
            HW(2, 6) = HW(2, 6) - RDN * ((1.0 - 2.0 * XNU) + 2.0 * RD2**2) * OME(I) * XJA * F3 / (RA * DE)
        end do

    end subroutine EXTINEQ

    !-------------------------------------------------------------------
    !> @brief Computes the GW matrix when the collocation point is one
    !>        of the nodes of the integration element (singular/logarithmic
    !>        quadrature).
    !-------------------------------------------------------------------
    subroutine LOCINEQ(XG1, YG1, XG2, YG2, XG3, YG3, GW, NODO)
        implicit none
        real, intent(in)    :: XG1, YG1, XG2, YG2, XG3, YG3
        integer, intent(in) :: NODO
        real, intent(out)   :: GW(2, 6)

        ! Data for standard Gauss quadrature
        real, parameter :: GI(10) = [ &
            0.9739065285, -0.9739065285,  0.8650633666, -0.8650633666, &
            0.6794095682, -0.6794095682,  0.4333953941, -0.4333953941, &
            0.1488743389, -0.1488743389 ]
        real, parameter :: OME(10) = [ &
            0.0666713443, 0.0666713443,  0.1494513491, 0.1494513491, &
            0.2190863625, 0.2190863625,  0.2692667193, 0.2692667193, &
            0.2955242247, 0.2955242247 ]

        ! Data for logarithmic quadrature
        real, parameter :: GIL(10) = [ &
            0.0090426309, 0.0539712662, 0.1353118246, 0.2470524162, &
            0.3802125396, 0.5237923179, 0.6657752055, 0.7941904160, &
            0.8981610912, 0.9688479887 ]
        real, parameter :: OMEL(10) = [ &
            0.1209551319, 0.1863635425, 0.1956608732, 0.1735771421, &
            0.1356956729, 0.0936467585, 0.0557877273, 0.0271598109, &
            0.0095151826, 0.0016381576 ]

        real :: X1, Y1, X2, Y2, X3, Y3
        real :: A1, B1, A2, B2, A, B, C, CONT1, CONT2
        real :: T1, T2, T3, F1, F2, F3, FL1, FL2, FL3, FLN1, FLN2, FLN3
        real :: XJA1, XJA11, XJA2, XLO, S1, S2, S3
        integer :: I

        X1 = 0.0; Y1 = 0.0; X2 = 0.0; Y2 = 0.0; X3 = 0.0; Y3 = 0.0
        A1 = 0.0; B1 = 0.0; A2 = 0.0; B2 = 0.0
        XJA1 = 0.0; XJA11 = 0.0; XJA2 = 0.0; XLO = 0.0
        S1 = 0.0; S2 = 0.0; S3 = 0.0

        ! Set local coordinate system
        select case (NODO)
        case (1)
            X3 = XG3 - XG1
            Y3 = YG3 - YG1
            X2 = XG2 - XG1
            Y2 = YG2 - YG1
            A1 = (X3 - 2.0 * X2) * 0.5
            B1 = X2
            A2 = (Y3 - 2.0 * Y2) * 0.5
            B2 = Y2
        case (2)
            X3 = XG3 - XG2
            Y3 = YG3 - YG2
            X1 = XG1 - XG2
            Y1 = YG1 - YG2
            A1 = X1 + X3
            B1 = X3 - X1
            A2 = Y1 + Y3
            B2 = Y3 - Y1
        case (3)
            X2 = XG2 - XG3
            Y2 = YG2 - YG3
            X1 = XG1 - XG3
            Y1 = YG1 - YG3
            A1 = (X1 - 2.0 * X2) * 0.5
            B1 = -X2
            A2 = (Y1 - 2.0 * Y2) * 0.5
            B2 = -Y2
        end select

        GW = 0.0
        A = A1**2 + A2**2
        B = 2.0 * (A1 * B1 + A2 * B2)
        C = B1**2 + B2**2
        CONT1 = (3.0 - 4.0 * XNU) / (8.0 * 3.1415926 * GE * (1.0 - XNU))
        CONT2 = CONT1 / (3.0 - 4.0 * XNU)

        do I = 1, 10
            T1 = ((A1 * GI(I) + B1)**2) / ((A2 * GI(I) + B2)**2 + (A1 * GI(I) + B1)**2)
            T2 = ((A2 * GI(I) + B2)**2) / ((A2 * GI(I) + B2)**2 + (A1 * GI(I) + B1)**2)
            T3 = (A1 * GI(I) + B1) * (A2 * GI(I) + B2) / ((A1 * GI(I) + B1)**2 + (A2 * GI(I) + B2)**2)

            ! Compute shape functions for numerical integration
            F3   = 0.5 * GI(I) * (GI(I) + 1.0)
            F2   = 1.0 - GI(I)**2
            F1   = 0.5 * GI(I) * (GI(I) - 1.0)
            FL3  = GIL(I) * (2.0 * GIL(I) - 1.0)
            FL2  = 4.0 * GIL(I) * (1.0 - GIL(I))
            FL1  = (GIL(I) - 1.0) * (2.0 * GIL(I) - 1.0)
            FLN3 = 0.5 * GIL(I) * (GIL(I) + 1.0)
            FLN2 = 1.0 - GIL(I)**2
            FLN1 = 0.5 * GIL(I) * (GIL(I) - 1.0)

            ! Compute GW coefficients
            select case (NODO)
            case (1)
                XJA1 = sqrt((4.0 * A1 * GIL(I) - 2.0 * A1 + 0.5 * X3)**2 + &
                            (4.0 * A2 * GIL(I) - 2.0 * A2 + 0.5 * Y3)**2) * 2.0
                XJA2 = sqrt((A1 * GI(I) * 2.0 + 0.5 * X3)**2 + (A2 * GI(I) * 2.0 + 0.5 * Y3)**2)
                XLO  = -log(2.0 * sqrt((GI(I) * A1 + B1)**2 + (GI(I) * A2 + B2)**2))
                S3   = CONT1 * (FL3 * XJA1 * OMEL(I) + F3 * XJA2 * XLO * OME(I))
                S2   = CONT1 * (FL2 * XJA1 * OMEL(I) + F2 * XJA2 * XLO * OME(I))
                S1   = CONT1 * (FL1 * XJA1 * OMEL(I) + F1 * XJA2 * XLO * OME(I))
            case (2)
                XJA1  = sqrt((0.5 * B1 - A1 * GIL(I))**2 + (0.5 * B2 - A2 * GIL(I))**2)
                XJA11 = sqrt((0.5 * B1 + A1 * GIL(I))**2 + (0.5 * B2 + A2 * GIL(I))**2)
                XJA2  = sqrt((0.5 * B1 + A1 * GI(I))**2  + (0.5 * B2 + A2 * GI(I))**2)
                XLO   = -0.5 * log((GI(I) * A1 * 0.5 + B1 * 0.5)**2 + (GI(I) * A2 * 0.5 + B2 * 0.5)**2)
                S3    = CONT1 * ((FLN1 * XJA1 + FLN3 * XJA11) * OMEL(I) + F3 * XJA2 * XLO * OME(I))
                S2    = CONT1 * (FLN2 * (XJA1 + XJA11) * OMEL(I) + F2 * XJA2 * XLO * OME(I))
                S1    = CONT1 * ((FLN3 * XJA1 + FLN1 * XJA11) * OMEL(I) + F1 * XJA2 * XLO * OME(I))
            case (3)
                XJA1 = sqrt((2.0 * A1 - 4.0 * A1 * GIL(I) - 0.5 * X1)**2 + &
                            (2.0 * A2 - 4.0 * A2 * GIL(I) - 0.5 * Y1)**2) * 2.0
                XJA2 = sqrt((2.0 * A1 * GI(I) - 0.5 * X1)**2 + (2.0 * A2 * GI(I) - 0.5 * Y1)**2)
                XLO  = -log(2.0 * sqrt((A1 * GI(I) + B1)**2 + (A2 * GI(I) + B2)**2))
                S3   = CONT1 * (FL1 * XJA1 * OMEL(I) + F3 * XJA2 * XLO * OME(I))
                S2   = CONT1 * (FL2 * XJA1 * OMEL(I) + F2 * XJA2 * XLO * OME(I))
                S1   = CONT1 * (FL3 * XJA1 * OMEL(I) + F1 * XJA2 * XLO * OME(I))
            end select

            GW(1, 5) = GW(1, 5) + S3 + CONT2 * F3 * T1 * XJA2 * OME(I)
            GW(1, 3) = GW(1, 3) + S2 + CONT2 * F2 * T1 * XJA2 * OME(I)
            GW(1, 1) = GW(1, 1) + S1 + CONT2 * F1 * T1 * XJA2 * OME(I)
            GW(1, 6) = GW(1, 6) + CONT2 * T3 * XJA2 * F3 * OME(I)
            GW(1, 4) = GW(1, 4) + CONT2 * T3 * XJA2 * F2 * OME(I)
            GW(1, 2) = GW(1, 2) + CONT2 * T3 * XJA2 * F1 * OME(I)
            GW(2, 5) = GW(1, 6)
            GW(2, 3) = GW(1, 4)
            GW(2, 1) = GW(1, 2)
            GW(2, 6) = GW(2, 6) + S3 + CONT2 * T2 * XJA2 * F3 * OME(I)
            GW(2, 4) = GW(2, 4) + S2 + CONT2 * T2 * XJA2 * F2 * OME(I)
            GW(2, 2) = GW(2, 2) + S1 + CONT2 * T2 * XJA2 * F1 * OME(I)
        end do

    end subroutine LOCINEQ

    !-------------------------------------------------------------------
    !> @brief Rearranges FI and DFI arrays using CONEC to store boundary
    !>        displacements and tractions, computes boundary stresses via
    !>        shape function derivatives, and evaluates displacements and
    !>        stresses at internal points in parallel using OpenMP.
    !-------------------------------------------------------------------
    subroutine INTEREQ()
        implicit none
        real :: HW(2, 6), GW(2, 6)
        real :: D11(6), D12(6), D22(6), S11(6), S12(6), S22(6)
        real :: CH, SX, TXY, SY
        integer :: EL, J, K
        integer :: NO1, NO2, NO3, NODE_LOC, DIR, NO_GLOB, COL_H, COL_G
        integer :: FIRST_G_DOF(NN)

        ! Rearrange FI and DFI arrays using element connectivity CONEC
        ! so FI stores all boundary displacements and DFI stores all tractions
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
                            CH         = FI(COL_H) * GE
                            FI(COL_H)  = DFI(COL_G)
                            DFI(COL_G) = CH
                        else
                            DFI(COL_G) = DFI(FIRST_G_DOF(COL_H))
                        end if
                    end if
                end do
            end do
        end do

        ! Compute boundary nodal stresses using shape function derivatives and CONEC
        call BOUND_STRESS()

        ! Compute the values of stresses and displacements at internal points
        if (L == 0) return

        !$omp parallel do default(shared) private(K, EL, NO1, NO2, NO3, NODE_LOC, NO_GLOB, DIR, J, COL_G, COL_H, HW, GW, D11, D12, D22, S11, S12, S22, SX, TXY, SY) schedule(dynamic)
        do K = 1, L
            SSOL(3*K-2) = 0.0
            SSOL(3*K-1) = 0.0
            SSOL(3*K)   = 0.0
            DSOL(2*K-1) = 0.0
            DSOL(2*K)   = 0.0
            do EL = 1, NE
                NO1 = CONEC(EL, 1)
                NO2 = CONEC(EL, 2)
                NO3 = CONEC(EL, 3)

                call EXTINEQ(CX(K), CY(K), X(NO1), Y(NO1), X(NO2), Y(NO2), &
                             X(NO3), Y(NO3), HW, GW)
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
                    end do
                end do
            end do

            SX  = SSOL(3*K-2)
            TXY = SSOL(3*K-1)
            SY  = SSOL(3*K)
            SVM_INT(K) = sqrt(SX**2 - SX * SY + SY**2 + 3.0 * TXY**2)
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
    !>        sigma_y) and von Mises equivalent stress at each boundary
    !>        node using element connectivity CONEC and shape function
    !>        derivatives (Section 3.6 & calcula_tens_cont.m).
    !-------------------------------------------------------------------
    subroutine BOUND_STRESS()
        implicit none
        integer :: EL, NO, NO1, NO2, NO3, NO_GLOBAL, I
        integer :: NODE_COUNT(N)
        real :: X1, Y1, X2, Y2, X3, Y3
        real :: QSI, DN1DQSI, DN2DQSI, DN3DQSI
        real :: DXDQSI, DYDQSI, DGAMADQSI
        real :: SX_TAN, SY_TAN, NX_VEC, NY_VEC
        real :: LIJ(2, 2), LJI(2, 2), S_LOC(2, 2), SIGMA_GLOB(2, 2)
        real :: TX, TY, S11, S12, S22, E22, DVDQSI
        real :: UNO1(2), UNO2(2), UNO3(2), V1(2), V2(2), V3(2)
        real :: E_MOD, SX, TXY, SY

        SBOU       = 0.0
        NODE_COUNT = 0
        E_MOD      = 2.0 * GE * (1.0 + XNU)

        do EL = 1, NE
            ! Global node indices of the 3-node quadratic element from CONEC
            NO1 = CONEC(EL, 1)
            NO2 = CONEC(EL, 2)
            NO3 = CONEC(EL, 3)

            ! Nodal coordinates
            X1 = X(NO1); Y1 = Y(NO1)
            X2 = X(NO2); Y2 = Y(NO2)
            X3 = X(NO3); Y3 = Y(NO3)

            ! Nodal displacements in global coordinates
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

                ! Transformation matrix (global -> local normal-tangential)
                LIJ(1, 1) =  NX_VEC; LIJ(1, 2) = NY_VEC
                LIJ(2, 1) = -NY_VEC; LIJ(2, 2) = NX_VEC

                ! Global node index for current local node
                NO_GLOBAL = CONEC(EL, NO)

                ! Element traction at current local node in global coordinates
                TX = DFI(6*EL - 6 + 2*NO - 1)
                TY = DFI(6*EL - 6 + 2*NO)

                ! Traction in local coordinates (s11 = normal stress, s12 = shear stress)
                S11 = LIJ(1, 1) * TX + LIJ(1, 2) * TY
                S12 = LIJ(2, 1) * TX + LIJ(2, 2) * TY

                ! Nodal displacements in local coordinates
                V1 = matmul(LIJ, UNO1)
                V2 = matmul(LIJ, UNO2)
                V3 = matmul(LIJ, UNO3)

                ! Tangential displacement derivative and membrane stress s22
                DVDQSI = DN1DQSI * V1(2) + DN2DQSI * V2(2) + DN3DQSI * V3(2)
                E22    = DVDQSI / DGAMADQSI
                S22    = (E_MOD * E22) / (1.0 - XNU**2) + (XNU * (1.0 + XNU) / (1.0 - XNU**2)) * S11

                ! Transformation matrix (local -> global)
                LJI(1, 1) = NX_VEC; LJI(1, 2) = -NY_VEC
                LJI(2, 1) = NY_VEC; LJI(2, 2) =  NX_VEC

                S_LOC(1, 1) = S11; S_LOC(1, 2) = S12
                S_LOC(2, 1) = S12; S_LOC(2, 2) = S22

                ! Global stress tensor sigma = LJI * S_LOC * LJI^T
                SIGMA_GLOB = matmul(LJI, matmul(S_LOC, transpose(LJI)))

                ! Accumulate nodal stress components [sigma_x, tau_xy, sigma_y]
                SBOU(NO_GLOBAL, 1) = SBOU(NO_GLOBAL, 1) + SIGMA_GLOB(1, 1)
                SBOU(NO_GLOBAL, 2) = SBOU(NO_GLOBAL, 2) + SIGMA_GLOB(1, 2)
                SBOU(NO_GLOBAL, 3) = SBOU(NO_GLOBAL, 3) + SIGMA_GLOB(2, 2)
                NODE_COUNT(NO_GLOBAL) = NODE_COUNT(NO_GLOBAL) + 1
            end do
        end do

        ! Average at shared extreme nodes and compute von Mises equivalent stress
        do I = 1, N
            if (NODE_COUNT(I) > 0) then
                SBOU(I, 1:3) = SBOU(I, 1:3) / real(NODE_COUNT(I))
            end if
            SX  = SBOU(I, 1)
            TXY = SBOU(I, 2)
            SY  = SBOU(I, 3)
            SBOU(I, 4) = sqrt(SX**2 - SX * SY + SY**2 + 3.0 * TXY**2)
        end do

    end subroutine BOUND_STRESS

    !-------------------------------------------------------------------
    !> @brief Computes the S and D stress kernel matrices using Gauss
    !>        quadrature to evaluate stress at an internal point (XP, YP).
    !-------------------------------------------------------------------
    subroutine SIGMAEQ(XP, YP, X1, Y1, X2, Y2, X3, Y3, D11, D12, D22, S11, S12, S22)
        implicit none
        real, intent(in)  :: XP, YP, X1, Y1, X2, Y2, X3, Y3
        real, intent(out) :: D11(6), D12(6), D22(6), S11(6), S12(6), S22(6)

        real, parameter :: GI(10) = [ &
            0.9739065285, -0.9739065285,  0.8650633666, -0.8650633666, &
            0.6794095683, -0.6794095682,  0.4333953941, -0.4333953941, &
            0.1488743389, -0.1488743389 ]
        real, parameter :: OME(10) = [ &
            0.0666713443, 0.0666713443,  0.1494513491, 0.1494513491, &
            0.2190863625, 0.2190863625,  0.2692667193, 0.2692667193, &
            0.2955242247, 0.2955242247 ]

        real :: FA, AL, A, B, C, D, DE
        real :: F1, F2, F3, XCO, YCO, XJA, ETA1, ETA2, RA, RD1, RD2, RDN
        integer :: I

        D11 = 0.0
        D12 = 0.0
        D22 = 0.0
        S11 = 0.0
        S12 = 0.0
        S22 = 0.0

        FA = 1.0 - 4.0 * XNU
        AL = 1.0 - 2.0 * XNU
        A  = X3 - 2.0 * X2 + X1
        B  = (X3 - X1) / 2.0
        C  = Y3 - 2.0 * Y2 + Y1
        D  = (Y3 - Y1) / 2.0
        DE = 4.0 * 3.141592 * (1.0 - XNU)

        do I = 1, 10
            ! Shape functions at the integration point
            F1 = GI(I) * (GI(I) - 1.0) * 0.5
            F2 = 1.0 - GI(I)**2
            F3 = GI(I) * (GI(I) + 1.0) * 0.5

            ! Geometrical parameters
            XCO  = X1 * F1 + X2 * F2 + X3 * F3
            YCO  = Y1 * F1 + Y2 * F2 + Y3 * F3
            XJA  = sqrt((GI(I) * A + B)**2 + (GI(I) * C + D)**2)
            ETA1 =  (GI(I) * C + D) / XJA
            ETA2 = -(GI(I) * A + B) / XJA
            RA   = sqrt((XP - XCO)**2 + (YP - YCO)**2)
            RD1  = (XCO - XP) / RA
            RD2  = (YCO - YP) / RA
            RDN  = RD1 * ETA1 + RD2 * ETA2

            ! Compute D and S coefficients
            D11(1) = D11(1) + (AL * RD1 + 2.0 * RD1**3) * OME(I) * XJA * F1 / (DE * RA)
            D11(2) = D11(2) + (2.0 * RD1**2 * RD2 - AL * RD2) * OME(I) * XJA * F1 / (DE * RA)
            D11(3) = D11(3) + (AL * RD1 + 2.0 * RD1**3) * OME(I) * XJA * F2 / (DE * RA)
            D11(4) = D11(4) + (2.0 * RD1**2 * RD2 - AL * RD2) * OME(I) * XJA * F2 / (DE * RA)
            D11(5) = D11(5) + (AL * RD1 + 2.0 * RD1**3) * OME(I) * XJA * F3 / (DE * RA)
            D11(6) = D11(6) + (2.0 * RD1**2 * RD2 - AL * RD2) * OME(I) * XJA * F3 / (DE * RA)

            D12(1) = D12(1) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F1 / (DE * RA) * OME(I) * XJA
            D12(2) = D12(2) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F1 / (DE * RA) * OME(I) * XJA
            D12(3) = D12(3) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F2 / (DE * RA) * OME(I) * XJA
            D12(4) = D12(4) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F2 / (DE * RA) * OME(I) * XJA
            D12(5) = D12(5) + (AL * RD2 + 2.0 * RD1**2 * RD2) * F3 / (DE * RA) * OME(I) * XJA
            D12(6) = D12(6) + (AL * RD1 + 2.0 * RD1 * RD2**2) * F3 / (DE * RA) * OME(I) * XJA

            D22(1) = D22(1) + (2.0 * RD1 * RD2**2 - AL * RD1) * F1 / (DE * RA) * OME(I) * XJA
            D22(2) = D22(2) + (AL * RD2 + 2.0 * RD2**3) * F1 / (DE * RA) * OME(I) * XJA
            D22(3) = D22(3) + (2.0 * RD1 * RD2**2 - AL * RD1) * F2 / (DE * RA) * OME(I) * XJA
            D22(4) = D22(4) + (AL * RD2 + 2.0 * RD2**3) * F2 / (DE * RA) * OME(I) * XJA
            D22(5) = D22(5) + (2.0 * RD1 * RD2**2 - AL * RD1) * F3 / (DE * RA) * OME(I) * XJA
            D22(6) = D22(6) + (AL * RD2 + 2.0 * RD2**3) * F3 / (DE * RA) * OME(I) * XJA

            S11(1) = S11(1) + (2.0 * RDN * (AL * RD1 + XNU * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * XNU * ETA1 * &
                     RD1**2 + AL * (2.0 * ETA1 * RD1**2 + 2.0 * ETA1) - FA * ETA1) * 2.0 * GE * F1 / (DE * RA**2) * &
                     OME(I) * XJA
            S11(2) = S11(2) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * XNU * ETA1 * RD1 * RD2 + &
                     AL * 2.0 * ETA2 * RD1**2 - FA * ETA2) * 2.0 * GE * F1 / (DE * RA**2) * OME(I) * XJA
            S11(3) = S11(3) + (2.0 * RDN * (AL * RD1 + XNU * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * XNU * ETA1 * &
                     RD1**2 + AL * (2.0 * ETA1 * RD1**2 + 2.0 * ETA1) - FA * ETA1) * 2.0 * GE * F2 / (DE * RA**2) * &
                     OME(I) * XJA
            S11(4) = S11(4) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * XNU * ETA1 * RD1 * RD2 + &
                     AL * 2.0 * ETA2 * RD1**2 - FA * ETA2) * 2.0 * GE * F2 / (DE * RA**2) * OME(I) * XJA
            S11(5) = S11(5) + (2.0 * RDN * (AL * RD1 + XNU * 2.0 * RD1 - 4.0 * RD1**3) + 4.0 * XNU * ETA1 * &
                     RD1**2 + AL * (2.0 * ETA1 * RD1**2 + 2.0 * ETA1) - FA * ETA1) * 2.0 * GE * F3 / (DE * RA**2) * &
                     OME(I) * XJA
            S11(6) = S11(6) + (2.0 * RDN * (AL * RD2 - 4.0 * RD1**2 * RD2) + 4.0 * XNU * ETA1 * RD1 * RD2 + &
                     AL * 2.0 * ETA2 * RD1**2 - FA * ETA2) * 2.0 * GE * F3 / (DE * RA**2) * OME(I) * XJA

            S12(1) = S12(1) + (2.0 * RDN * (XNU * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * XNU * (ETA1 * RD2 * &
                     RD1 + ETA2 * RD1**2) + AL * (2.0 * ETA1 * RD1 * RD2 + ETA2)) * 2.0 * GE * F1 / &
                     (DE * RA**2) * OME(I) * XJA
            S12(2) = S12(2) + (2.0 * RDN * (XNU * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * XNU * (ETA1 * RD2**2 + &
                     ETA2 * RD1 * RD2) + AL * (2.0 * ETA2 * RD1 * RD2 + ETA1)) * 2.0 * GE * F1 / &
                     (DE * RA**2) * OME(I) * XJA
            S12(3) = S12(3) + (2.0 * RDN * (XNU * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * XNU * (ETA1 * RD2 * &
                     RD1 + ETA2 * RD1**2) + AL * (2.0 * ETA1 * RD1 * RD2 + ETA2)) * 2.0 * GE * F2 / &
                     (DE * RA**2) * OME(I) * XJA
            S12(4) = S12(4) + (2.0 * RDN * (XNU * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * XNU * (ETA1 * RD2**2 + &
                     ETA2 * RD1 * RD2) + AL * (2.0 * ETA2 * RD1 * RD2 + ETA1)) * 2.0 * GE * F2 / &
                     (DE * RA**2) * OME(I) * XJA
            S12(5) = S12(5) + (2.0 * RDN * (XNU * RD2 - 4.0 * RD1**2 * RD2) + 2.0 * XNU * (ETA1 * RD2 * &
                     RD1 + ETA2 * RD1**2) + AL * (2.0 * ETA1 * RD1 * RD2 + ETA2)) * 2.0 * GE * F3 / &
                     (DE * RA**2) * OME(I) * XJA
            S12(6) = S12(6) + (2.0 * RDN * (XNU * RD1 - 4.0 * RD1 * RD2**2) + 2.0 * XNU * (ETA1 * RD2**2 + &
                     ETA2 * RD1 * RD2) + AL * (2.0 * ETA2 * RD1 * RD2 + ETA1)) * 2.0 * GE * F3 / &
                     (DE * RA**2) * OME(I) * XJA

            S22(1) = S22(1) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * XNU * ETA2 * RD1 * RD2 + &
                     AL * 2.0 * ETA1 * RD2**2 - FA * ETA1) * 2.0 * GE * F1 / (DE * RA**2) * OME(I) * XJA
            S22(2) = S22(2) + (2.0 * RDN * (AL * RD2 + 2.0 * XNU * RD2 - 4.0 * RD2**3) + 4.0 * XNU * ETA2 * &
                     RD2**2 + AL * (2.0 * ETA2 * RD2**2 + 2.0 * ETA2) - FA * ETA2) * 2.0 * GE * F1 / &
                     (DE * RA**2) * OME(I) * XJA
            S22(3) = S22(3) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * XNU * ETA2 * RD1 * RD2 + &
                     AL * 2.0 * ETA1 * RD2**2 - FA * ETA1) * 2.0 * GE * F2 / (DE * RA**2) * OME(I) * XJA
            S22(4) = S22(4) + (2.0 * RDN * (AL * RD2 + 2.0 * XNU * RD2 - 4.0 * RD2**3) + 4.0 * XNU * ETA2 * &
                     RD2**2 + AL * (2.0 * ETA2 * RD2**2 + 2.0 * ETA2) - FA * ETA2) * 2.0 * GE * F2 / &
                     (DE * RA**2) * OME(I) * XJA
            S22(5) = S22(5) + (2.0 * RDN * (AL * RD1 - 4.0 * RD1 * RD2**2) + 4.0 * XNU * ETA2 * RD1 * RD2 + &
                     AL * 2.0 * ETA1 * RD2**2 - FA * ETA1) * 2.0 * GE * F3 / (DE * RA**2) * OME(I) * XJA
            S22(6) = S22(6) + (2.0 * RDN * (AL * RD2 + 2.0 * XNU * RD2 - 4.0 * RD2**3) + 4.0 * XNU * ETA2 * &
                     RD2**2 + AL * (2.0 * ETA2 * RD2**2 + 2.0 * ETA2) - FA * ETA2) * 2.0 * GE * F3 / &
                     (DE * RA**2) * OME(I) * XJA
        end do

    end subroutine SIGMAEQ

end module assembly_mod
