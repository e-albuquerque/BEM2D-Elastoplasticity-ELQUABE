!-----------------------------------------------------------------------
! MODULE: io_mod
!> @brief Input and Output routines for the 2D Elastoplastic BEM solver
!>        (INPUTEQ and OUTPTEQ).
!-----------------------------------------------------------------------
module io_mod
    use data_mod
    implicit none

    real, parameter :: PI_IO = 3.141592653589793

contains

    !-------------------------------------------------------------------
    !> @brief Reads mesh, boundary connectivity, material properties,
    !>        boundary conditions, internal points, internal cells (6-node
    !>        triangles or 8-node quadrilaterals), and plasticity parameters.
    !-------------------------------------------------------------------
    subroutine INPUTEQ()
        implicit none
        character(len=80) :: TITLE
        integer :: I, J, ICELL, IOS, TMP_ND
        real    :: SINF, X1, Y1, X2, Y2, X3, Y3, CROSS

        write(IPR, 100)
100     format(' ', 79('*'))

        ! 1. Read job title
        read(INP, '(A)') TITLE
        write(IPR, '(A)') TITLE

        ! 2. Read number of boundary elements, internal points, and elastic properties
        read(INP, *) NE, L, GE, XNU
        write(IPR, 210) NE, L, GE, XNU
210     format(//2X, 'DATA'/2X, 'NUMBER OF BOUNDARY ELEMENTS=', &
               I4/2X, 'NUMBER OF INTERNAL POINTS=', I4/2X, 'SHEAR MODULUS=', &
               E14.7/2X, 'POISSON RATIO=', E14.7)

        ! Dynamically allocate elastic arrays based on NE and L
        call allocate_arrays()

        ! 3. Read boundary nodes coordinates
        write(IPR, 500)
        read(INP, *) (X(I), Y(I), I = 1, N)
        do I = 1, N
            CD(1, I) = X(I)
            CD(2, I) = Y(I)
            write(IPR, 240) I, X(I), Y(I)
        end do
500     format(//2X, 'BOUNDARY NODES COORDINATES'///4X, &
               'NODE', 10X, 'X', 18X, 'Y'/)
240     format(4X, I4, 2(5X, E14.7))

        ! 4. Read boundary element connectivity matrix CONEC(NE, 3)
        write(IPR, 600)
600     format(//2X, 'ELEMENT CONNECTIVITY'//4X, &
               'ELE', 6X, 'NODE 1', 6X, 'NODE 2', 6X, 'NODE 3'/)
        do I = 1, NE
            read(INP, *) (CONEC(I, J), J = 1, 3)
            write(IPR, 650) I, (CONEC(I, J), J = 1, 3)
        end do
650     format(3X, I4, 3(6X, I6))

        ! 5. Read boundary conditions in DFI(I) vector:
        ! If KODE(I) == 0 -> DFI(I) is a known displacement
        ! If KODE(I) == 1 -> DFI(I) is a known traction
        write(IPR, 800)
800     format(//2X, 'BOUNDARY CONDITIONS'//35X, 'PRESCRIBED VALUES' / &
               7X, '-----FIRST NODE------', 3X, '-----SECOND NODE-----', 3X, &
               '-----THIRD NODE------'/ &
               1X, 'ELE', 3X, 'X DIR.', 2X, 'C', 3X, 'Y DIR.', 2X, 'C', 3X, &
               'X DIR.', 2X, 'C', 3X, 'Y DIR.', 2X, 'C', 3X, 'X DIR.', &
               2X, 'C', 3X, 'Y DIR.', 2X, 'C'/)
        do I = 1, NE
            read(INP, *) (KODE(6*I-6+J), DFI(6*I-6+J), J = 1, 6)
            write(IPR, 950) I, (DFI(6*I-6+J), KODE(6*I-6+J), J = 1, 6)
        end do
950     format(1X, I3, 6(F9.3, 2X, I1))

        ! 6. Read internal points coordinates
        if (L > 0) then
            read(INP, *) (CX(I), CY(I), I = 1, L)
            do I = 1, L
                CD(1, N + I) = CX(I)
                CD(2, N + I) = CY(I)
            end do
        end if

        ! 7. Read internal cell control parameters (if present)
        NCELL = 0
        NODC  = 6
        NSIG  = 4
        TOLGP = 1.0e-5
        NAUTO = -6
        read(INP, *, iostat=IOS) NCELL, NODC, NSIG, TOLGP, NAUTO
        if (IOS /= 0) then
            NCELL = 0
        end if

        ! Set up 2D elastic constants (Plane Strain NSIG=4 or Plane Stress NSIG=3)
        PR = XNU
        if (NSIG == 3) PR = XNU / (1.0 + XNU)
        PR1   = 1.0 - PR
        PR2   = 1.0 - 2.0 * PR
        PR3   = 3.0 - 3.0 * PR
        PR4   = 1.0 - 4.0 * PR
        PR5   = 3.0 - 4.0 * PR
        GAM   = 4.0
        CON   = -1.0 / (4.0 * PI_IO * PR1)
        ROOT3 = sqrt(3.0)
        CLAME = PR / PR2
        CJUP  = 2.0

        ! Allocate plasticity arrays
        call allocate_plastic_arrays()

        if (NCELL > 0) then
            write(IPR, 1000) NCELL, NODC, NSIG, TOLGP, NAUTO
1000        format(///2X, '2D ELASTOPLASTIC PART (INTERNAL CELLS)'// &
                   2X, 'NUMBER OF INTERNAL CELLS (NCELL) =', I5/ &
                   2X, 'NODES PER INTERNAL CELL  (NODC)  =', I5/ &
                   2X, 'STRESS COMPONENTS        (NSIG)  =', I5/ &
                   2X, 'GAUSS TOLERANCE          (TOLGP) =', E13.6/ &
                   2X, 'GAUSS CONTROL            (NAUTO) =', I5// &
                   4X, 'CELL', 15X, 'CONNECTIVITY NODES')

            ! 8. Read internal cells connectivity LNDC(NODC, NCELL)
            do ICELL = 1, NCELL
                read(INP, *) (LNDC(J, ICELL), J = 1, NODC)

                ! Ensure counter-clockwise orientation of each 2D cell
                X1 = CD(1, LNDC(1, ICELL)); Y1 = CD(2, LNDC(1, ICELL))
                X2 = CD(1, LNDC(2, ICELL)); Y2 = CD(2, LNDC(2, ICELL))
                X3 = CD(1, LNDC(3, ICELL)); Y3 = CD(2, LNDC(3, ICELL))
                CROSS = (X2 - X1) * (Y3 - Y1) - (Y2 - Y1) * (X3 - X1)
                if (CROSS < 0.0) then
                    if (NODC == 6) then
                        TMP_ND = LNDC(2, ICELL); LNDC(2, ICELL) = LNDC(3, ICELL); LNDC(3, ICELL) = TMP_ND
                        TMP_ND = LNDC(4, ICELL); LNDC(4, ICELL) = LNDC(6, ICELL); LNDC(6, ICELL) = TMP_ND
                    else if (NODC == 8) then
                        TMP_ND = LNDC(2, ICELL); LNDC(2, ICELL) = LNDC(4, ICELL); LNDC(4, ICELL) = TMP_ND
                        TMP_ND = LNDC(5, ICELL); LNDC(5, ICELL) = LNDC(8, ICELL); LNDC(8, ICELL) = TMP_ND
                        TMP_ND = LNDC(6, ICELL); LNDC(6, ICELL) = LNDC(7, ICELL); LNDC(7, ICELL) = TMP_ND
                    end if
                end if

                write(IPR, '(2X, I5, 2X, 8I6)') ICELL, (LNDC(J, ICELL), J = 1, NODC)
            end do

            ! 9. Read material plasticity properties
            read(INP, *) NCRIT, UNIAX, FRICT, HARDI, HARDK, EM
            write(IPR, 1020) NCRIT, UNIAX, FRICT, HARDI, HARDK, EM
1020        format(/2X, 'NCRIT', 6X, 'UNIAX', 8X, 'FRICT', 8X, 'HARDI', 8X, 'HARDK', &
                   9X, 'EM'/2X, I3, 3X, 1P5E13.5)

            ! 10. Read incremental-iterative parameters
            read(INP, *) NINCS, FINCR, MITER, TOLER
            write(IPR, 1030) NINCS, FINCR, MITER, TOLER
1030        format(/2X, 'NINCS     FINCR    MITER     TOLER'/ &
                   2X, I4, 4X, F8.6, I7, 1PE14.4)

            ! Normalize plasticity parameters by 2*GE (as in BEMECH)
            FRICT = FRICT * (PI_IO / 180.0)
            SINF  = sin(FRICT)
            ALFA  = 0.0
            CEQS  = 1.0
            if (NCRIT == 4) ALFA = 2.0 * SINF / (ROOT3 * (3.0 - SINF))
            if (NCRIT == 3) CEQS = 0.5 * (1.0 - SINF)
            if (NCRIT == 4) CEQS = 1.0 / ROOT3 - ALFA
            if (NCRIT == 3) UNIAX = UNIAX * cos(FRICT) / CEQS
            if (NCRIT == 4) UNIAX = (6.0 * UNIAX * cos(FRICT) / (ROOT3 * (3.0 - SINF))) / CEQS
            CEQP = sqrt(2.0 / 3.0)
            if (NCRIT == 3) CEQP = sqrt(2.0) * (1.0 - SINF) / sqrt(3.0 - 2.0 * SINF + 3.0 * SINF**2)
            if (NCRIT == 4) CEQP = (1.0 / ROOT3 - ALFA) / sqrt(3.0 * ALFA * ALFA + 0.5)

            UNIAX = UNIAX / (2.0 * GE)
            HARDI = EM * HARDI / (2.0 * GE)
            HARDK = (1.0 - EM) * HARDK / (2.0 * GE)
            HARDS = HARDI + HARDK
        end if

    end subroutine INPUTEQ

    !-------------------------------------------------------------------
    !> @brief Prints displacements, tractions, stresses, and equivalent
    !>        plastic strains at boundary nodes and internal points.
    !-------------------------------------------------------------------
    subroutine OUTPTEQ()
        implicit none
        integer :: I, K, EL
        real :: TX_BEF, TY_BEF, TX_AFT, TY_AFT

        write(IPR, 100)
100     format(' ', 79('*')//1X, 'RESULTS'//2X, 'BOUNDARY NODES'///10X, &
               'X', 17X, 'Y', 11X, 'DISPLACEMENT X', 4X, 'DISPLACEMENT Y'/)
        do I = 1, N
            write(IP1, *) FI(2*I-1), FI(2*I), sqrt(FI(2*I-1)**2 + FI(2*I)**2)
            write(IPR, 200) X(I), Y(I), FI(2*I-1), FI(2*I)
        end do
200     format(4(4X, E14.7))

        write(IPR, 150)
150     format(///28X, 'TRACTION X', 3X, 'TRACTION Y', 3X, 'TRACTION X', &
               3X, 'TRACTION Y'/7X, 'X', 12X, 'Y', 7X, 'BEFORE NODE', 2X, 'BEFORE NODE', &
               2X, 'AFTER NODE', 3X, 'AFTER NODE'/)

        do I = 1, N
            TX_BEF = 0.0; TY_BEF = 0.0
            TX_AFT = 0.0; TY_AFT = 0.0
            do EL = 1, NE
                if (CONEC(EL, 2) == I) then
                    TX_BEF = DFI(6*EL - 3)
                    TY_BEF = DFI(6*EL - 2)
                    TX_AFT = DFI(6*EL - 3)
                    TY_AFT = DFI(6*EL - 2)
                else
                    if (CONEC(EL, 3) == I) then
                        TX_BEF = DFI(6*EL - 1)
                        TY_BEF = DFI(6*EL)
                    end if
                    if (CONEC(EL, 1) == I) then
                        TX_AFT = DFI(6*EL - 5)
                        TY_AFT = DFI(6*EL - 4)
                    end if
                end if
            end do
            write(IPR, 250) X(I), Y(I), TX_BEF, TY_BEF, TX_AFT, TY_AFT
        end do
250     format(6(1X, E12.5))

        ! Write boundary stresses [sigma_x, tau_xy, sigma_y, sigma_vm, epstn] to stbou.dat and eq12cb.out
        write(IPR, 260)
260     format(///2X, 'BOUNDARY NODES STRESSES AND PLASTIC STRAINS'//8X, 'X', 15X, 'Y', 12X, &
               'SIGMA X', 10X, 'TAU XY', 9X, 'SIGMA Y', 8X, 'VON MISES', 6X, 'EQ-P-STRAIN')
        do I = 1, N
            write(IP3, *) SBOU(I, 1), SBOU(I, 2), SBOU(I, 3), SBOU(I, 4), EPSTN(I)
            write(IPR, 270) X(I), Y(I), SBOU(I, 1), SBOU(I, 2), SBOU(I, 3), SBOU(I, 4), EPSTN(I)
        end do
270     format(7(2X, E14.7))

        if (L > 0) then
            write(IPR, 300)
300         format(///2X, 'INTERNAL POINTS DISPLACEMENTS'//8X, 'X', 15X, 'Y', 9X, &
                   'DISPLACEMENT X', 2X, 'DISPLACEMENT Y')
            do K = 1, L
                write(IPR, 400) CX(K), CY(K), DSOL(2*K-1), DSOL(2*K)
                write(IP2, *) DSOL(2*K-1), DSOL(2*K), &
                              sqrt(DSOL(2*K-1)**2 + DSOL(2*K)**2)
            end do

            write(IPR, 350)
350         format(//2X, 'INTERNAL POINTS STRESSES AND PLASTIC STRAINS'//8X, 'X', 15X, 'Y', 12X, &
                   'SIGMA X', 10X, 'TAU XY', 9X, 'SIGMA Y', 8X, 'VON MISES', 6X, 'EQ-P-STRAIN')
            do K = 1, L
                write(IP4, *) SSOL(3*K-2), SSOL(3*K-1), SSOL(3*K), SVM_INT(K), EPSTN(N + K)
                write(IPR, 450) CX(K), CY(K), SSOL(3*K-2), SSOL(3*K-1), SSOL(3*K), SVM_INT(K), EPSTN(N + K)
            end do
        end if
400     format(4(2X, E14.7))
450     format(7(2X, E14.7))

        write(IPR, 500)
500     format(' ', 79('*'))

    end subroutine OUTPTEQ

end module io_mod
