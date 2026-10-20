!-----------------------------------------------------------------------
! MODULE: io_mod
!> @brief Input and Output routines for ELQUABE (INPUTEQ and OUTPTEQ).
!-----------------------------------------------------------------------
module io_mod
    use data_mod
    implicit none

contains

    !-------------------------------------------------------------------
    !> @brief Reads mesh, material properties, boundary conditions,
    !>        and internal points from input file, allocating arrays.
    !-------------------------------------------------------------------
    subroutine INPUTEQ()
        implicit none
        character(len=80) :: TITLE
        integer :: I, J

        write(IPR, 100)
100     format(' ', 79('*'))

        ! Read job title
        read(INP, '(A)') TITLE
        write(IPR, '(A)') TITLE

        ! Read number of boundary elements, internal points, and material properties
        read(INP, *) NE, L, GE, XNU
        write(IPR, 210) NE, L, GE, XNU
210     format(//2X, 'DATA'/2X, 'NUMBER OF BOUNDARY ELEMENTS=', &
               I3/2X, 'NUMBER OF INTERNAL POINTS=', I3/2X, 'SHEAR MODULUS=', &
               E14.7/2X, 'POISSON RATIO=', E14.7)

        ! Dynamically allocate arrays based on NE and L
        call allocate_arrays()

        ! Read boundary nodes coordinates
        write(IPR, 500)
        read(INP, *) (X(I), Y(I), I = 1, N)
        do I = 1, N
            write(IPR, 240) I, X(I), Y(I)
        end do
500     format(//2X, 'BOUNDARY NODES COORDINATES'///4X, &
               'NODE', 10X, 'X', 18X, 'Y'/)
240     format(5X, I3, 2(5X, E14.7))

        ! Read boundary conditions in DFI(I) vector:
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

        ! Read internal points coordinates
        if (L > 0) then
            read(INP, *) (CX(I), CY(I), I = 1, L)
        end if

    end subroutine INPUTEQ

    !-------------------------------------------------------------------
    !> @brief Prints displacements, tractions, and stresses at boundary
    !>        nodes, and displacements and stresses at internal points.
    !-------------------------------------------------------------------
    subroutine OUTPTEQ()
        implicit none
        integer :: I, K

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
        write(IPR, 250) X(1), Y(1), DFI(6*NE-1), DFI(6*NE), DFI(1), DFI(2)
        write(IPR, 250) X(2), Y(2), DFI(3), DFI(4), DFI(3), DFI(4)

        do I = 2, NE
            write(IPR, 250) X(2*I-1), Y(2*I-1), DFI(6*I-7), DFI(6*I-6), &
                            DFI(6*I-5), DFI(6*I-4)
            write(IPR, 250) X(2*I), Y(2*I), DFI(6*I-3), DFI(6*I-2), &
                            DFI(6*I-3), DFI(6*I-2)
        end do
250     format(6(1X, E12.5))

        ! Write boundary stresses [sigma_x, tau_xy, sigma_y, sigma_vm] to stbou.dat and eq12cb.out
        write(IPR, 260)
260     format(///2X, 'BOUNDARY NODES STRESSES'//8X, 'X', 15X, 'Y', 12X, &
               'SIGMA X', 10X, 'TAU XY', 9X, 'SIGMA Y', 8X, 'VON MISES')
        do I = 1, N
            write(IP3, *) SBOU(I, 1), SBOU(I, 2), SBOU(I, 3), SBOU(I, 4)
            write(IPR, 270) X(I), Y(I), SBOU(I, 1), SBOU(I, 2), SBOU(I, 3), SBOU(I, 4)
        end do
270     format(6(2X, E14.7))

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
350         format(//2X, 'INTERNAL POINTS STRESSES'//8X, 'X', 15X, 'Y', 12X, &
                   'SIGMA X', 10X, 'TAU XY', 9X, 'SIGMA Y', 8X, 'VON MISES')
            do K = 1, L
                write(IP4, *) SSOL(3*K-2), SSOL(3*K-1), SSOL(3*K), SVM_INT(K)
                write(IPR, 450) CX(K), CY(K), SSOL(3*K-2), SSOL(3*K-1), SSOL(3*K), SVM_INT(K)
            end do
        end if
400     format(4(2X, E14.7))
450     format(6(2X, E14.7))

        write(IPR, 500)
500     format(' ', 79('*'))

    end subroutine OUTPTEQ

end module io_mod
