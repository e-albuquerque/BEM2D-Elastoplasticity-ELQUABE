!-----------------------------------------------------------------------
! PROGRAM: ELQUABE (Modernized Fortran 90 Version)
!> @brief Solves two-dimensional elastic problems using quadratic
!>        boundary elements (Brebbia & Dominguez).
!-----------------------------------------------------------------------
program ELQUABE
    use data_mod
    use io_mod
    use assembly_mod
    use solver_mod
    implicit none

    character(len=64) :: FILEIN, FILEOUT
    real :: D

    ! Assign unit numbers for input and output files
    INP = 5
    IPR = 6
    IP1 = 7
    IP2 = 8
    IP3 = 9
    IP4 = 10

    ! Open files for input and output
    FILEIN  = 'eq12cb.dat'
    FILEOUT = 'eq12cb.out'

    open(INP, file=trim(FILEIN),  status='OLD')
    open(IPR, file=trim(FILEOUT), status='UNKNOWN') ! General output report
    open(IP1, file='dpbou.dat',   status='UNKNOWN') ! Boundary displacements
    open(IP2, file='dpint.dat',   status='UNKNOWN') ! Internal displacements
    open(IP3, file='stbou.dat',   status='UNKNOWN') ! Boundary stresses/tractions
    open(IP4, file='stint.dat',   status='UNKNOWN') ! Internal stresses

    ! 1. Read data and dynamically allocate arrays
    call INPUTEQ()

    ! 2. Compute G and H matrices and form system (A * X = F)
    call GHMATEQ()

    ! 3. Solve system of equations
    call SLNPD(H, FI, D, NN, NX, IPR)

    ! 4. Compute stress and displacement values at internal points
    call INTEREQ()

    ! 5. Print results at boundary nodes and internal points
    call OUTPTEQ()

    ! 6. Deallocate dynamic memory and close files
    call deallocate_arrays()

    close(INP)
    close(IPR)
    close(IP1)
    close(IP2)
    close(IP3)
    close(IP4)

end program ELQUABE
