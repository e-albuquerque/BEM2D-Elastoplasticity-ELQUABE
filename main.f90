!-----------------------------------------------------------------------
! PROGRAM: ELQUABE (Modernized 2D Elastoplastic BEM Version)
!> @brief Solves two-dimensional elastic and elastoplastic problems using
!>        quadratic boundary elements and 2D domain cells (6-node triangles
!>        and 8-node quadrilaterals).
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
    INP = 15
    IPR = 16
    IP1 = 17
    IP2 = 18
    IP3 = 19
    IP4 = 20

    ! Open files for input and output
    FILEIN  = 'eq12cb.dat'
    FILEOUT = 'eq12cb.out'

    open(INP, file=trim(FILEIN),  status='OLD')
    open(IPR, file=trim(FILEOUT), status='UNKNOWN') ! General output report
    open(IP1, file='dpbou.dat',   status='UNKNOWN') ! Boundary displacements
    open(IP2, file='dpint.dat',   status='UNKNOWN') ! Internal displacements
    open(IP3, file='stbou.dat',   status='UNKNOWN') ! Boundary stresses & plastic strain
    open(IP4, file='stint.dat',   status='UNKNOWN') ! Internal stresses & plastic strain

    ! 1. Read data and dynamically allocate arrays
    call INPUTEQ()

    ! 2. Compute G and H matrices and form system (A * X = F)
    call GHMATEQ()

    ! 3. Solve elastic boundary system of equations
    call SLNPD(H, FI, D, NN, NX, IPR)

    ! 4. Compute elastic stress and displacement values at boundary & internal points
    call INTEREQ()

    ! 5. If internal cells exist (NCELL > 0), assemble domain cell matrices
    !    and run the 2D incremental-iterative elastoplastic solver
    if (NCELL > 0) then
        call NL_COEFS()
        call NL_SOLVE()
    end if

    ! 6. Print results at boundary nodes and internal points
    call OUTPTEQ()

    ! 7. Deallocate dynamic memory and close files
    call deallocate_arrays()

    close(INP)
    close(IPR)
    close(IP1)
    close(IP2)
    close(IP3)
    close(IP4)

end program ELQUABE
