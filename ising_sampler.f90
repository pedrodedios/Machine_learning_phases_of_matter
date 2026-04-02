!====================================================================================
! ising_sampler.f90
!
! Stripped-down Metropolis sampler for the 2D Ising model.
! Purpose : generate decorrelated spin configurations for ML training.
!
! Output  : one binary file per temperature  configs_L{L}_T{T:.4f}.bin
!           each file holds  N_CONFIGS × (XS*YS)  signed bytes  (±1 as int8)
!           plus a human-readable index file  configs_index.dat
!
! Python loader (one line):
!   cfgs = np.fromfile('configs_L20_T2.0000.bin', dtype=np.int8).reshape(N, L*L)
!====================================================================================






module ising_utils
    implicit none
    contains

    pure integer function shiftup(num, NS)
        integer, intent(in) :: num, NS
        shiftup = mod(num, NS) + 1          ! cleaner than if-branch
    end function shiftup

    pure integer function shiftdown(num, NS)
        integer, intent(in) :: num, NS
        shiftdown = num - 1
        if (shiftdown == 0) shiftdown = NS
    end function shiftdown

end module ising_utils









!------------------------------------------------------------------------------------
! Core subroutine: thermalize once, then collect n_configs decorrelated snapshots.
! Each snapshot is written as a flat row of XS*YS int8 values (±1) via stream I/O,
! so Python can read it with np.fromfile — no Fortran record markers in the file.
!------------------------------------------------------------------------------------

subroutine sample_configs(temperature, XS, YS, n_therm, n_between, n_configs, out_unit)

    use ising_utils
    implicit none

    real(8),    intent(in) :: temperature
    integer,    intent(in) :: XS, YS, n_therm, n_between, n_configs, out_unit

    integer(1), allocatable :: lattice(:,:)   ! int8 directly — no conversion on write
    integer :: i, j, s, c
    integer :: up_x, down_x, up_y, down_y
    real(8) :: beta, delta, nb_sum, p, rnd
    integer(1) :: flipped

    allocate(lattice(XS, YS))
    beta = 1.0d0 / temperature


    !── Hot start ─────────────────────────────────────────────────────────
    do j = 1, YS
        do i = 1, XS
            call random_number(rnd)
            lattice(i,j) = merge(1_1, -1_1, rnd >= 0.5d0)
        end do
    end do


    !── Thermalization (discard these sweeps entirely) ────────────────────
    do s = 1, n_therm
        call sweep(lattice, XS, YS, beta)
    end do

    !── Collect n_configs decorrelated configurations ─────────────────────


    do c = 1, n_configs

        ! Advance chain n_between sweeps between each saved sample
        do s = 1, n_between
            call sweep(lattice, XS, YS, beta)
        end do

        ! Write flat row of XS*YS int8 values 
        write(out_unit) lattice   ! Fortran writes column-major; Python reshape handles this

    end do


    deallocate(lattice)




contains

    !── Single Metropolis sweep over all sites 

    subroutine sweep(lat, NX, NY, bt)

        integer(1), intent(inout) :: lat(NX, NY)
        integer,    intent(in)    :: NX, NY
        real(8),    intent(in)    :: bt
        integer :: ii, jj, ux, dx, uy, dy
        real(8) :: dE, nb, pp, rr

        do jj = 1, NY
            do ii = 1, NX
                ux = shiftup(ii, NX);   dx = shiftdown(ii, NX)
                uy = shiftup(jj, NY);   dy = shiftdown(jj, NY)

                nb  = real(lat(ii,uy) + lat(ii,dy) + lat(ux,jj) + lat(dx,jj), 8)
                dE  = 2.0d0 * real(lat(ii,jj), 8) * nb

                if (dE <= 0.0d0) then
                    lat(ii,jj) = -lat(ii,jj)
                else
                    pp = exp(-bt * dE)
                    call random_number(rr)
                    if (rr < pp) lat(ii,jj) = -lat(ii,jj)
                end if
            end do
        end do
    end subroutine sweep


end subroutine sample_configs









!====================================================================================
program run_sampler

    use ising_utils
    implicit none

    !── Lattice ─────────────────────────────────────────────────────────────────────
    integer, parameter :: XS = 30, YS = 30   ! system size (must be square for the NN)


    !── Sampling parameters ─────────────────────────────────────────────────────────
    integer, parameter :: N_THERM   = 40000   ! thermalization sweeps (increase near Tc)
    integer, parameter :: N_BETWEEN = 160     ! decorrelation sweeps between saved configs
    integer, parameter :: N_CONFIGS = 1000    ! configs to save per temperature



    real(8), parameter :: T_START = 1.0d0
    real(8), parameter :: T_END   = 3.6d0
    real(8), parameter :: DELTA_T = 0.08d0




    !── Local vars ──────────────────────────────────────────────────────────────────
    real(8)       :: temperature
    integer       :: n, T_steps, file_unit, idx_unit
    character(64) :: filename

   
     call random_seed()   ! uncomment to use system time seed

    T_steps  = nint((T_END - T_START) / DELTA_T)
    idx_unit = 99
    open(unit=idx_unit, file='configs_index.dat', status='replace')
    write(idx_unit, '(A)') '# filename   T   XS   YS   n_configs   n_therm   n_between'

    do n = 0, T_steps

        temperature = T_START + n * DELTA_T

        !── Build filename and open stream binary file ──────────────────────────────
        write(filename, '("configs_L",I0,"_T",F6.4,".bin")') XS, temperature
        file_unit = 10 + n

        open(unit   = file_unit,      &
             file   = trim(filename), &
             status = 'replace',      &
             access = 'stream',       &   ! no Fortran record markers
             form   = 'unformatted')

        print '(A,F5.3,A)', '  Sampling T = ', temperature, ' ...'

        call sample_configs(temperature, XS, YS, &
                            N_THERM, N_BETWEEN, N_CONFIGS, file_unit)

        close(file_unit)

        !── Write index entry ───────────────────────────────────────────────────────
        write(idx_unit, '(A,1X,F6.4,5(1X,I6))') &
            trim(filename), temperature, XS, YS, N_CONFIGS, N_THERM, N_BETWEEN

        print '(A,A)', '  Saved: ', trim(filename)

    end do

    close(idx_unit)
    print *, ''
    print *, 'Done. See configs_index.dat for file listing.'



end program run_sampler
