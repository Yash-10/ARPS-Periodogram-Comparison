!!! 
!!! Calculate TCF Filter 
!!! 
!!!

module tcf_tools
  use m_median
  implicit none

contains 
  subroutine tcf(y,na,inper,outpow,outdepth,outphase,outdur,outmad)

    implicit none

    ! Arguments:
    ! Inputs vectors
    double precision, intent(in), dimension(:) :: y         ! time series; NAs imputed with zeros
    integer, intent(in), dimension(:) :: na                 ! 0-1 label of NA values (NA = 1); same size as `y`
    double precision, intent(in), dimension(:) :: inper     ! periods to test
    ! Output vectors
    ! each same size as `inper`, provide estimated values for each period tested 
    double precision, intent(out), dimension(:) :: outpow   ! periodogram power
    double precision, intent(out), dimension(:) :: outdepth ! transit depth
    double precision, intent(out), dimension(:) :: outphase ! transit phase
    double precision, intent(out), dimension(:) :: outdur   ! transit duration
    double precision, intent(out), dimension(:) :: outmad   ! 

    ! Internal variables
    double precision, dimension(:), allocatable :: phdepth,phdepthma
    double precision, dimension(:), allocatable :: phdur
    integer, dimension(:), allocatable :: phdepthcnt

    double precision prod, depth, pow_best, depth_best
    double precision rind, rcnt
    integer ind
    integer nper, ny
    double precision per, ph
    integer i
    integer iper, iph, idur, ispk, nacnt, intper
    integer dur_best, ph_best
    integer ing, egr, ing_best, egr_best
    integer durmem, durmax, qmax

    ! *********
    ! Set parameters:
    durmax = 50   !maximum transit duration in units of cadence
    qmax = 6      !maximum period/duration ratio (to set max possible duration)
    durmem = 5    !spike width; ad-hoc value for ingress/egress "memory"
    ! *********

    nper = size(inper)   !number of input periods
    ny = size(y)         !time series length (# obs)

    ! Loop over periods
    Period: do iper = 1, nper
       per = inper(iper)        !current period
       intper =  ceiling(per)   !round up to integer cadence
       allocate( phdepth(intper) )
       allocate( phdepthcnt(intper) )
       allocate( phdur(intper) )
       allocate( phdepthma(intper) )

       ! Loop over phases; integer-valued in cadence units, up to period
       Phase: do iph = 1, intper
          ph = dble(iph)         !change current phase type for calculations

          ! Initialize some counters
          ispk = 0              !comb spike count
          nacnt = 0             !NA count
          ind = iph             !spike index in array; starts at current phase          
          prod = 0.d0           !dot product of filter and time series

          ! Loop over comb spikes to calculate dot product
          ! since spike height is 1, this means just summing over time series values at each spike
          Comb: do while (ind < ny) !stop if next spike would be beyond length of time series
             prod = prod + y(ind)
             nacnt = nacnt + na(ind)
             ispk = ispk + 1
             rind = per*dble(ispk) + ph !time of next spike
             ind = ceiling(rind)        !round off to the next cadence
          end do Comb
          phdepthcnt(iph) = max((ispk - nacnt),1) !count of non-NA spikes
          phdepth(iph) = prod/phdepthcnt(iph)     !average depth of non-NA spikes
       end do Phase

       ! Now we have values for all (integer) phases of the current period
       ! next we need to test transit durations and combine down-up spikes

       ! Some temp variables
       depth = 0.d0
       phdepthma = phdepth      !used for moving average of depth over spike width ("memory")

       ! Initialize estimates
       depth_best = 0.d0        !optimal depth*** sum (not average, so technically not real "depth")
       ph_best = 0              !optimal phase (in cadences)
       dur_best = 0             !optimal duration (in cadences)
       ing_best = 0             !optimal ingress location (equal to `ph_best`, separate for convenience)

       ! Loop over (integer) transit durations, up to maximum period/duration ratio
       Duration: do idur = 1, max(min(intper/qmax,durmax),1)

          ! Average values of spike width, up to `durmem`, as long as not longer than current duration tested 
          if ( (idur>1) .and. (idur<=durmem) ) then
             phdepthma = phdepthma + (cshift(phdepth,idur)-phdepthma)/dble(idur)
          end if

          ! Combine down-up spikes using `idur` duration for all phases
          phdur(1:(intper-idur)) =  -phdepthma(1:(intper-idur)) + phdepthma((1+idur):intper)
          phdur((intper-idur+1):intper) = -phdepthma((intper-idur+1):intper) + phdepthma(1:idur)

          ! Get optimal values
          ing = maxloc(phdur, dim=1) !phase index where filter is maximized for test duration
          depth = phdur(ing)         !depth of combined down-up spikes at best phase (for current duration)

          ! Save value if better than previous ones
          if(depth > depth_best) then
             depth_best = depth
             ph_best = ing
             dur_best = idur
             ing_best = ing
          end if
       end do Duration

       ! Get egress index from ingress and duration
       egr_best = mod(ing_best+dur_best,intper)
       if(egr_best==0) egr_best = intper

       ! Loop to calculate periodogram "power"
       ! needs to be run separately because we need to use the number of observations per spike
       ! in principle can be incorporated into duration loop, kept here for simplicity
       depth = 0.d0  !calling this "depth" here technically is a misnomer
       rcnt = 0.d0
       do i = 0, min(dur_best,durmem)-1
          ing = mod(ing_best+i,intper)
          egr = mod(egr_best+i,intper)
          if(ing==0) ing = intper
          if(egr==0) egr = intper
          depth = depth - phdepth(ing)*dble(phdepthcnt(ing))
          depth = depth + phdepth(egr)*dble(phdepthcnt(egr))
          rcnt = rcnt + phdepthcnt(ing)+phdepthcnt(egr)
       end do
       pow_best = depth/rcnt*depth  !least squares power
       depth_best = depth/rcnt*2.d0 !optimal depth corresponding to max power

       ! Save best estimates at current period
       outpow(iper) = pow_best
       outdepth(iper) = depth_best
       outphase(iper) = dble(ph_best)
       outdur(iper) = dble(dur_best)
       outmad(iper) = median(abs(phdepth))

       ! Clean up temporary arrays
       deallocate(phdepth)
       deallocate(phdepthcnt)
       deallocate(phdur)
       deallocate(phdepthma)

    end do Period

  end subroutine tcf

end module tcf_tools


