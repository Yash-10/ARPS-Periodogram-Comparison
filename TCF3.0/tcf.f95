!!! Calculate TCF Filter 
!!! 
!!! uses ad-hoc mod for longer memory of spikes
!!!

module tcf_tools
  use m_median
  implicit none

contains 
  ! subroutine tcf(y,n,na,minp,maxp,splitper_temp,vsn,vdepth,vphase,vdur,vper)
!!!
  subroutine tcf(y,na,inper,outpow,outdepth,outphase,outdur,outmad)
!!!
    implicit none
    double precision, intent(in), dimension(:) :: y
    integer, intent(in), dimension(:) :: na
    double precision, intent(in), dimension(:) :: inper
    double precision, intent(out), dimension(:) :: outpow, outdepth,&
         & outphase, outdur,outmad

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

    durmax = 50
    durmem = 5
    qmax = 6                    !divide intper by qmax for dur
    nper = size(inper) 
    ny = size(y)
    Period: do iper = 1, nper  !loop over periods 
       per = inper(iper)        !current period
       intper =  ceiling(per)
       allocate( phdepth(intper) )
       allocate( phdepthcnt(intper) )
       allocate( phdur(intper) )
       allocate( phdepthma(intper) )
       Phase: do iph = 1, intper !nph   !loop over phases
          ph = dble(iph)        !current phase
          ispk = 0              !i-th comb spike
          nacnt = 0             !NA count
          rind = ph          
          ind = ceiling(rind)   !array index for i-th spike
          prod = 0.d0
          Comb: do while (ind < ny) !loop over comb spikes
             prod = prod + y(ind) !+ y(ind+1)+ y(ind+2)!+ y(ind+3)
             nacnt = nacnt + na(ind) !+ na(ind+1)+ na(ind+2)!+ na(ind+3) !count NAs
             ispk = ispk + 1
             rind = per*dble(ispk) + ph
             ind = ceiling(rind)
          end do Comb
          phdepthcnt(iph) = max((ispk - nacnt),1)
          phdepth(iph) = prod/phdepthcnt(iph)
       end do Phase
       !!
       depth = 0.d0
       depth_best = 0.d0
       ph_best = 0
       dur_best = 0
       ing_best = 0
       phdepthma = phdepth
       Duration: do idur = 1, max(min(intper/qmax,durmax),1) !loop over transit duration
          !! idur = 1
          if ( (idur>1) .and. (idur<=durmem) ) then
!!$             phdur(1:(intper-idur)) =  phdepth((1+idur):intper)
!!$             phdur((intper-idur+1):intper) = phdepth(1:idur)
!!$             phdepthma = phdepthma + (phdur-phdepthma)/dble(idur)
             phdepthma = phdepthma + (cshift(phdepth,idur)-phdepthma)/dble(idur)
         end if
          !Below: a tad faster than the cshift implementation
          phdur(1:(intper-idur)) =  -phdepthma(1:(intper-idur)) + phdepthma((1+idur):intper)
          phdur((intper-idur+1):intper) = -phdepthma((intper-idur+1):intper) + phdepthma(1:idur)
!!$          phdur=  -phdepthma + cshift(phdepthma,idur)
          !!
          ing = maxloc(phdur, dim=1) 
          depth = phdur(ing)
          if(depth > depth_best) then
             depth_best = depth !estimated transit depth
             ph_best = ing      !estimated transit phase
             dur_best = idur    !estimated transit duration
             ing_best = ing     !estimated transit ingress
          end if
       end do Duration
       egr_best = mod(ing_best+dur_best,intper)
       if(egr_best==0) egr_best = intper
       !!
       depth = 0.d0
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
       pow_best = depth/rcnt*depth !least squares power
       depth_best = depth/rcnt*2.d0
       !!
       outpow(iper) = pow_best
       outdepth(iper) = depth_best
       outphase(iper) = dble(ph_best)
       outdur(iper) = dble(dur_best)
       outmad(iper) = median(abs(phdepth))
       !!
       deallocate(phdepth)
       deallocate(phdepthcnt)
       deallocate(phdur)
       deallocate(phdepthma)
    end do Period

  end subroutine tcf

end module tcf_tools


