MODULE sbccpl_cesm
   !!======================================================================
   !!                       ***  MODULE  sbccpl_cesm  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in
   !!                               NCAR CESM infrastructure
   !!======================================================================
   !! History :  3.3  ! 2012 (P.G. Fogli, CMCC) Original code
   !!----------------------------------------------------------------------
#if defined CCSMCOUPLED
   !!----------------------------------------------------------------------
   !!   'CCSMCOUPLED'    CESM coupled formulation
   !!----------------------------------------------------------------------
   !!   sbc_cpl_cesm_init  : allocate space for the received fields
   !!   sbc_cpl_cesm_rcv   : receive fields from the atmosphere/sea ice over the ocean
   !!   sbc_cpl_cesm_finalize: deallocate space for the received fields
   !!
   USE par_kind
   USE par_oce
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! SBC variables
   USE sbcrnf
   USE phycst          ! physical constants
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE geo2ocean       ! 
   USE in_out_manager
   USE lib_fortran
   USE sbcapr, ONLY : apr
   !!
   IMPLICIT NONE
   PRIVATE

   LOGICAL, PARAMETER, PUBLIC :: lk_cesm = .true.

   ! all fields received from the atmosphere / sea ice
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: taux_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: tauy_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: evap_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: rain_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: snow_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: roff_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: ioff_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: meltw_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: salt_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: swnet_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: sen_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: lat_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: lwup_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: lwdn_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: melth_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: ifrac_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: pslv_x2o
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: duu10n_x2o
#if defined key_cpl_carbon_cycle
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: co2_x2o
#endif
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) :: sicew_x2o

   LOGICAL, SAVE, PUBLIC :: lrecv

   PUBLIC   sbc_cpl_cesm_init      ! routine called by ocn_init_mct (ocn_comp_mct.F90)
   PUBLIC   sbc_cpl_cesm_rcv       ! routine called by sbc (sbcmod.F90)
   PUBLIC   sbc_cpl_cesm_finalize  ! routine called by ocn_final_mct (ocn_comp_mct.F90)
   
   REAL(wp), SAVE :: garea

   !! * Substitutions
   ! for DO macro
#  include "do_loop_substitute.h90"

CONTAINS

   SUBROUTINE sbc_cpl_cesm_init

      INTEGER :: istat, ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( taux_x2o(jpi,jpj),  tauy_x2o(jpi,jpj),                       &
                evap_x2o(jpi,jpj),  rain_x2o(jpi,jpj),  snow_x2o(jpi,jpj),   &
                roff_x2o(jpi,jpj),  ioff_x2o(jpi,jpj),  meltw_x2o(jpi,jpj),  &
                salt_x2o(jpi,jpj),  swnet_x2o(jpi,jpj), sen_x2o(jpi,jpj),    &
                lat_x2o(jpi,jpj),   lwup_x2o(jpi,jpj),  lwdn_x2o(jpi,jpj),   &
                melth_x2o(jpi,jpj), ifrac_x2o(jpi,jpj), pslv_x2o(jpi,jpj),   &
                duu10n_x2o(jpi,jpj),  apr(jpi,jpj),     sicew_x2o(jpi,jpj),  &
                STAT=ierr(1) )
      !

#if defined key_cpl_carbon_cycle
      ALLOCATE( co2_x2o(jpi,jpj), STAT=ierr(2) )
#endif
      !
      istat = MAXVAL( ierr )
      IF( lk_mpp    )   CALL mpp_sum ( 'sbccpl_cesm', istat )
      IF( istat > 0 )   CALL ctl_stop('sbc_cpl_cesm_init: allocation of arrays failed')

      ! Initialization

      garea = glob_sum('sbccpl_cesm', e1e2t(:,:))

      taux_x2o   = 0.0_wp
      tauy_x2o   = 0.0_wp
      evap_x2o   = 0.0_wp
      rain_x2o   = 0.0_wp
      snow_x2o   = 0.0_wp
      roff_x2o   = 0.0_wp
      ioff_x2o   = 0.0_wp
      meltw_x2o  = 0.0_wp
      salt_x2o   = 0.0_wp
      swnet_x2o  = 0.0_wp
      sen_x2o    = 0.0_wp
      lat_x2o    = 0.0_wp
      lwup_x2o   = 0.0_wp
      lwdn_x2o   = 0.0_wp
      melth_x2o  = 0.0_wp
      ifrac_x2o  = 0.0_wp
      pslv_x2o   = 0.0_wp
      duu10n_x2o = 0.0_wp
      apr        = 0.0_wp
      utau       = 0.0_wp
      vtau       = 0.0_wp
      taum       = 0.0_wp
      emp        = 0.0_wp
      sfx        = 0.0_wp
      rnf        = 0.0_wp
      qsr        = 0.0_wp
      qns        = 0.0_wp
      fr_i       = 0.0_wp
      wndm       = 0.0_wp
      apr        = 0.0_wp
#if defined key_cpl_carbon_cycle
      co2_x2o    = 0.0_wp
      atm_co2    = 0.0_wp
#endif

   END SUBROUTINE sbc_cpl_cesm_init

   SUBROUTINE sbc_cpl_cesm_rcv( kt, k_fsbc, k_ice )     
      !! 
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
!
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lat    -- latent heat flux                         (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction
!    o  roff   -- river runoff flux                        (kg/m2/s)
!    o  ioff   -- ice runoff flux                          (kg/m2/s)
!
!    The following fields are sometimes received from the coupler,
!      depending on model options:
!
!    o  pslv   -- sea-level pressure                       (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
!
!-----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt       ! ocean model time step index
      INTEGER, INTENT(in) ::   k_fsbc   ! frequency of sbc (-> ice model) computation 
      INTEGER, INTENT(in) ::   k_ice    ! ice management in the sbc (=0/1/2/3)
      !!
      ! local data
      REAL(wp) :: zrnfex   ! excess runoff to be redistributed over E-P
      REAL(wp) :: tmpvar   ! excess runoff to be redistributed over E-P
      INTEGER  :: i, j     ! dummy loop index
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   ztx, zty, zcptn

      !IF (lwp) WRITE(numout,*) 'sbc_cpl_rcv: called ', kt, lrecv

      ! Return if this isn't a coupling time step
      IF (.NOT. lrecv) RETURN

      !  1. distribute wind stress
      !     rotate components to local coordinates
      !     shift from T points to U,V points

      taux_x2o(:,:) = taux_x2o(:,:)*tmask(:,:,1)
      tauy_x2o(:,:) = tauy_x2o(:,:)*tmask(:,:,1)

      ! rotate true zonal/meridional wind stress into local coordinates
      call rot_rep(taux_x2o(:,:), tauy_x2o(:,:), 'T', 'en->i', ztx(:,:))
      call rot_rep(taux_x2o(:,:), tauy_x2o(:,:), 'T', 'en->j', zty(:,:))

      ! LBC (halo) update
      call lbc_lnk( 'sbccpl_cesm', ztx(:,:), 'T', -1._wp, zty(:,:), 'T', -1._wp )

      ! and shift T to U,V grid
      utau(:,:) = 0.0_wp
      vtau(:,:) = 0.0_wp
      DO j = 1, jpj-1                                          ! T ==> (U,V)
         DO i = 1, jpi-1   ! vector opt.
            IF ((tmask(i+1,j,1)+tmask(i,j,1))/=0.0_wp) THEN
               utau(i,j) = (ztx(i+1,j)*tmask(i+1,j,1)+ztx(i,j)*tmask(i,j,1)) / &
                           (tmask(i+1,j,1)+tmask(i,j,1))
            ENDIF
            IF ((tmask(i,j+1,1)+tmask(i,j,1))/=0.0_wp) THEN
               vtau(i,j) = (zty(i,j+1)*tmask(i,j+1,1)+zty(i,j)*tmask(i,j,1)) / &
                           (tmask(i,j+1,1)+tmask(i,j,1))
            ENDIF
         END DO
      END DO
      utau(:,:) = utau(:,:)*umask(:,:,1)
      vtau(:,:) = vtau(:,:)*vmask(:,:,1)

      ! Compute wind stress module at T points
      taum(:,:) = SQRT( ztx(:,:)*ztx(:,:) + zty(:,:)*zty(:,:) )

!     2. distribute fresh water, salt and heat fluxes

      ! Water budget (E-P-R)
      ! NEMO: emp>0 ==> d(SSH)/dt<0 (evaporative water budget ==> ocean loses water)
      !       Sign convention: positive upward
      !       Unit system: MKS

      ! Put runoff in rnf so it's taken into account in sbc_rnf_div
      ! when ln_rnf = .TRUE.
      rnf(:,:) = (roff_x2o(:,:)+ioff_x2o(:,:)) * tmask(:,:,1)
      ! Set a threshold on the runoff (rn_rnf_bnd) and redistribute excess
      ! runoff on the global E-P to avoid SSS<0.0 near river mouths
      zrnfex = 0.0_wp
      IF ( rn_rnf_bnd > 0.0_wp ) THEN
         ! Compute excess runoff
         ztx(:,:) = max(rnf(:,:)-rn_rnf_bnd,0.0_wp)
         rnf(:,:) = min(rnf(:,:),rn_rnf_bnd)
         zrnfex   = glob_sum('sbccpl_cesm', ztx(:,:)*e1e2t(:,:))
      END IF

      IF ( ln_rnf ) THEN
         emp(:,:) = -(evap_x2o(:,:)+rain_x2o(:,:)+snow_x2o(:,:)+meltw_x2o(:,:))
      ELSE
!         emp(:,:) = emp(:,:)-rnf(:,:)
         emp(:,:) = -(evap_x2o(:,:)+rain_x2o(:,:)+snow_x2o(:,:)+meltw_x2o(:,:)+ &
                      roff_x2o(:,:)+ioff_x2o(:,:))
         rnf(:,:) = 0._wp
      ENDIF
      !
      IF ( zrnfex > 0.0_wp ) THEN
         ! Redistribute excess runoff on the global E-P
         zrnfex    = zrnfex/garea
         emp(:,:)  = emp(:,:) -zrnfex
!         IF (lwp) WRITE(numout,*) ' RNF REDIST ',zrnfex
      END IF

      ! Salt Flux over the ocean from cpl  [(Kg*psu)/(m^2 s)]
      sfx(:,:) = salt_x2o(:,:)
      ! Compute ice only water flux from virtual salt flux
      sicew_x2o(:,:)=0._wp
      WHERE (sss_m(:,:) .GT. 0._wp) sicew_x2o(:,:)= salt_x2o(:,:) / sss_m(:,:)

      ! Heat budget
      ! NEMO: [qsr, qns]>0 ==> d(SST)/dt>0
      !       Sign convention: positive downward
      !       Units: MKS
      ! in NEMO 3.6 heat flux related to emp has to be subtracted from qns
      zcptn(:,:) = rcp * sst_m(:,:)
      ! Solar radiation 
      qsr(:,:)  = swnet_x2o(:,:)
      ! Non solar radiation minus emp
      qns(:,:)  = lwup_x2o(:,:)+lwdn_x2o(:,:)+sen_x2o(:,:)+lat_x2o(:,:)+         &
                  melth_x2o(:,:)-(snow_x2o(:,:)+ioff_x2o(:,:))*rLfus              &
                  - emp(:,:) * zcptn(:,:)

      emp(:,:) = emp(:,:)*tmask(:,:,1)
      sfx(:,:) = sfx(:,:)*tmask(:,:,1)
      qsr(:,:) = qsr(:,:)*tmask(:,:,1)
      qns(:,:) = qns(:,:)*tmask(:,:,1)

!     3. distribute ice cover, wind speed module
      fr_i(:,:) = ifrac_x2o(:,:)*tmask(:,:,1)
      ! m2/s2 -> m/s
      wndm(:,:) = 0.0_wp
      WHERE (duu10n_x2o(:,:) > 0.0_wp)
        wndm(:,:) = SQRT(duu10n_x2o(:,:))
      END WHERE
      wndm(:,:) = wndm(:,:)*tmask(:,:,1)

!     4. distribute atmospheric pressure and convert from Pa to hPa
      apr(:,:) = pslv_x2o(:,:)*0.01*tmask(:,:,1)

#if defined key_cpl_carbon_cycle
!     5. distribute atmospheric co2
      atm_co2(:,:) = co2_x2o(:,:)*tmask(:,:,1)
#endif

!     update ghost cells for fluxes received from the coupler
      call lbc_lnk( 'sbccpl_cesm', utau(:,:), 'U', -1._wp, vtau(:,:), 'V', -1._wp, &
        &         apr(:,:) , 'T', 1._wp, emp(:,:), 'T', 1._wp, sfx(:,:) , 'T', 1._wp, &
        &         qsr(:,:) , 'T', 1._wp, qns(:,:), 'T', 1._wp, fr_i(:,:), 'T', 1._wp, &
        &         wndm(:,:), 'T', 1._wp)

      IF ( ln_rnf ) call lbc_lnk('sbccpl_cesm', rnf(:,:),  'T', 1._wp)

#if defined key_cpl_carbon_cycle
      call lbc_lnk('sbccpl_cesm', atm_co2(:,:), 'T', 1._wp)
#endif

      ! Reset coupling time step flag
      lrecv = .FALSE.

   END SUBROUTINE sbc_cpl_cesm_rcv

   SUBROUTINE sbc_cpl_cesm_finalize

   IF (ALLOCATED(taux_x2o))   DEALLOCATE(taux_x2o)
   IF (ALLOCATED(tauy_x2o))   DEALLOCATE(tauy_x2o)
   IF (ALLOCATED(evap_x2o))   DEALLOCATE(evap_x2o)
   IF (ALLOCATED(rain_x2o))   DEALLOCATE(rain_x2o)
   IF (ALLOCATED(snow_x2o))   DEALLOCATE(snow_x2o)
   IF (ALLOCATED(roff_x2o))   DEALLOCATE(roff_x2o)
   IF (ALLOCATED(ioff_x2o))   DEALLOCATE(ioff_x2o)
   IF (ALLOCATED(meltw_x2o))  DEALLOCATE(meltw_x2o)
   IF (ALLOCATED(salt_x2o))   DEALLOCATE(salt_x2o)
   IF (ALLOCATED(swnet_x2o))  DEALLOCATE(swnet_x2o)
   IF (ALLOCATED(lat_x2o))    DEALLOCATE(lat_x2o)
   IF (ALLOCATED(sen_x2o))    DEALLOCATE(sen_x2o)
   IF (ALLOCATED(lwup_x2o))   DEALLOCATE(lwup_x2o)
   IF (ALLOCATED(lwdn_x2o))   DEALLOCATE(lwdn_x2o)
   IF (ALLOCATED(melth_x2o))  DEALLOCATE(melth_x2o)
   IF (ALLOCATED(ifrac_x2o))  DEALLOCATE(ifrac_x2o)
   IF (ALLOCATED(pslv_x2o))   DEALLOCATE(pslv_x2o)
   IF (ALLOCATED(duu10n_x2o)) DEALLOCATE(duu10n_x2o)
#if defined key_cpl_carbon_cycle
   IF (ALLOCATED(co2_x2o))    DEALLOCATE(co2_x2o)
#endif
   IF (ALLOCATED(sicew_x2o))   DEALLOCATE(sicew_x2o)
   END SUBROUTINE sbc_cpl_cesm_finalize

#else
   IMPLICIT NONE
   PRIVATE

   LOGICAL, PARAMETER, PUBLIC :: lk_cesm = .false.
#endif

   !!======================================================================
END MODULE sbccpl_cesm
