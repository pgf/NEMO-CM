#if defined CCSMCOUPLED
MODULE nemogcm
   !!======================================================================
   !!                       ***  MODULE nemogcm   ***
   !! Ocean system   : NEMO GCM (ocean dynamics, on-line tracers, biochemistry and sea-ice)
   !!======================================================================
   !! History :  OPA  ! 1990-10  (C. Levy, G. Madec)  Original code
   !!            7.0  ! 1991-11  (M. Imbard, C. Levy, G. Madec)
   !!            7.1  ! 1993-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, C. Perigaud, G. Caniaux, B. Colot, C. Maes) release 7.1
   !!             -   ! 1992-06  (L.Terray)  coupling implementation
   !!             -   ! 1993-11  (M.A. Filiberti) IGLOO sea-ice
   !!            8.0  ! 1996-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, L.Terray, M.A. Filiberti, J. Vialar, A.M. Treguier, M. Levy) release 8.0
   !!            8.1  ! 1997-06  (M. Imbard, G. Madec)
   !!            8.2  ! 1999-11  (M. Imbard, H. Goosse)  sea-ice model
   !!                 ! 1999-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP
   !!                 ! 2000-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and modules
   !!             -   ! 2004-06  (R. Redler, NEC CCRLE, Germany) add OASIS[3/4] coupled interfaces
   !!             -   ! 2004-08  (C. Talandier) New trends organization
   !!             -   ! 2005-06  (C. Ethe) Add the 1D configuration possibility
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!             -   ! 2006-03  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   ! 2006-04  (G. Madec, R. Benshila)  Step reorganization
   !!             -   ! 2007-07  (J. Chanut, A. Sellar) Unstructured open boundaries (BDY)
   !!            3.2  ! 2009-08  (S. Masson)  open/write in the listing file in mpp
   !!            3.3  ! 2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!             -   ! 2011-11  (C. Harris) decomposition changes for running with CICE
   !!            3.6  ! 2012-05  (C. Calone, J. Simeon, G. Madec, C. Ethe) Add grid coarsening
   !!             -   ! 2014-12  (G. Madec) remove KPP scheme and cross-land advection (cla)
   !!            4.0  ! 2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   nemo_gcm      : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
   !!   nemo_init     : initialization of the NEMO system
   !!   nemo_ctl      : initialisation of the contol print
   !!   nemo_closefile: close remaining open files
   !!   nemo_alloc    : dynamical allocation
   !!----------------------------------------------------------------------
   USE step_oce       ! module used in the ocean time stepping module (step.F90)
   !
   USE phycst         ! physical constant                  (par_cst routine)
   USE domain         ! domain initialization   (dom_init & dom_cfg routines)
   USE wet_dry        ! Wetting and drying setting   (wad_init routine)
   USE usrdef_nam     ! user defined configuration namelist
   USE tide_mod, ONLY : tide_init   ! tidal components initialization   (tide_init routine)
   USE bdyini  , ONLY : bdy_init    ! open boundary cond. setting       (bdy_init routine)
   USE istate         ! initial state setting          (istate_init routine)
   USE trdini         ! dyn/tra trends initialization     (trd_init routine)
   USE icbini         ! handle bergs, initialisation
   USE icbstp  , ONLY : icb_end     ! handle bergs, close iceberg files
   USE cpl_oasis3     ! OASIS3 coupling
   USE dyndmp         ! Momentum damping (C1D only)
   USE step_diu       ! diurnal bulk SST timestepping (called from here if run offline)
   USE crsini         ! initialise grid coarsening utility
   USE dia25h  , ONLY : dia_25h_init   ! 25h mean output (initialisation)
   USE c1d            ! 1D configuration
#if defined key_top
   USE trcini         ! passive tracer initialisation
#endif
#if defined key_nemocice_decomp
   USE ice_domain_size, only: nx_global, ny_global
#endif
#if defined key_qco   ||   defined key_linssh
   USE stpmlf         ! NEMO time-stepping               (stp_MLF   routine)
#else
   USE step           ! NEMO time-stepping                 (stp     routine)
#endif
   !
   USE lib_mpp        ! distributed memory computing
   USE mppini         ! shared/distributed memory setting (mpp_init routine)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE halo_mng       ! halo manager
   USE sbccpl_cesm,   ONLY: lk_cesm

   USE shr_sys_mod,   ONLY: shr_sys_abort
#if defined NUOPC_IN_NEMO
   USE ESMF
   use nuopc_shr_methods , only : chkerr
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   nemo_gcm    ! called by model.F90
   PUBLIC   nemo_init   ! needed by AGRIF
   PUBLIC   nemo_alloc  ! needed by TAM
   PUBLIC   nemo_closefile   ! needed by ocn_comp_mct.F90
   PUBLIC   cform_aaa        ! needed by ocn_comp_mct.F90

   CHARACTER(lc) :: diri
   CHARACTER(lc) :: diro
   CHARACTER(lc),public :: logfile
   namelist / modelio / diri, diro, logfile

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

#if ! defined key_mpi_off
   ! need MPI_Wtime
   INCLUDE 'mpif.h'
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: nemogcm.F90 15267 2021-09-17 09:04:34Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_gcm
      !!----------------------------------------------------------------------
      !!                     ***  DUMMY ROUTINE nemo_gcm  ***
      !!
      !!----------------------------------------------------------------------
      write(*,*) 'nemo_gcm: You should not see this with CESM!!'

   END SUBROUTINE nemo_gcm


   SUBROUTINE nemo_init(ilocal_comm)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_init  ***
      !!
      !! ** Purpose :   initialization of the NEMO GCM
      !!----------------------------------------------------------------------
      INTEGER ::   ios                ! local integers
      INTEGER, INTENT(INOUT) ::   ilocal_comm   ! local MPI comm received from CPL driver
      !
      INTEGER :: inum          ! tmp Fortran unit number
      INTEGER :: rcode         ! return status
      LOGICAL :: exists        ! file existance logical
      CHARACTER(len=*), PARAMETER :: func = 'nemo_cesm_init'
      CHARACTER(len=80) :: nmlfile   ! modelio namelist file name

#if defined NUOPC_IN_NEMO
      type(ESMF_Config) :: config
      integer :: rc
      character(*),parameter :: u_FILE_u = &
       __FILE__
#endif
      !!
      NAMELIST/namctl/ sn_cfctl, ln_timing, ln_diacfl, nn_isplt, nn_jsplt , nn_ictls,   &
         &                                             nn_ictle, nn_jctls , nn_jctle
      NAMELIST/namcfg/ ln_read_cfg, cn_domcfg, ln_closea, ln_write_cfg, cn_domcfg_out, ln_use_jattr
      !!----------------------------------------------------------------------
      !
      cxios_context = 'nemo'
      !
      !                             !-------------------------------------------------!
      !                             !     set communicator & select the local rank    !
      !                             !  must be done as soon as possible to get narea  !
      !                             !-------------------------------------------------!
      !
#if defined key_xios
      IF( Agrif_Root() ) THEN
         IF( lk_oasis ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )                               ! nemo local communicator given by oasis
            CALL xios_initialize( "not used"       , local_comm =ilocal_comm )   ! send nemo communicator to xios
         ELSEIF ( lk_cesm ) THEN
            CALL xios_initialize( "not used",local_comm=ilocal_comm )    ! send nemo communicator given by cesm to xios

         ELSE
            CALL xios_initialize( "for_xios_mpi_id", return_comm=ilocal_comm )   ! nemo local communicator given by xios
         ENDIF
      ENDIF
      CALL mpp_start( ilocal_comm )
#else
      IF( lk_oasis ) THEN
         IF( Agrif_Root() ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )          ! nemo local communicator given by oasis
         ENDIF
         CALL mpp_start( ilocal_comm )
      ELSE
         CALL mpp_start( )
      ENDIF
#endif
      !
      narea = mpprank + 1               ! mpprank: the rank of proc (0 --> mppsize -1 )
      lwm = (narea == 1)                ! control of output namelists
      !
      !                             !---------------------------------------------------------------!
      !                             ! Open output files, reference and configuration namelist files !
      !                             !---------------------------------------------------------------!
      !
      ! open ocean.output as soon as possible to get all output prints (including errors messages)
      IF( lwm )   CALL ctl_opn(     numout,        'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open reference and configuration namelist files
                  CALL load_nml( numnam_ref,        'namelist_ref',                                           -1, lwm )
                  CALL load_nml( numnam_cfg,        'namelist_cfg',                                           -1, lwm )
      IF( lwm )   CALL ctl_opn(     numond, 'output.namelist.dyn', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open /dev/null file to be able to supress output write easily
      IF( Agrif_Root() ) THEN
                  CALL ctl_opn(     numnul,           '/dev/null', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
#ifdef key_agrif
      ELSE
                  numnul = Agrif_Parent(numnul)
#endif
      ENDIF
      !                             !--------------------!
      !                             ! Open listing units !  -> need sn_cfctl from namctl to define lwp
      !                             !--------------------!
      !
      READ  ( numnam_ref, namctl, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namctl in reference namelist' )
      READ  ( numnam_cfg, namctl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namctl in configuration namelist' )
      !
      ! finalize the definition of namctl variables
      IF( narea < sn_cfctl%procmin .OR. narea > sn_cfctl%procmax .OR. MOD( narea - sn_cfctl%procmin, sn_cfctl%procincr ) /= 0 )   &
         &   CALL nemo_set_cfctl( sn_cfctl, .FALSE. )
      !
      lwp = (narea == 1) .OR. sn_cfctl%l_oceout    ! control of all listing output print
      !
      diri = '.'
      diro = '.'
      logfile = ''

#if defined NUOPC_IN_NEMO
      nmlfile = 'nuopc.runconfig'   ! TODO: multi-instance version
      INQUIRE(FILE=TRIM(nmlfile),EXIST=exists)

      IF (.NOT. exists) THEN
         ! Fall back to the default ocean.output file
         IF(lwp) THEN                      ! open listing units
            !
            IF( .NOT. lwm )   &            ! alreay opened for narea == 1
               &            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE., narea )
         ENDIF
      ELSE
         ! Open CESM style log file (ocn.log.*)
         config = ESMF_ConfigCreate(rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         call ESMF_ConfigLoadFile(config, "nuopc.runconfig", rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         
         call ESMF_ConfigGetAttribute(config,value=diro,label="diro = ", rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         call ESMF_ConfigGetAttribute(config,value=logfile,label="logfile = ocn.", rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         IF(lwp) THEN                            ! open listing units
         !
            IF (LEN_TRIM(logfile) > 0) THEN
               numout = get_unit()
               OPEN(numout,FILE=TRIM(diro)//'/ocn.'//TRIM(logfile),STATUS='REPLACE', &
                   ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=rcode)
               IF (rcode /= 0) THEN
                  WRITE(*,FMT='(3A,I6)') 'ERROR: opening ',TRIM(logfile),': iostat=',rcode
                  CALL shr_sys_abort(func//': ERROR opening '//TRIM(logfile) )
               ENDIF
            ELSE
               ! Fall back to the default ocean.output file
               ! WRITE(numout,FMT='(A)') 'logfile not opened'
               CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', &
                   'SEQUENTIAL', -1, 6, .FALSE., narea )
            ENDIF
         ENDIF
      ENDIF
#else
      nmlfile = 'ocn_modelio.nml'   ! TODO: multi-instance version
      INQUIRE(FILE=TRIM(nmlfile),EXIST=exists)

      IF (.NOT. exists) THEN
         ! Fall back to the default ocean.output file
         IF(lwp) THEN                      ! open listing units
            !
            IF( .NOT. lwm )   &            ! alreay opened for narea == 1
               &            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE., narea )
         ENDIF
      ELSE
         ! Open CESM style log file (ocn.log.*)
         inum = get_unit()
         rcode = 0
         !
         OPEN (inum,FILE=nmlfile,ACTION='READ',IOSTAT=rcode)
         IF (rcode /= 0) THEN
            WRITE(*,FMT='(3A,I6)') 'ERROR: opening ',TRIM(nmlfile),': iostat=',rcode
            CALL shr_sys_abort(func//': ERROR opening '//TRIM(nmlfile) )
         ENDIF
         !
         READ (inum,NML=modelio,IOSTAT=rcode)
         IF (rcode /= 0) THEN
            WRITE(*,FMT='(3A,I6)') 'ERROR: reading ',TRIM(nmlfile),': iostat=',rcode
            CALL shr_sys_abort(func//': ERROR reading '//TRIM(nmlfile) )
         ENDIF
         !
         CLOSE(inum)
         !
         IF(lwp) THEN                            ! open listing units
         !
            IF (LEN_TRIM(logfile) > 0) THEN
               IF (sn_cfctl%l_prtctl) THEN
                 WRITE(logfile,FMT='(A,"_",I4.4)') TRIM(logfile), narea-1
               ENDIF
               numout = get_unit()
               OPEN(numout,FILE=TRIM(diro)//'/'//TRIM(logfile),STATUS='REPLACE', &
                   ACCESS='SEQUENTIAL',FORM='FORMATTED',IOSTAT=rcode)
               IF (rcode /= 0) THEN
                  WRITE(*,FMT='(3A,I6)') 'ERROR: opening ',TRIM(logfile),': iostat=',rcode
                  CALL shr_sys_abort(func//': ERROR opening '//TRIM(logfile) )
               ENDIF
            ELSE
               ! Fall back to the default ocean.output file
               ! WRITE(numout,FMT='(A)') 'logfile not opened'
               CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', &
                   'SEQUENTIAL', -1, 6, .FALSE., narea )
            ENDIF
         ENDIF
      ENDIF
#endif
      !
      IF(lwp) THEN                      ! open listing units
         !
         IF( .NOT. lwm )   &            ! alreay opened for narea == 1
            &            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE., narea )
         !
         WRITE(numout,*)
         WRITE(numout,*) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC'
         WRITE(numout,*) '                       NEMO team'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '                NEMO version 4.0  (2020) '
         WRITE(numout,*)
         WRITE(numout,*) "           ._      ._      ._      ._      ._    "
         WRITE(numout,*) "       _.-._)`\_.-._)`\_.-._)`\_.-._)`\_.-._)`\_ "
         WRITE(numout,*)
         WRITE(numout,*) "           o         _,           _,             "
         WRITE(numout,*) "            o      .' (        .-' /             "
         WRITE(numout,*) "           o     _/..._'.    .'   /              "
         WRITE(numout,*) "      (    o .-'`      ` '-./  _.'               "
         WRITE(numout,*) "       )    ( o)           ;= <_         (       "
         WRITE(numout,*) "      (      '-.,\\__ __.-;`\   '.        )      "
         WRITE(numout,*) "       )  )       \) |`\ \)  '.   \      (   (   "
         WRITE(numout,*) "      (  (           \_/       '-._\      )   )  "
         WRITE(numout,*) "       )  ) jgs                     `    (   (   "
         WRITE(numout,*) "     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "
         WRITE(numout,*)
         WRITE(numout,*) "     NEMO is running in the NCAR CESM framework  "
         WRITE(numout,*)
         !
         WRITE(numout,cform_aaa)    ! Flag AAAAAAA
         !
         !                          ! Control print of the working precision
         WRITE(numout,*)
         IF( wp == dp ) THEN   ;   WRITE(numout,*) "par_kind : wp = Working precision = dp = double-precision"
         ELSE                  ;   WRITE(numout,*) "par_kind : wp = Working precision = sp = single-precision"
         ENDIF
                                   WRITE(numout,*) "~~~~~~~~                                 ****************"
                                   WRITE(numout,*)
         !
      ENDIF
      !
      IF(lwm) WRITE( numond, namctl )
      !
      !                             !------------------------------------!
      !                             !  Set global domain size parameters !
      !                             !------------------------------------!
      !
      READ  ( numnam_ref, namcfg, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcfg in reference namelist' )
      READ  ( numnam_cfg, namcfg, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcfg in configuration namelist' )
      !
      IF( ln_read_cfg ) THEN            ! Read sizes in domain configuration file
         CALL domain_cfg ( cn_cfg, nn_cfg, Ni0glo, Nj0glo, jpkglo, l_Iperio, l_Jperio, l_NFold, c_NFtype )
      ELSE                              ! user-defined namelist
         CALL usr_def_nam( cn_cfg, nn_cfg, Ni0glo, Nj0glo, jpkglo, l_Iperio, l_Jperio, l_NFold, c_NFtype )
      ENDIF
      !
      IF(lwm)   WRITE( numond, namcfg )
      !
      !                             !-----------------------------------------!
      !                             ! mpp parameters and domain decomposition !
      !                             !-----------------------------------------!
      CALL mpp_init

#if defined key_loop_fusion
      IF( nn_hls == 1 ) THEN
         CALL ctl_stop( 'STOP', 'nemogcm : Loop fusion can be used only with extra-halo' )
      ENDIF
#endif

      CALL halo_mng_init()
      ! Now we know the dimensions of the grid and numout has been set: we can allocate arrays
      CALL nemo_alloc()

      ! Initialise time level indices
      Nbb = 1   ;   Nnn = 2   ;   Naa = 3   ;   Nrhs = Naa
#if defined key_agrif
      Kbb_a = Nbb   ;   Kmm_a = Nnn   ;   Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
#endif
      !                             !-------------------------------!
      !                             !  NEMO general initialization  !
      !                             !-------------------------------!

      CALL nemo_ctl                          ! Control prints of namctl and namcfg
      !
      !                                      ! General initialization
      IF( ln_timing    )   CALL timing_init     ! timing
      IF( ln_timing    )   CALL timing_start( 'nemo_init')
      !
                           CALL     phy_cst         ! Physical constants
                           CALL     eos_init        ! Equation of state
                           CALL     wad_init        ! Wetting and drying options

#if defined key_agrif
     CALL Agrif_Declare_Var_ini   !  "      "   "   "      "  DOM
#endif
                           CALL     dom_init( Nbb, Nnn, Naa )   ! Domain
      IF( ln_crs       )   CALL     crs_init(      Nnn      )   ! coarsened grid: domain initialization
      IF( sn_cfctl%l_prtctl )   &
         &                 CALL prt_ctl_init        ! Print control

                           CALL diurnal_sst_bulk_init       ! diurnal sst
      IF( ln_diurnal   )   CALL diurnal_sst_coolskin_init   ! cool skin
      !
      IF( ln_diurnal_only ) THEN                    ! diurnal only: a subset of the initialisation routines
         CALL  istate_init( Nbb, Nnn, Naa )         ! ocean initial state (Dynamics and tracers)
         CALL     sbc_init( Nbb, Nnn, Naa )         ! Forcings : surface module
         CALL tra_qsr_init                          ! penetrative solar radiation qsr
         IF( ln_diaobs ) THEN                       ! Observation & model comparison
            CALL dia_obs_init( Nnn )                ! Initialize observational data
            CALL dia_obs( nit000 - 1, Nnn )         ! Observation operator for restart
         ENDIF
         IF( lk_asminc )   CALL asm_inc_init( Nbb, Nnn, Nrhs )   ! Assimilation increments
         !
         RETURN                                       ! end of initialization
      ENDIF
      !

                           CALL  istate_init( Nbb, Nnn, Naa )    ! ocean initial state (Dynamics and tracers)

      !                                      ! external forcing
                           CALL    tide_init                     ! tidal harmonics
                           CALL     sbc_init( Nbb, Nnn, Naa )    ! surface boundary conditions (including sea-ice)
                           CALL     bdy_init                     ! Open boundaries initialisation

      !                                      ! Ocean physics
                           CALL zdf_phy_init( Nnn )    ! Vertical physics

      !                                         ! Lateral physics
                           CALL ldf_tra_init      ! Lateral ocean tracer physics
                           CALL ldf_eiv_init      ! eddy induced velocity param. must be done after ldf_tra_init
                           CALL ldf_dyn_init      ! Lateral ocean momentum physics

      !                                      ! Active tracers
      IF( ln_traqsr    )   CALL tra_qsr_init      ! penetrative solar radiation qsr
                           CALL tra_bbc_init      ! bottom heat flux
                           CALL tra_bbl_init      ! advective (and/or diffusive) bottom boundary layer scheme
                           CALL tra_dmp_init      ! internal tracer damping
                           CALL tra_adv_init      ! horizontal & vertical advection
                           CALL tra_ldf_init      ! lateral mixing

      !                                      ! Dynamics
      IF( ln_c1d       )   CALL dyn_dmp_init         ! internal momentum damping
                           CALL dyn_adv_init         ! advection (vector or flux form)
                           CALL dyn_vor_init         ! vorticity term including Coriolis
                           CALL dyn_ldf_init         ! lateral mixing
                           CALL dyn_hpg_init( Nnn )  ! horizontal gradient of Hydrostatic pressure
                           CALL dyn_spg_init         ! surface pressure gradient

      !                                      ! Icebergs
                           CALL icb_init( rn_Dt, nit000)   ! initialise icebergs instance

                                                ! ice shelf
                           CALL isf_init( Nbb, Nnn, Naa )
#if defined key_top
      !                                      ! Passive tracers
                           CALL     trc_init( Nbb, Nnn, Naa )
#endif
      IF( l_ldfslp     )   CALL ldf_slp_init    ! slope of lateral mixing

      !                                      ! Misc. options
                           CALL sto_par_init    ! Stochastic parametrization
      IF( ln_sto_eos   )   CALL sto_pts_init    ! RRandom T/S fluctuations

      !                                      ! Diagnostics
                           CALL     flo_init( Nnn )    ! drifting Floats
      IF( ln_diacfl    )   CALL dia_cfl_init    ! Initialise CFL diagnostics
                           CALL dia_dct_init    ! Sections tranports
                           CALL dia_hsb_init( Nnn )    ! heat content, salt content and volume budgets
                           CALL     trd_init( Nnn )    ! Mixed-layer/Vorticity/Integral constraints trends
                           CALL dia_obs_init( Nnn )    ! Initialize observational data
                           CALL dia_25h_init( Nbb )    ! 25h mean  outputs
                           CALL dia_detide_init ! Weights computation for daily detiding of model diagnostics
      IF( ln_diaobs    )   CALL dia_obs( nit000-1, Nnn )   ! Observation operator for restart
                           CALL dia_mlr_init    ! Initialisation of IOM context management for multiple-linear-regression analysis

      !                                      ! Assimilation increments
      IF( lk_asminc    )   CALL asm_inc_init( Nbb, Nnn, Nrhs )   ! Initialize assimilation increments
      !
      IF(lwp) WRITE(numout,cform_aaa)           ! Flag AAAAAAA
      !
      IF( ln_timing    )   CALL timing_stop( 'nemo_init')
      !
   END SUBROUTINE nemo_init


   SUBROUTINE nemo_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_ctl  ***
      !!
      !! ** Purpose :   control print setting
      !!
      !! ** Method  : - print namctl and namcfg information and check some consistencies
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_ctl: Control prints'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '                              sn_cfctl%l_runstat = ', sn_cfctl%l_runstat
         WRITE(numout,*) '                              sn_cfctl%l_trcstat = ', sn_cfctl%l_trcstat
         WRITE(numout,*) '                              sn_cfctl%l_oceout  = ', sn_cfctl%l_oceout
         WRITE(numout,*) '                              sn_cfctl%l_layout  = ', sn_cfctl%l_layout
         WRITE(numout,*) '                              sn_cfctl%l_prtctl  = ', sn_cfctl%l_prtctl
         WRITE(numout,*) '                              sn_cfctl%l_prttrc  = ', sn_cfctl%l_prttrc
         WRITE(numout,*) '                              sn_cfctl%l_oasout  = ', sn_cfctl%l_oasout
         WRITE(numout,*) '                              sn_cfctl%procmin   = ', sn_cfctl%procmin
         WRITE(numout,*) '                              sn_cfctl%procmax   = ', sn_cfctl%procmax
         WRITE(numout,*) '                              sn_cfctl%procincr  = ', sn_cfctl%procincr
         WRITE(numout,*) '                              sn_cfctl%ptimincr  = ', sn_cfctl%ptimincr
         WRITE(numout,*) '      timing by routine               ln_timing  = ', ln_timing
         WRITE(numout,*) '      CFL diagnostics                 ln_diacfl  = ', ln_diacfl
      ENDIF
      !
      IF( .NOT.ln_read_cfg )   ln_closea = .false.   ! dealing possible only with a domcfg file
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      read domain configuration file                ln_read_cfg      = ', ln_read_cfg
         WRITE(numout,*) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
         WRITE(numout,*) '         keep closed seas in the domain (if exist)     ln_closea     = ', ln_closea
         WRITE(numout,*) '      create a configuration definition file        ln_write_cfg     = ', ln_write_cfg
         WRITE(numout,*) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
         WRITE(numout,*) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
      ENDIF
#if defined key_xios
      IF(lwp) THEN                  ! print I/O setting
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_init  : key_xios in use'
         WRITE(numout,*) 'Nemo I/O is handled with XIOS 2.5'
         WRITE(numout,*)
      ENDIF
#endif
      !
      IF( 1._wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nemo_ctl: The intrinsec SIGN function follows f2003 standard.',  &
         &                                                'Compile with key_nosignedzero enabled:',   &
         &                                                '--> add -Dkey_nosignedzero to the definition of %CPP in your arch file' )
      !
#if defined key_agrif
      IF( ln_timing )   CALL ctl_stop( 'AGRIF not implemented with ln_timing = true')
#endif
      !
   END SUBROUTINE nemo_ctl


   SUBROUTINE nemo_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!----------------------------------------------------------------------
      !
      IF( lk_mpp )   CALL mppsync
      !
      CALL iom_close                                 ! close all input/output files managed by iom_*
      !
      IF( numstp          /= -1 )   CLOSE( numstp          )   ! time-step file
      IF( numrun          /= -1 )   CLOSE( numrun          )   ! run statistics file
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond          )   ! oce output namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni          )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice      )   ! ice variables (temp. evolution)
      IF( numdct_vol      /= -1 )   CLOSE( numdct_vol      )   ! volume transports
      IF( numdct_heat     /= -1 )   CLOSE( numdct_heat     )   ! heat transports
      IF( numdct_salt     /= -1 )   CLOSE( numdct_salt     )   ! salt transports
      !
      IF (lwp) write(numout,*) ' '
      IF (lwp) write(numout,*) '         ******* END OF MODEL RUN *******'
      IF (lwp) write(numout,*) ' '
      !
      IF( numout          /= -1 )   CLOSE( numout          )   ! standard model output file
      !
      !numout = 6                                     ! redefine numout in case it is used after this point...
      !
   END SUBROUTINE nemo_closefile


   SUBROUTINE nemo_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OCE modules
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE diawri    , ONLY : dia_wri_alloc
      USE dom_oce   , ONLY : dom_oce_alloc
      USE trc_oce   , ONLY : trc_oce_alloc
      USE bdy_oce   , ONLY : bdy_oce_alloc
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        oce_alloc    ()    ! ocean
      ierr = ierr + dia_wri_alloc()
      ierr = ierr + dom_oce_alloc()    ! ocean domain
      ierr = ierr + zdf_oce_alloc()    ! ocean vertical physics
      ierr = ierr + trc_oce_alloc()    ! shared TRC / TRA arrays
      ierr = ierr + bdy_oce_alloc()    ! bdy masks (incl. initialization)
      !
      CALL mpp_sum( 'nemogcm', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nemo_alloc: unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE nemo_alloc


   SUBROUTINE nemo_set_cfctl(sn_cfctl, setto )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_set_cfctl  ***
      !!
      !! ** Purpose :   Set elements of the output control structure to setto.
      !!
      !! ** Method  :   Note this routine can be used to switch on/off some
      !!                types of output for selected areas.
      !!----------------------------------------------------------------------
      TYPE(sn_ctl), INTENT(inout) :: sn_cfctl
      LOGICAL     , INTENT(in   ) :: setto
      !!----------------------------------------------------------------------
      sn_cfctl%l_runstat = setto
      sn_cfctl%l_trcstat = setto
      sn_cfctl%l_oceout  = setto
      sn_cfctl%l_layout  = setto
      sn_cfctl%l_prtctl  = setto
      sn_cfctl%l_prttrc  = setto
      sn_cfctl%l_oasout  = setto
   END SUBROUTINE nemo_set_cfctl

   !!======================================================================
END MODULE nemogcm
#endif
