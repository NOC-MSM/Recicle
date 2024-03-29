MODULE dtatsd
   !!======================================================================
   !!                     ***  MODULE  dtatsd  ***
   !! Ocean data  :  read ocean Temperature & Salinity Data from gridded data
   !!======================================================================
   !! History :  OPA  ! 1991-03  ()  Original code
   !!             -   ! 1992-07  (M. Imbard)
   !!            8.0  ! 1999-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module 
   !!            3.3  ! 2010-10  (C. Bricaud, S. Masson)  use of fldread
   !!            3.4  ! 2010-11  (G. Madec, C. Ethe) Merge of dtatem and dtasal + suppression of CPP keys
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dta_tsd      : read and time interpolated ocean Temperature & Salinity Data
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_tsd_init   ! called by opa.F90
   PUBLIC   dta_tsd        ! called by istate.F90 and tradmp.90

   LOGICAL , PUBLIC ::   ln_tsd_init      !: T & S data flag
   LOGICAL , PUBLIC ::   ln_tsd_interp    !: vertical interpolation flag
   LOGICAL , PUBLIC ::   ln_tsd_tradmp    !: internal damping toward input data flag

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tsd   ! structure of input SST (file informations, fields read)
   INTEGER                                 ::   jpk_init , inum_dta
   INTEGER                                 ::   id ,linum   ! local integers
   INTEGER                                 ::   zdim(4)

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dtatsd.F90 7753 2017-03-03 11:46:59Z mocavero $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_tsd_init( ld_tradmp )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd_init  ***
      !!                    
      !! ** Purpose :   initialisation of T & S input data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates T & S data structure 
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in), OPTIONAL ::   ld_tradmp   ! force the initialization when tradp is used
      !
      INTEGER ::   ios, ierr0, ierr1, ierr2, ierr3, ierr4, ierr5   ! local integers
      !!
      CHARACTER(len=100)            ::   cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(jpts+2)::   slf_i           ! array of namelist informations on the fields to read
      TYPE(FLD_N)                   ::   sn_tem, sn_sal, sn_dep, sn_msk
      
      !!
      NAMELIST/namtsd/   ln_tsd_init, ln_tsd_interp, ln_tsd_tradmp, cn_dir, sn_tem, sn_sal, sn_dep, sn_msk
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_tsd_init')
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0  ; ierr4 = 0  ;  ierr5 = 0 
      !
      REWIND( numnam_ref )              ! Namelist namtsd in reference namelist : 
      READ  ( numnam_ref, namtsd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namtsd in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namtsd, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtsd )

      IF( PRESENT( ld_tradmp ) )   ln_tsd_tradmp = .TRUE.     ! forces the initialization when tradmp is used
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dta_tsd_init : Temperature & Salinity data '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtsd'
         WRITE(numout,*) '      Initialisation of ocean T & S with T &S input data   ln_tsd_init   = ', ln_tsd_init
         WRITE(numout,*) '      iInterpolation of initial conditions in the vertical ln_tsd_interp = ', ln_tsd_interp
         WRITE(numout,*) '      damping of ocean T & S toward T &S input data        ln_tsd_tradmp = ', ln_tsd_tradmp
         WRITE(numout,*)
         IF( .NOT.ln_tsd_init .AND. .NOT.ln_tsd_tradmp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   T & S data not used'
         ENDIF
      ENDIF
      !
      IF( ln_rstart .AND. ln_tsd_init ) THEN
         CALL ctl_warn( 'dta_tsd_init: ocean restart and T & S data intialisation, ',   &
            &           'we keep the restart T & S values and set ln_tsd_init to FALSE' )
         ln_tsd_init = .FALSE.
      ENDIF
      IF( ln_tsd_interp .AND. ln_tsd_tradmp ) THEN
            CALL ctl_stop( 'dta_tsd_init: Tracer damping and vertical interpolation not yet configured' )   ;   RETURN
      ENDIF
      IF( ln_tsd_interp .AND. LEN(TRIM(sn_msk%wname)) > 0 ) THEN
            CALL ctl_stop( 'dta_tsd_init: Using vertical interpolation and weights files not recommended' )   ;   RETURN
      ENDIF
      !
      !                             ! allocate the arrays (if necessary)
      IF(  ln_tsd_init .OR. ln_tsd_tradmp  ) THEN
         !
         IF( ln_tsd_interp ) THEN
           ALLOCATE( sf_tsd(jpts+2), STAT=ierr0 ) ! to carry the addtional depth information
         ELSE
           ALLOCATE( sf_tsd(jpts  ), STAT=ierr0 ) 
         ENDIF 
         IF( ierr0 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd_init: unable to allocate sf_tsd structure' )   ;   RETURN
         ENDIF
         !
         slf_i(jp_tem) = sn_tem   ;   slf_i(jp_sal) = sn_sal
         IF( ln_tsd_interp ) slf_i(jp_dep) = sn_dep   ;   slf_i(jp_msk) = sn_msk
         CALL fld_fill( sf_tsd, slf_i, cn_dir, 'dta_tsd', 'Temperature & Salinity data', 'namtsd' )

         IF( ln_tsd_interp ) THEN
            !CALL iom_open ( trim(cn_dir) // trim(sn_dep%clname), inum_dta ) 
            CALL fld_clopn ( sf_tsd(jp_dep) ) 
            IF(lwp) WRITE(numout,*) 'INFO: ', sf_tsd(jp_dep)%num, sn_dep%clvar
            id = iom_varid( sf_tsd(jp_dep)%num, sn_dep%clvar, zdim )
            jpk_init = zdim(3)
            IF(lwp) WRITE(numout,*) 'Dimension of veritcal coordinate in ICs: ', jpk_init
            !CALL iom_close( inum_dta )   ! Close the input file
            !
                                 ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr0 )
            IF( sn_tem%ln_tint ) ALLOCATE( sf_tsd(jp_tem)%fdta(jpi,jpj,jpk_init,2) , STAT=ierr1 )
                                 ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr2 )
            IF( sn_sal%ln_tint ) ALLOCATE( sf_tsd(jp_sal)%fdta(jpi,jpj,jpk_init,2) , STAT=ierr3 )  
                                 ALLOCATE( sf_tsd(jp_dep)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr4 )
                                 ALLOCATE( sf_tsd(jp_msk)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr5 )
         ELSE
                                 ALLOCATE( sf_tsd(jp_tem)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
            IF( sn_tem%ln_tint ) ALLOCATE( sf_tsd(jp_tem)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                 ALLOCATE( sf_tsd(jp_sal)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
            IF( sn_sal%ln_tint ) ALLOCATE( sf_tsd(jp_sal)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )  
         ENDIF ! ln_tsd_interp

         !
         IF( ierr0 + ierr1 + ierr2 + ierr3 + ierr4 + ierr5 > 0 ) THEN
            CALL ctl_stop( 'dta_tsd : unable to allocate T & S data arrays' )   ;   RETURN
         ENDIF
         !                         ! fill sf_tsd with sn_tem & sn_sal and control print
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_tsd_init')
      !
   END SUBROUTINE dta_tsd_init


   SUBROUTINE dta_tsd( kt, ptsd )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_tsd  ***
      !!                    
      !! ** Purpose :   provides T and S data at kt
      !! 
      !! ** Method  : - call fldread routine
      !!              - ORCA_R2: add some hand made alteration to read data  
      !!              - 'key_orca_lev10' interpolates on 10 times more levels
      !!              - s- or mixed z-s coordinate: vertical interpolation on model mesh
      !!              - ln_tsd_tradmp=F: deallocates the T-S data structure
      !!                as T-S data are no are used
      !!
      !! ** Action  :   ptsd   T-S data on medl mesh and interpolated at time-step kt
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt     ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   ptsd   ! T & S data
      !
      INTEGER ::   ji, jj, jk, jl, jk_init   ! dummy loop indicies
      INTEGER ::   ik, il0, il1, ii0, ii1, ij0, ij1        ! local integers
      REAL(wp)::   zl, zi
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_tsd')
      !
      CALL fld_read( kt, 1, sf_tsd )      !==   read T & S data at kt time step   ==!
      !
      !
!!gm  This should be removed from the code   ===>>>>  T & S files has to be changed
      !
      !                                   !==   ORCA_R2 configuration and T & S damping   ==! 
!!gm end
      !
      IF( kt == nit000 .AND. lwp )THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dta_tsd: interpolates T & S data onto current mesh'
      ENDIF
      !
      IF( ln_tsd_interp ) THEN                 ! probably should use pointers in the following to make more readable
      !
         DO jk = 1, jpk                        ! determines the intepolated T-S profiles at each (i,j) points
            DO jj= 1, jpj
               DO ji= 1, jpi
                  zl = gdept_0(ji,jj,jk)
                  IF( zl < sf_tsd(jp_dep)%fnow(ji,jj,1) ) THEN                     ! above the first level of data
                     ptsd(ji,jj,jk,jp_tem) = sf_tsd(jp_tem)%fnow(ji,jj,1) 
                     ptsd(ji,jj,jk,jp_sal) = sf_tsd(jp_sal)%fnow(ji,jj,1)
                  ELSEIF( zl > sf_tsd(jp_dep)%fnow(ji,jj,jpk_init) ) THEN          ! below the last level of data
                     ptsd(ji,jj,jk,jp_tem) = sf_tsd(jp_tem)%fnow(ji,jj,jpk_init)
                     ptsd(ji,jj,jk,jp_sal) = sf_tsd(jp_sal)%fnow(ji,jj,jpk_init)
                  ELSE                                                             ! inbetween : vertical interpolation between jk_init & jk_init+1
                     DO jk_init = 1, jpk_init-1                                    ! when  gdept(jk_init) < zl < gdept(jk_init+1)
                        IF( sf_tsd(jp_msk)%fnow(ji,jj,jk_init+1) == 0 ) THEN       ! if there is no data fill down
                           sf_tsd(jp_tem)%fnow(ji,jj,jk_init+1) = sf_tsd(jp_tem)%fnow(ji,jj,jk_init)
                           sf_tsd(jp_sal)%fnow(ji,jj,jk_init+1) = sf_tsd(jp_sal)%fnow(ji,jj,jk_init)
                        ENDIF
                        IF( (zl-sf_tsd(jp_dep)%fnow(ji,jj,jk_init)) * (zl-sf_tsd(jp_dep)%fnow(ji,jj,jk_init+1)) <= 0._wp ) THEN
                           zi = ( zl - sf_tsd(jp_dep)%fnow(ji,jj,jk_init) ) / &
                        &       (sf_tsd(jp_dep)%fnow(ji,jj,jk_init+1)-sf_tsd(jp_dep)%fnow(ji,jj,jk_init))
                           ptsd(ji,jj,jk,jp_tem) = sf_tsd(jp_tem)%fnow(ji,jj,jk_init) + &
                        &                          (sf_tsd(jp_tem)%fnow(ji,jj,jk_init+1)-sf_tsd(jp_tem)%fnow(ji,jj,jk_init)) * zi
                           ptsd(ji,jj,jk,jp_sal) = sf_tsd(jp_sal)%fnow(ji,jj,jk_init) + &
                        &                          (sf_tsd(jp_sal)%fnow(ji,jj,jk_init+1)-sf_tsd(jp_sal)%fnow(ji,jj,jk_init)) * zi
                        ENDIF
                     END DO
                  ENDIF
               ENDDO
            ENDDO
         END DO
         !
         ptsd(:,:,:,jp_tem) = ptsd(:,:,:,jp_tem) *tmask(:,:,:)
         ptsd(:,:,:,jp_sal) = ptsd(:,:,:,jp_sal) *tmask(:,:,:)
      ELSE                                !==   z- or zps- coordinate   ==!
         !                             
         ptsd(:,:,:,jp_tem) = sf_tsd(jp_tem)%fnow(:,:,:)  * tmask(:,:,:)  ! Mask
         ptsd(:,:,:,jp_sal) = sf_tsd(jp_sal)%fnow(:,:,:)  * tmask(:,:,:)
         !
         IF( ln_zps ) THEN                      ! zps-coordinate (partial steps) interpolation at the last ocean level
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ik = mbkt(ji,jj) 
                  IF( ik > 1 ) THEN
                     zl = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                     ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik-1,jp_tem)
                     ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik-1,jp_sal)
                  ENDIF
                  ik = mikt(ji,jj)
                  IF( ik > 1 ) THEN
                     zl = ( gdept_0(ji,jj,ik) - gdept_1d(ik) ) / ( gdept_1d(ik+1) - gdept_1d(ik) ) 
                     ptsd(ji,jj,ik,jp_tem) = (1.-zl) * ptsd(ji,jj,ik,jp_tem) + zl * ptsd(ji,jj,ik+1,jp_tem)
                     ptsd(ji,jj,ik,jp_sal) = (1.-zl) * ptsd(ji,jj,ik,jp_sal) + zl * ptsd(ji,jj,ik+1,jp_sal)
                  END IF
               END DO
            END DO
         ENDIF
         !
      ENDIF
      !
      IF( .NOT.ln_tsd_tradmp ) THEN                   !==   deallocate T & S structure   ==! 
         !                                              (data used only for initialisation)
         IF(lwp) WRITE(numout,*) 'dta_tsd: deallocte T & S arrays as they are only use to initialize the run'
                                        DEALLOCATE( sf_tsd(jp_tem)%fnow )     ! T arrays in the structure
         IF( sf_tsd(jp_tem)%ln_tint )   DEALLOCATE( sf_tsd(jp_tem)%fdta )
                                        DEALLOCATE( sf_tsd(jp_sal)%fnow )     ! S arrays in the structure
         IF( sf_tsd(jp_sal)%ln_tint )   DEALLOCATE( sf_tsd(jp_sal)%fdta )
         IF( ln_tsd_interp )            DEALLOCATE( sf_tsd(jp_dep)%fnow )     ! T arrays in the structure
         IF( ln_tsd_interp )            DEALLOCATE( sf_tsd(jp_msk)%fnow )     ! T arrays in the structure
                                        DEALLOCATE( sf_tsd              )     ! the structure itself
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_tsd')
      !
   END SUBROUTINE dta_tsd

   !!======================================================================
END MODULE dtatsd
