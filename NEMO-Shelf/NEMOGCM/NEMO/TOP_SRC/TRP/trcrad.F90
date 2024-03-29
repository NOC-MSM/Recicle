MODULE trcrad
   !!======================================================================
   !!                       ***  MODULE  trcrad  ***
   !! Ocean passive tracers:  correction of negative concentrations
   !!======================================================================
   !! History :   -   !  01-01  (O. Aumont & E. Kestenare)  Original code
   !!            1.0  !  04-03  (C. Ethe)  free form F90
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_rad    : correction of negative concentrations
   !!----------------------------------------------------------------------
   USE oce_trc             ! ocean dynamics and tracers variables
   USE trc                 ! ocean passive tracers variables
   USE trd_oce
   USE trdtra
   USE prtctl_trc          ! Print control for debbuging
#if defined key_tracer_budget
   USE iom
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_rad         ! routine called by trcstp.F90

   !! * Substitutions
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE trc_rad( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad  ***
      !!
      !! ** Purpose :   "crappy" routine to correct artificial negative
      !!              concentrations due to isopycnal scheme
      !!
      !! ** Method  : - PISCES or LOBSTER: Set negative concentrations to zero
      !!                while computing the corresponding tracer content that
      !!                is added to the tracers. Then, adjust the tracer 
      !!                concentration using a multiplicative factor so that 
      !!                the total tracer concentration is preserved.
      !!              - CFC: simply set to zero the negative CFC concentration
      !!                (the total CFC content is not strictly preserved)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index      
      CHARACTER (len=22) :: charout
      ! +++>>> FABM
      INTEGER :: jn
      ! FABM <<<+++
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_rad')
      !
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_rad : Correct artificial negative concentrations '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      IF( lk_cfc     )   CALL trc_rad_sms( kt, trb, trn, jp_cfc0 , jp_cfc1               )  ! CFC model
      IF( lk_c14b    )   CALL trc_rad_sms( kt, trb, trn, jp_c14b0, jp_c14b1              )  ! bomb C14
      IF( lk_pisces  )   CALL trc_rad_sms( kt, trb, trn, jp_pcs0 , jp_pcs1, cpreserv='Y' )  ! PISCES model
      IF( lk_my_trc  )   CALL trc_rad_sms( kt, trb, trn, jp_myt0 , jp_myt1               )  ! MY_TRC model
      ! +++>>> FABM
      IF( lk_fabm  )   THEN
        DO jn=1,jp_fabm ! state variable loop
          IF (lk_rad_fabm(jn)) THEN
           CALL trc_rad_sms( kt, trb, trn, jn+jp_fabm_m1 , jn+jp_fabm_m1 )
          ENDIF
        END DO
      END IF
      ! FABM <<<+++
      !
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rad')")
         CALL prt_ctl_trc_info( charout )
         CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_rad')
      !
   END SUBROUTINE trc_rad

   SUBROUTINE trc_rad_sms( kt, ptrb, ptrn, jp_sms0, jp_sms1, cpreserv )
      !!-----------------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad_sms  ***
      !!
      !! ** Purpose :   "crappy" routine to correct artificial negative
      !!              concentrations due to isopycnal scheme
      !!
      !! ** Method  : 2 cases :
      !!                - Set negative concentrations to zero while computing
      !!                  the corresponding tracer content that is added to the
      !!                  tracers. Then, adjust the tracer concentration using
      !!                  a multiplicative factor so that the total tracer 
      !!                  concentration is preserved.
      !!                - simply set to zero the negative CFC concentration
      !!                  (the total content of concentration is not strictly preserved)
      !!--------------------------------------------------------------------------------
      !! Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER  , INTENT( in ) ::  &
         jp_sms0, &       !: First index of the passive tracer model
         jp_sms1          !: Last  index of  the passive tracer model

      REAL(wp), DIMENSION (jpi,jpj,jpk,jptra), INTENT( inout )  :: &
         ptrb, ptrn       !: before and now traceur concentration

      CHARACTER( len = 1) , INTENT(in), OPTIONAL  :: &
         cpreserv          !: flag to preserve content or not
      
      ! Local declarations
      INTEGER  :: ji, jj, jk, jn     ! dummy loop indices
      REAL(wp) :: ztrcorb, ztrmasb   ! temporary scalars
      REAL(wp) :: zcoef, ztrcorn, ztrmasn   !    "         "
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrtrdb, ztrtrdn   ! workspace arrays
#if defined key_tracer_budget
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  ztrtrdb_m1 ! slwa 
#endif
      REAL(wp) :: zs2rdt
      LOGICAL ::   lldebug = .FALSE.
      !!----------------------------------------------------------------------

 
      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrtrdb, ztrtrdn )
#if defined key_tracer_budget
      IF( kt == nittrc000 .AND. l_trdtrc) THEN
         IF (.not. ALLOCATED(ztrtrdb_m1)) ALLOCATE( ztrtrdb_m1(jpi,jpj,jpk,jptra) )  ! slwa
            DO jn = jp_sms0, jp_sms1
               IF( ln_rsttr .AND.    &                     ! Restart: read in restart  file
                  iom_varid( numrtr, 'rdb_trend_'//TRIM(ctrcnm(jn)), ldstop = .FALSE. ) > 0 ) THEN
                  IF(lwp) WRITE(numout,*) '          nittrc000-nn_dttrc RDB tracer trend read for',TRIM(ctrcnm(jn))
                  CALL iom_get( numrtr, jpdom_autoglo, 'rdb_trend_'//TRIM(ctrcnm(jn)), ztrtrdb_m1(:,:,:,jn) )   ! before tracer trend for rdb
               ELSE
                  IF(lwp) WRITE(numout,*) '          no nittrc000-nn_dttrc RDB tracer trend for',TRIM(ctrcnm(jn)),', setting to 0.'
                  ztrtrdb_m1(:,:,:,jn)=0.0
               ENDIF
            END DO
      ENDIF
#endif
      
      IF( PRESENT( cpreserv )  ) THEN   !  total tracer concentration is preserved 
      
         DO jn = jp_sms0, jp_sms1
            !                                                        ! ===========
            ztrcorb = 0.e0   ;   ztrmasb = 0.e0
            ztrcorn = 0.e0   ;   ztrmasn = 0.e0

            IF( l_trdtrc ) THEN
               ztrtrdb(:,:,:) = ptrb(:,:,:,jn)                        ! save input trb for trend computation
               ztrtrdn(:,:,:) = ptrn(:,:,:,jn)                        ! save input trn for trend computation
            ENDIF
            !                                                         ! sum over the global domain 
            ztrcorb = glob_sum( MIN( 0., ptrb(:,:,:,jn) ) * cvol(:,:,:) )
            ztrcorn = glob_sum( MIN( 0., ptrn(:,:,:,jn) ) * cvol(:,:,:) )

            ztrmasb = glob_sum( MAX( 0., ptrb(:,:,:,jn) ) * cvol(:,:,:) )
            ztrmasn = glob_sum( MAX( 0., ptrn(:,:,:,jn) ) * cvol(:,:,:) )

            IF( ztrcorb /= 0 ) THEN
               zcoef = 1. + ztrcorb / ztrmasb
               DO jk = 1, jpkm1
                  ptrb(:,:,jk,jn) = MAX( 0., ptrb(:,:,jk,jn) )
                  ptrb(:,:,jk,jn) = ptrb(:,:,jk,jn) * zcoef * tmask(:,:,jk)
               END DO
            ENDIF

            IF( ztrcorn /= 0 ) THEN
               zcoef = 1. + ztrcorn / ztrmasn
               DO jk = 1, jpkm1
                  ptrn(:,:,jk,jn) = MAX( 0., ptrn(:,:,jk,jn) )
                  ptrn(:,:,jk,jn) = ptrn(:,:,jk,jn) * zcoef * tmask(:,:,jk)
               END DO
            ENDIF
            !
            IF( l_trdtrc ) THEN
               !
               zs2rdt = 1. / ( 2. * rdt )
               ztrtrdb(:,:,:) = ( ptrb(:,:,:,jn) - ztrtrdb(:,:,:) ) * zs2rdt
               ztrtrdn(:,:,:) = ( ptrn(:,:,:,jn) - ztrtrdn(:,:,:) ) * zs2rdt 
#if defined key_tracer_budget
! slwa budget code
               DO jk = 1, jpkm1
                  ztrtrdb(:,:,jk) = ztrtrdb(:,:,jk) * e1t(:,:) * e2t(:,:) * fse3t(:,:,jk)
                  ztrtrdn(:,:,jk) = ztrtrdn(:,:,jk) * e1t(:,:) * e2t(:,:) * fse3t(:,:,jk)
               END DO
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb_m1(:,:,:,jn) )
               ztrtrdb_m1(:,:,:,jn)=ztrtrdb(:,:,:)
#else
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb )       ! Asselin-like trend handling
#endif
               CALL trd_tra( kt, 'TRC', jn, jptra_radn, ztrtrdn )       ! standard     trend handling
              !
            ENDIF

         END DO
         !
         !
      ELSE  ! total CFC content is not strictly preserved

         DO jn = jp_sms0, jp_sms1  

           IF( l_trdtrc ) THEN
              ztrtrdb(:,:,:) = ptrb(:,:,:,jn)                        ! save input trb for trend computation
              ztrtrdn(:,:,:) = ptrn(:,:,:,jn)                        ! save input trn for trend computation
           ENDIF

            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     ptrn(ji,jj,jk,jn) = MAX( 0. , ptrn(ji,jj,jk,jn) )
                     ptrb(ji,jj,jk,jn) = MAX( 0. , ptrb(ji,jj,jk,jn) )
                  END DO
               END DO
            END DO
         
            IF( l_trdtrc ) THEN
               !
               zs2rdt = 1. / ( 2. * rdt * FLOAT( nn_dttrc ) )
               ztrtrdb(:,:,:) = ( ptrb(:,:,:,jn) - ztrtrdb(:,:,:) ) * zs2rdt
               ztrtrdn(:,:,:) = ( ptrn(:,:,:,jn) - ztrtrdn(:,:,:) ) * zs2rdt 
#if defined key_tracer_budget
! slwa budget code
               DO jk = 1, jpkm1
                  ztrtrdb(:,:,jk) = ztrtrdb(:,:,jk) * e1t(:,:) * e2t(:,:) * fse3t(:,:,jk)
                  ztrtrdn(:,:,jk) = ztrtrdn(:,:,jk) * e1t(:,:) * e2t(:,:) * fse3t(:,:,jk)
               END DO
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb_m1(:,:,:,jn) )
               ztrtrdb_m1(:,:,:,jn)=ztrtrdb(:,:,:)
#else
               CALL trd_tra( kt, 'TRC', jn, jptra_radb, ztrtrdb )       ! Asselin-like trend handling
#endif
               CALL trd_tra( kt, 'TRC', jn, jptra_radn, ztrtrdn )       ! standard     trend handling
              !
            ENDIF
            !
         ENDDO

      ENDIF

#if defined key_tracer_budget
      !                                           Write in the tracer restart file
      !                                          *******************************
      IF( lrst_trc ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc : RDB trend at last time step for tracer budget written in tracer restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         DO jn = jp_sms0, jp_sms1
            CALL iom_rstput( kt, nitrst, numrtw, 'rdb_trend_'//TRIM(ctrcnm(jn)), ztrtrdb_m1(:,:,:,jn) )
         END DO
      ENDIF
#endif

      IF( l_trdtrc )  CALL wrk_dealloc( jpi, jpj, jpk, ztrtrdb, ztrtrdn )

   END SUBROUTINE trc_rad_sms
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                         NO TOP model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rad( kt )              ! Empty routine
      INTEGER, INTENT(in) ::   kt
      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rad
#endif
   
   !!======================================================================
END MODULE trcrad
