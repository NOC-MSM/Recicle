MODULE trcwri
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    TOP :   Output of passive tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_top'                                           TOP models
   !!----------------------------------------------------------------------
   !! trc_wri_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE dom_oce     ! ocean space and time domain variables
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE dianam      ! Output file name
   USE trcwri_pisces
   USE trcwri_cfc
   USE trcwri_c14b
   USE trcwri_my_trc
   ! +++>>> FABM
   USE trcwri_fabm
   ! FABM <<<+++

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri      

   !! * Substitutions
#  include "top_substitute.h90"

CONTAINS

#if defined key_tracer_budget
   SUBROUTINE trc_wri( kt , fl)  !slwa
#else
   SUBROUTINE trc_wri( kt )
#endif
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri  ***
      !! 
      !! ** Purpose :   output passive tracers fields and dynamical trends
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )     :: kt
      ! +++>>>FABM
#if defined key_tracer_budget
      INTEGER, INTENT( in ), OPTIONAL     :: fl  ! slwa
#endif
      ! FABM <<<+++
      !
      INTEGER                   :: jn
      CHARACTER (len=20)        :: cltra
      CHARACTER (len=40)        :: clhstnam
      INTEGER ::   inum = 11            ! temporary logical unit
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_wri')
      !
      IF( lk_offline .AND. kt == nittrc000 .AND. lwp ) THEN    ! WRITE root name in date.file for use by postpro
         CALL dia_nam( clhstnam, nn_writetrc,' ' )
         CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         WRITE(inum,*) clhstnam
         CLOSE(inum)
      ENDIF
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( lk_pisces  )   CALL trc_wri_pisces     ! PISCES 
      IF( lk_cfc     )   CALL trc_wri_cfc        ! surface fluxes of CFC
      IF( lk_c14b    )   CALL trc_wri_c14b       ! surface fluxes of C14
      ! +++>>>FABM
#if defined key_tracer_budget
      IF( PRESENT(fl) ) THEN
         IF( lk_fabm    )   CALL trc_wri_fabm (kt, fl) ! MY_TRC  tracers for budget
         IF( lk_my_trc ) CALL trc_wri_my_trc (kt, fl)    ! MY_TRC  tracers for budget
      ELSE
         IF( lk_fabm    )   CALL trc_wri_fabm (kt) ! FABM  tracers for budget
         IF( lk_my_trc  )   CALL trc_wri_my_trc (kt) ! MY_TRC  tracers
      ENDIF
#else
      IF( lk_fabm  )   CALL trc_wri_fabm (kt)     ! FABM  tracers
      IF( lk_my_trc  )   CALL trc_wri_my_trc(kt)     ! MY_TRC  tracers
#endif
      ! FABM <<<+++
      !

      IF( nn_timing == 1 )  CALL timing_stop('trc_wri')
      !
   END SUBROUTINE trc_wri

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri
CONTAINS
   SUBROUTINE trc_wri( kt )                     ! Empty routine   
   INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_wri
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri
