      program main

      use class_opacity
      implicit none

      type(Opacity) :: op
      real(kind=8), parameter :: T = 4e4
      real(kind=8), parameter :: Pgas = 4e5
      real(kind=8) :: Rho, log10Rho, dlnRho_dlnPgas_const_T, &
            dlnRho_dlnT_const_Pgas
      real(kind=8) :: res(species), d_dlnRho_const_T(species), &
            d_dlnT_const_Rho(species), d_dabar_const_TRho(species), &
            d_dzbar_const_TRho(species)
      real(kind=8) :: kappa, dlnkap_dlnRho, dlnkap_dlnT
      integer :: ierr
      
      op = init_Opacity()
      call eos_PT(op, Pgas, T, &
            Rho, log10Rho, dlnRho_dlnPgas_const_T, &
            dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_const_TRho, d_dzbar_const_TRho, &
            ierr)
      if(ierr /= 0) stop 'Error in eos_PT'
      call kap_log10DT(op, log10(Rho), log10(T), &
            kappa, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      if(ierr /= 0) stop 'Error in lap_log10DT'
      call shutdown_Opacity(op)
 
      print *, 'temperature', T
      print *, 'pressure', Pgas
      print *, 'density', Rho
      print *, 'log10 density', log10Rho
      print *, 'kappa', kappa

      end program main
