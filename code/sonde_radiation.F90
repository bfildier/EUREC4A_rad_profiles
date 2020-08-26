subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rte_rrtmgp_clouds stopping"
    stop
  end if
end subroutine stop_on_err
! -------------------------------------------------------------
program sonde_radiation
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use netcdf
  use mo_simple_netcdf,      only: get_dim_size, read_field, write_field
  implicit none
  ! -------------------------------------------------------------
  integer, parameter                   :: n_gases = 8
  character(len=3), dimension(n_gases) :: gases = ["h2o", "co2", "ch4", "n2o", "o3 ", "o2 ", "n2 ", "co "]
  integer :: igas, nlay, nargs
  integer :: ncid
  logical :: top_at_1
  character(len=256) :: file_name, rrtmgp_dir, lw_coeffs_file, sw_coeffs_file
  real(wp), dimension(:,:), allocatable :: play, plev, tlay
  real(wp), dimension(:),   allocatable :: gas_conc
  real(wp), dimension(:),   allocatable :: sfc_t, mu0
  real(wp), dimension(:,:), allocatable, &
                                 target :: flux_up, flux_dn, flux_net
  real(wp), dimension(:,:), allocatable :: sfc_emis, sfc_alb, sw_toa

  type(ty_gas_concs)          :: gas_concs
  type(ty_gas_optics_rrtmgp)  :: gas_optics_lw, gas_optics_sw
  type(ty_source_func_lw)     :: lw_sources
  type(ty_optical_props_1scl) :: lw_optical_props
  type(ty_optical_props_2str) :: sw_optical_props
  type(ty_fluxes_broadband)   :: fluxes

  ! -------------------------------------------------------------
  !
  ! Default absorption coefficent files
  !
  
  rrtmgp_dir     = "../rte-rrtmgp/rrtmgp/data"
  lw_coeffs_file = "rrtmgp-data-lw-g256-2018-12-04.nc"
  sw_coeffs_file = "rrtmgp-data-sw-g224-2018-12-04.nc"
  print *, "Usage: sonde_radiation combo_file [rrtmgp_dir] [lw_coeffs_file] [sw_coeffs_file]"
  nargs = command_argument_count()
  if(nargs < 1) call stop_on_err("sonde_radiation: have to specify at least which file you want radiation computed for")
  call get_command_argument(1, file_name)
  if(nargs >= 2) call get_command_argument(2, rrtmgp_dir)
  if(nargs >= 3) call get_command_argument(3, lw_coeffs_file)
  if(nargs >= 4) call get_command_argument(4, sw_coeffs_file)
  if(nargs >  4) print *, "ignoring arguments beyond #4"

  !
  ! Read data from the profile file
  !
  if(nf90_open(trim(file_name), NF90_WRITE, ncid) /= NF90_NOERR) &
    call stop_on_err("sonde_radiation: can't find file " // trim(file_name))
  nlay = get_dim_size(ncid, "play")
  if(get_dim_size(ncid, "plev") /= nlay+1) &
    call stop_on_err("sonde_radiation: layer and level dimensions in file inconsistent")

  allocate(play(1,nlay),  tlay(1,nlay), plev(1,nlay+1))
  play(1,:) = read_field(ncid, "play", nlay)
  tlay(1,:) = read_field(ncid, "tlay", nlay)
  plev(1,:) = read_field(ncid, "plev", nlay+1)
  top_at_1 = play(1, 1) < play(1,nlay)

  !
  ! Gas concentrations (volume mixing ratios) go in a type
  !
  allocate(gas_conc(nlay))
  call stop_on_err(gas_concs%init(gases))
  do igas = 1, n_gases
    gas_conc(:) = read_field(ncid, trim(gases(igas)), nlay)
    call stop_on_err(gas_concs%set_vmr(trim(gases(igas)), gas_conc))
  end do
  allocate(sfc_t(1))
  sfc_t(:) = read_field(ncid, "sfc_t")

  allocate(flux_up(1, nlay+1), flux_dn(1, nlay+1), flux_net(1, nlay+1))
  fluxes%flux_dn  => flux_dn
  fluxes%flux_up  => flux_up
  fluxes%flux_net => flux_net

  ! -------------------------------------------------------------
  !
  ! Radiation, finally
  !

  !
  ! Longwave
  !
  call load_and_init(gas_optics_lw, trim(rrtmgp_dir) // "/" // trim(lw_coeffs_file), gas_concs)
  call stop_on_err(lw_sources%alloc           (1, nlay, gas_optics_lw))
  call stop_on_err(lw_optical_props%alloc_1scl(1, nlay, gas_optics_lw))
  allocate(sfc_emis(gas_optics_lw%get_nband(), 1))
  sfc_emis(:,1) = read_field(ncid, "sfc_emis")
  !
  ! Gas optics
  !
  call stop_on_err(gas_optics_lw%gas_optics(play, plev, tlay, &
                                            sfc_t, gas_concs, &
                                            lw_optical_props, lw_sources))
  !
  ! Radiative transfer
  !
  call stop_on_err(rte_lw(lw_optical_props,  top_at_1, &
                          lw_sources, sfc_emis,        &
                          fluxes, n_gauss_angles = 3))
  call stop_on_err(write_field(ncid, "lw_dn",  flux_dn (1,:)))
  call stop_on_err(write_field(ncid, "lw_up",  flux_up (1,:)))
  call stop_on_err(write_field(ncid, "lw_net", flux_net(1,:)))

  !
  ! Shortwave
  !
  call load_and_init(gas_optics_sw, trim(rrtmgp_dir) // "/" // trim(sw_coeffs_file), gas_concs)
  allocate(mu0(1), &
           sfc_alb(gas_optics_sw%get_nband(), 1), &
           sw_toa(1, gas_optics_sw%get_ngpt()))
  sfc_alb(:,:) = read_field(ncid, "sfc_alb")
  mu0(:)       = read_field(ncid, "cos_sza")
  call stop_on_err(sw_optical_props%alloc_2str(1, nlay, gas_optics_sw))
  !
  ! Gas optics
  !
  call stop_on_err(gas_optics_sw%gas_optics(play, plev, tlay, &
                                            gas_concs,        &
                                            sw_optical_props, sw_toa))
  !
  ! Radiative transfer
  !
  call stop_on_err(rte_sw(sw_optical_props, top_at_1,      &
                          mu0, sw_toa, sfc_alb, sfc_alb, &
                          fluxes))
  call stop_on_err(write_field(ncid, "sw_dn",  flux_dn (1,:)))
  call stop_on_err(write_field(ncid, "sw_up",  flux_up (1,:)))
  call stop_on_err(write_field(ncid, "sw_net", flux_net(1,:)))

  ncid = nf90_close(ncid)

end program sonde_radiation
