

platform = 'BCO'


import xarray as xr

data = xr.open_dataset('../output/rad_profiles_proxy_sondes_fixedQ.nc')

data_new = data.where(data.platform == platform,drop=True)

data_new.to_netcdf('../output/rad_profiles_proxy_sondes_fixedQ_wo_%s.nc'%platform)
