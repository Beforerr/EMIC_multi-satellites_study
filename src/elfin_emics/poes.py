import pyspedas
from pytplot import get_data
import xarray as xr
import pandas as pd
import proplot as pplt

def get_poes_mep_pro(probe = "noaa19", trange = ["2021-04-17 01:47:00", "2021-04-17 01:52:00"]):

    # NOTE: POES flux data is `support data` in POES SEM data file
    # NOTE: `xarray` does not support cdf backend
    pyspedas.poes.sem(probe = probe, trange=trange, get_support_data=True, time_clip=True, varformat='mep*')
    pyspedas.poes.sem(probe = probe, trange=trange, time_clip=True, varformat='mlt*')
    pyspedas.poes.sem(probe = probe, trange=trange, time_clip=True, varformat='l_igrf*')

    # assign coordinates mlt and l_igrf
    mep_pro_flux: xr.DataArray = get_data('mep_pro_flux',xarray=True)
    mep_pro_flux = mep_pro_flux.assign_coords({
        "mlt": ("time", get_data('mlt', xarray=True).data),
        "l_igrf": ("time", get_data('l_igrf', xarray=True).data, {"long_name": "L-shell"}),
        })

    mep_pro_flux.attrs['long_name'] = mep_pro_flux.attrs["CDF"]["VATT"]["LABLAXIS"]
    mep_pro_flux.attrs['units'] = r"$\#/cm^2/s/str/keV$"

    ds=mep_pro_flux.to_dataset(name="mep_pro_flux")
    return ds


def plot_poes_mep_pro(probe = "noaa19", trange = ["2021-04-17 01:47:00", "2021-04-17 01:52:00"], save=False, fmts=("pdf")):

    ds = get_poes_mep_pro(probe = probe, trange = trange)
    temp_ds = ds.sel(time=slice(*trange))

    start_time = pd.to_datetime(trange[0]).strftime("%H:%M")
    end_time = pd.to_datetime(trange[1]).strftime("%H:%M")
    leftlabels=('80 - 240 keV', '30-80 keV') # info copied from Themis Summary Plots

    n = 2
    fig, axs = pplt.subplots(nrows=n, refaspect=n, refwidth=2.5 * n, space=0)
    axs: list[pplt.Axes]
    axs[0].plot(x="l_igrf", y="mep_pro_flux", data=temp_ds.sel(v1_dim=0,v2_dim=1), label=r"$\parallel$")
    axs[0].plot(x="l_igrf", y="mep_pro_flux", data=temp_ds.sel(v1_dim=1,v2_dim=1), label=r"$\perp$")

    axs[1].plot(x="l_igrf", y="mep_pro_flux", data=temp_ds.sel(v1_dim=0,v2_dim=0), label=r"$\parallel$")
    axs[1].plot(x="l_igrf", y="mep_pro_flux", data=temp_ds.sel(v1_dim=1,v2_dim=0), label=r"$\perp$")

    axs[0].set_yscale("log")
    axs[1].set_yscale("log")

    axs[0].text(
        0.72, 0.92, f'Time : {start_time} - {end_time}',
        transform="axes", size='large'
    )
    axs[0].text(
        0.82, 0.82, f'<MLT> = {int(temp_ds.mlt.mean())}',
        transform="axes", size='large'
    )

    axs.format(
        leftlabels=leftlabels,
        bottomlabels=(f'{probe.upper()}',)
    )

    if save:
        fname = f"../figures/POES/{probe}_mep_pro_flux_{start_time} - {end_time}"
        for fmt in fmts:
            fig.savefig(fname+'.'+fmt)
    
    return fig