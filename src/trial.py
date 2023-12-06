
# %%
import pdpipe as pdp
import altair as alt

l_fitler = pdp.PdPipeline([
    pdp.drop_rows_where['l_igrf'] > 8,
    pdp.drop_rows_where['l_igrf'] < 2,
])

mlt_dawn_fitler = pdp.PdPipeline([
    (pdp.keep_rows_where['mlt'] < 6) | (pdp.keep_rows_where['mlt'] > 21) ,
])

mlt_dusk_fitler = pdp.PdPipeline([
    pdp.keep_rows_where['mlt'] > 12 ,
    pdp.keep_rows_where['mlt'] < 22 ,
])

pro_ratio_filter = pdp.PdPipeline([
    pdp.drop_rows_where['mep_pro_flux_ratio'] > 1 ,
])

ele_ratio_filter = pdp.PdPipeline([
    pdp.drop_rows_where['mep_ele_flux_ratio'] > 1 ,
])


def plot_mep_ele_flux_ratio(df):
    return alt.Chart(df).mark_rect().encode(
        alt.X("hours(time):O"),
        alt.Y("l_igrf:Q", bin=True),
        alt.Color("mean(mep_ele_flux_ratio):Q",
                scale=alt.Scale(type="log"),title="")
    ).properties(
        title="Electron flux ratio (parallel/perpendicular) at Dawn sector (MLT < 6 or MLT > 21)"
    )

def plot_mep_pro_flux_ratio(df):
    return alt.Chart(df).mark_rect().encode(
        alt.X("hours(time):O"),
        alt.Y("l_igrf:Q", bin=True),
        alt.Color("mean(mep_pro_flux_ratio):Q",
                scale=alt.Scale(type="log"),title="")
    ).properties(
        title="Proton flux ratio (parallel/perpendicular) at Dusk sector (MLT > 12 and MLT < 22)"
    )

# file_list = glob.glob("./data/*.sav")
# dfs = list(map(read_sav, file_list))

# df = dfs[0]
# Select time between "01:15" and "01:20", time is index
# df = df[(df.index >= pd.to_datetime("2021-04-17 01:15:00")) & (df.index <= pd.to_datetime("2021-04-17 01:20:00"))]
# df["mep_pro_flux_tel0_0"].hvplot(logy=True) * df["mep_pro_flux_tel90_0"].hvplot(logy=True)

# for df in dfs[0:1]:
#     p_chart_1 = plot_mep_pro_flux_ratio(l_fitler(df))
#     p_chart_2 = plot_mep_pro_flux_ratio(l_fitler(mlt_dusk_fitler(df)))
#     display(p_chart_1 | p_chart_2)

# for df in dfs[0:1]:
#     e_chart_0 = plot_mep_ele_flux_ratio(l_fitler(df))
#     e_chart_1 = plot_mep_ele_flux_ratio(l_fitler(mlt_dawn_fitler(df)))
#     display(e_chart_0 | e_chart_1)

# p_chart_0 = plot_mep_pro_flux_ratio(l_fitler(pd.concat(dfs))) 
# p_chart_1 = plot_mep_pro_flux_ratio(mlt_dusk_fitler(l_fitler(pd.concat(dfs))))
# p_chart_2 = plot_mep_pro_flux_ratio(pro_ratio_filter(mlt_dusk_fitler(l_fitler(pd.concat(dfs)))))
# e_chart_0 = plot_mep_ele_flux_ratio(l_fitler(pd.concat(dfs)))
# e_chart_1 = plot_mep_ele_flux_ratio(mlt_dawn_fitler(l_fitler(pd.concat(dfs))))
# e_chart_2 = plot_mep_ele_flux_ratio(ele_ratio_filter(mlt_dawn_fitler(l_fitler(pd.concat(dfs)))))

# display(
#     p_chart_0 | p_chart_1 , 
#     e_chart_0 | e_chart_1 , 
#     p_chart_2 ,  e_chart_2
# )



#| code-summary: Load the data using `xarray`
def load_goesr(trange, probe=16) -> xr.Dataset:

    # download the data files using `pyspedas`
    # can not load file directly from `pyspedas` because `OSError: [Errno -101] NetCDF: HDF error: '/data/goes/goes16/l2/data/mpsh-l2-avg1m/2021/04/sci_mpsh-l2-avg1m_g16_d20210417_v1-0-3.nc'`
    
    files = pyspedas.goes.mpsh(trange=trange, probe = probe, downloadonly=True)
    ds = xr.open_mfdataset(files,combine='nested',concat_dim='record_number')

    coords_dict = {
        "time": ds["L2_SciData_TimeStamp"],
    }

    return ds.assign_coords(coords_dict)


#| code-summary: Plot average differential electron flux for each channel by telescopes
def plot_goes_group(group):
    name, data = group
    return (data.sum(dim='telescopes').hvplot(label='sum',logy=True) * data.hvplot(by='telescopes',logy=True))

def plot_goes_eflux_all_channels(ds):
    return hv.Layout(
        list(map(plot_goes_group, ds["AvgDiffElectronFlux"].groupby("electron_diff_chans")))
    ).opts(shared_axes=False)

# plot_goes_eflux_all_channels(ds)


#| code-summary: Plot average differential electron flux for each channel summed over telescopes
def plot_goes_eflux(ds,backend='hvplot'):

    da = ds["AvgDiffElectronFlux"].sum(dim='telescopes')

    match backend:
        case "hvplot":
            return da.hvplot(by='electron_diff_chans',logy=True)
        
        case "proplot":
            # proplot can not deal with Dask arrays, so we have to convert to numpy
            da = xr.DataArray(
                da.to_numpy(), coords=[da.time, da.electron_diff_chans], name='Time-averaged electron fluxes', attrs={'units': da.attrs['units']}
            )

            fig = pplt.figure(refaspect=2,refwidth=5)
            ax = fig.subplot()

            ax.plot(da)

            ax.legend(loc='r', ncols=1, label='Energy channels', frame=False)
            ax.format(
                xlabel='Time',xformatter="concise",
                yscale="log",yformatter="log",ylim=(1, 1e7),
            )

            return fig
            # fig.savefig('../figures/goes_eflux.svg')
        case _:
            raise ValueError(f"Unknown backend: {backend}")

# plot_goes_eflux(ds, backend='proplot')


#%%

#| code-summary: Create a new dataframe with column name 'B_int' which are the magnetic field intensity and column name 'freq' which is the average frequency

def process(raw_xr):
    raw_df = raw_xr.to_dataframe()
    _df = (
        raw_df.spec_bins[1] - raw_df.spec_bins[0]
    ) * 1000  # difference between each frequency bin in Hz
    # get the magnetic field intensity
    df = (
        (_df * raw_df["erg_pwe_ofa_l2_spec_B_spectra_132"]).groupby("time").sum()
    ).to_frame(name="wave_intensity")

    # get the average frequency
    df["wave_freq"] = (
        (1000 * raw_df["spec_bins"])
        * (_df * raw_df["erg_pwe_ofa_l2_spec_B_spectra_132"])
    ).groupby("time").sum() / df["wave_intensity"]
    return df

# pipeline = pdp.PdPipeline(
#     [
#         pdp.df["Lm"] << orb_Lm_itrp.data[:, 0],
#         pdp.df["mlat"] << orb_mlat_itrp,
#         pdp.df["B_eq"] << orb_beq_itrp,
#         pdp.DropNa(),
#     ]
# )

#| code-summary: Filter the summed data to only include values where the summed intensity is greater than 10 pT^2/Hz, |mlat|<15
# pipeline_filter = pdp.PdPipeline(
#     [
#         # Filter the summed data to only include values where the summed intensity is greater than 10 pT^2/Hz
#         pdp.drop_rows_where["wave_intensity"] < 10,
#         # Filter the data to only include data with |mlat|<15
#         pdp.RowDrop({"mlat": lambda x: np.abs(x) > 15}),
#     ]
# )

# pwe_B_int = pipeline(process(pwe_spec_filtered))
# pwe_B_filtered = pipeline_filter(pwe_B_int)
# pwe_B_filtered.hvplot.scatter(x="wave_freq",y="wave_intensity")