import pandas as pd
import proplot as pplt
import altair as alt
import pdpipe as pdp
import numpy as np
import xarray as xr
import glob
import pyspedas
from pytplot import get_data
import unyt as un

DEFAULT_FMTS = ('svg', 'pdf')

def label_formatter(da: xr.DataArray):
    if "units" in da.attrs.keys() and da.units != "":
        return f"{da.long_name} ({da.units})"
    else:
        return da.long_name


def plot_erg_orb(df: pd.DataFrame, backend="proplot"):

    if backend == "altair":
        chart1 = (
            alt.Chart(df.reset_index())
            .mark_line()
            .encode(
                x="time",
                y="Lm",
            )
        )

        chart2 = (
            alt.Chart(df.reset_index())
            .mark_line()
            .encode(
                x="time",
                y="mlat",
            )
        )

        chart3 = (
            alt.Chart(df)
            .mark_point()
            .encode(
                x="Lm",
                y="mlat",
            )
        )

        return chart1 & chart2 & chart3

    elif backend == "proplot":
        # df = df.reset_index()
        fig, axs = pplt.subplots(ncols=3, nrows=1, sharey=False)
        
        axs[0].plot('Lm', data=df)
        axs[1].plot('mlat', data=df)
        axs[2].scatter(df['Lm'],df['mlat'])
        # return fig

def plot_erg_wave(df: pd.DataFrame):
    return df.hvplot.scatter(
        x="Lm", y="wave_intensity", title="Intensity [pT^2]", logy=True
    ) + df.hvplot.scatter(
        x="Lm", y="wave_freq", title="Frequency [kHz]", datashade=True
    )


# Plot the intensity (log scale) vs Lm using altair
def plot_erg_wave_intensity(df: pd.DataFrame):
    chart1 = (
        alt.Chart(df)
        .mark_point()
        .encode(
            x=alt.X("Lm", scale=alt.Scale(zero=False)),
            y=alt.Y(
                "wave_intensity",
                scale=alt.Scale(type="log"),
                axis=alt.Axis(title="Intensity [pT^2]"),
            ),
        )
    )

    chart2 = (
        alt.Chart(df)
        .mark_line()
        .encode(
            x=alt.X("Lm", bin=alt.Bin(maxbins=40)),
            y=alt.Y(
                "mean(wave_intensity)",
                scale=alt.Scale(type="log"),
                axis=alt.Axis(title="Mean Intensity [pT^2]"),
            ),
        )
    )

    band2 = (
        alt.Chart(df)
        .mark_errorband(extent="ci")
        .encode(
            x=alt.X("Lm", bin=alt.Bin(maxbins=40)),
            y=alt.Y(
                "wave_intensity",
                scale=alt.Scale(type="log"),
                axis=alt.Axis(title="Mean Intensity [pT^2]"),
            ),
        )
    )

    (chart1 | chart2 + band2).display()


def plot_erg_wave_Lm_mlat(df: pd.DataFrame, backend="altair"):

    """
    Parameters
    ----------
    df (pd.DataFrame): dataframe containing Lm and mlat data
    Returns: altair chart

    Notes
    -----
    backend can be either "altair" or "hvplot", but "hvplot" is not yet implemented completely
    """
    if backend == "altair":
        heatmap = (
            alt.Chart(df.dropna())
            .mark_rect()
            .encode(
                alt.X("Lm", bin=alt.Bin(maxbins=50)),
                alt.Y("mlat", bin=alt.Bin(maxbins=50)),
                alt.Color(
                    "mean(wave_intensity)",
                    scale=alt.Scale(type="log"),
                    legend=alt.Legend(title="Mean Intensity [pT^2]"),
                ),
            )
        )

        points = (
            alt.Chart(df.dropna())
            .mark_circle(
                color="black",
                size=1,
            )
            .encode(
                x="Lm",
                y="mlat",
            )
        )
        return heatmap + points

    elif backend == "hvplot":
        return df.hvplot.heatmap(
            x="Lm",
            y="mlat",
            C="wave_intensity",
            logz=True,
            cmap="viridis",
            xlabel="Lm",
            ylabel="mlat",
            clabel="wave_intensity",
        )


# PLOT: Plot the wave statistics
# (plot_erg_wave(pwe_B_int) + plot_erg_wave(pwe_B_filtered)).cols(2) if PLOT else None

# PLOT: Plot wave intensity vs Lm
# plot_erg_wave_intensity(pwe_B_filtered) if PLOT else None

# PLOT: Plot a log heatmap from binned Lm and mlat data using altair
# plot_erg_wave_Lm_mlat(pwe_B_int).display() if PLOT else None
# plot_erg_wave_Lm_mlat(pwe_B_filtered).display() if PLOT else None



# TODO
# if PLOT:
#     _point = (
#         alt.Chart(pwe_B_energy)
#         .mark_circle(opacity=0.2)
#         .encode(
#             x=alt.X("Lm", scale=alt.Scale(zero=False)),
#             y=alt.Y(
#                 "D_whistler", scale=alt.Scale(type="log"), axis=alt.Axis(format="e")
#             ),
#         )
#     )

#     _line = (
#         alt.Chart(pwe_B_energy)
#         .mark_line(color="red")
#         .encode(
#             x=alt.X("Lm", bin=alt.Bin(maxbins=40), title="L_m"),
#             y=alt.Y("mean(D_whistler)", title="D_ch [rad/s]"),
#         )
#     )

#     (_point + _line).display()

# if PLOT:
#     _point1 = (
#         alt.Chart(pwe_B_energy)
#         .mark_circle(opacity=0.2)
#         .encode(
#             x=alt.X("Lm", scale=alt.Scale(zero=False)),
#             y=alt.Y("B_ratio", scale=alt.Scale(type="log"), axis=alt.Axis(format="e")),
#         )
#     )

#     _point2 = (
#         alt.Chart(pwe_B_energy)
#         .mark_circle(opacity=0.2)
#         .encode(
#             x=alt.X("Lm", scale=alt.Scale(zero=False)),
#             y=alt.Y("Ω_ratio"),
#         )
#     )
#     _point1 | _point2


# %%

from functools import partial

pipe_time = pdp.PdPipeline(
    [pdp.AggByCols("time", partial(pd.to_datetime, unit="s")), pdp.SetIndex("time")]
)

pipeline_filter_j = pdp.PdPipeline(
    [
        pdp.df.reset_index("energy"),
        # "J_ratio" should be in the range of [0.11, 1], fitler out the rest
        # NOTE: should we filter out the larger than one part?
        (pdp.keep_rows_where["J_ratio"] >= 0.11)
        # & (pdp.keep_rows_where["J_ratio"] <= 1),
    ]
)

def unit_df(raw_df: pd.DataFrame) -> pd.DataFrame:
    r"""
    Converts the DataFrame to unyt units
    >>> unit_df(df1.head())
    """
    type_dic = {
        "alpha_LC": un.deg,
        "elx_Btot": un.nT,
        "elx_foot_Btot": un.nT,
        "energy": un.keV,
        "velocity": un.m / un.s,
        "tau_B": un.s,
        "D_ch": un.s**-1,
    }

    df = raw_df.copy()
    for col, unit in type_dic.items():
        if col in df.columns:
            df[col] = df[col] * unit
    return df


# NOTE: `pdpipe` does not work
# unit_df = pdp.PdPipeline([
#     pdp.ApplyByCols("alpha_LC", lambda x: x * u.deg),
#     pdp.df["alpha_LC (unit)"] << pdp.df["alpha_LC"] * u.deg,
#     pdp.df["elx_Btot"] << pdp.df["elx_Btot"] * u.nT,
#     pdp.df["elx_foot_Btot"] << pdp.df["elx_foot_Btot"] * u.nT,
#     pdp.df["L"] << pdp.df["L"] * u.dimensionless,
# ])

# NOTE: `pint` does not work quite right with pandas
# import pint
# import pint_pandas
# def unit_df(raw_df: pd.DataFrame) -> pd.DataFrame:
#     type_dic = {
#         "alpha_LC": "pint[degree]",
#         "elx_Btot": "pint[nT]",
#         "elx_foot_Btot": "pint[nT]",
#         "L": "pint[dimensionless]",
#         "energy": "pint[keV]"
#     }
#     return raw_df.astype(type_dic)


# NOTE: `astropy.unit` does not work quite right with pandas
# from astropy.constants import c, m_e
# from astropy import units as u

# df1["alpha_LC"] = df1["alpha_LC"] * u.deg
# df1["elx_Btot"] = df1["elx_Btot"] * u.nT
# df1["elx_foot_Btot"] = df1["elx_foot_Btot"] * u.nT
# df1["L"] = df1["L"] * u.dimensionless_unscaled



#| code-summary: Load CSV files into DataFrame
def get_elf_df() -> pd.DataFrame:
    # Get a list of all CSV files in a directory
    file_list = glob.glob("../data/elx_orb_*.csv")

    # Create an empty list to store DataFrames
    dfs = []

    # Iterate over each file in the list and read it into a DataFrame
    for file in file_list:
        df = pd.read_csv(file)
        dfs.append(df)

    # Concatenate all DataFrames into a single DataFrame
    elf_orb_df = pd.concat(dfs, ignore_index=True)
    elf_orb_df = pipe_time(elf_orb_df)

    # Get a list of all CSV files in a directory
    file_list = glob.glob("../data/elx_j_ratio_*.csv")

    dfs = []
    for file in file_list:
        arr = np.genfromtxt(file, delimiter=",", skip_header=1)
        ebins = np.genfromtxt(file, delimiter=",", max_rows=1)[1:]
        time = pd.to_datetime(arr[:, 0], unit="s")
        # create a DataArray with coordinates 'time' and 'energy'
        da = xr.DataArray(
            arr[:, 1:],
            dims=["time", "energy"],
            coords={"time": time, "energy": ebins},
            name="J_ratio",
        )

        dfs.append(da.to_dataframe())

    # Concatenate all DataFrames into a single DataFrame
    elf_j_df = pd.concat(dfs)
    elf_j_df = pipeline_filter_j(elf_j_df)
    return pd.merge(elf_orb_df, elf_j_df, left_index=True, right_index=True)


def Ty(alpha):
    y = np.sin(alpha)
    C1 = 1.3809
    C2 = -0.1851
    C3 = -0.4559
    n = 0.863
    return C1 + C2 * y ** (0.5) + C3 * y**n


def tau_B(L, v, Ty):
    Re = 6371e3 * un.m  # unit: m
    return 4.0 * L * Re / v * Ty


def velocity(erg):
    return c * np.sqrt(1.0 - 1.0 / (erg / (un.me * un.c**2.0) + 1) ** 2.0)


def cal_df(raw_df: pd.DataFrame) -> pd.DataFrame:
    # deep copy
    df = unit_df(raw_df)
    df["velocity"] = velocity(df["energy"].values)
    df["alpha_eq_LC"] = np.arcsin(
        (df["elx_foot_Btot"].values / df["elx_Btot"].values) ** 0.5
        * np.sin(df["alpha_LC"].values)
    )
    df["Ty"] = Ty(df["alpha_eq_LC"].values)
    df["tau_B"] = tau_B(df["L"].values, df["velocity"].values, df["Ty"].values)
    df["z_0"] = 0.9 / df["J_ratio"].values
    df["D_EMIC"] = (
        4.0
        * (df["alpha_eq_LC"].values) ** 2.0
        / df["tau_B"].values
        / df["z_0"].values ** 2.0
    )

    return df


def plot_elf_energy_D_EMIC(df, **kwargs):
    points = (
        alt.Chart(df)
        .mark_point(
            opacity=0.2,
        )
        .encode(
            x=alt.X("energy", scale=alt.Scale(type="log"), title="Energy [keV]"),
            y=alt.Y("D_EMIC", scale=alt.Scale(type="log"), title="D_EMIC [s^-1]"),
            color=alt.Color(
                "L", bin=alt.Bin(nice=False, **kwargs), scale=alt.Scale(scheme="dark2")
            ),
        )
    )
    line = (
        alt.Chart(df)
        .mark_line()
        .encode(
            x=alt.X("energy", scale=alt.Scale(type="log")),
            y=alt.Y("mean(D_EMIC)", scale=alt.Scale(type="log")),
            color=alt.Color(
                "L", bin=alt.Bin(nice=False, **kwargs), scale=alt.Scale(scheme="dark2")
            ),
        )
    )
    return points + line

# %%
#| code-summary: Plot diffusion coefficients of EMIC wave over energy in different L shells

# elf_raw_df = get_elf_df()
# ebins = elf_raw_df["energy"].unique()
# elf_raw_df.describe()
# elf_df = cal_df(elf_raw_df)
# display(elf_df.head(),elf_df.describe())
# plot_elf_energy_D_EMIC(elf_df.where((elf_df.L > 3.5) & (elf_df.L < 6.5)), step=1)



# %%
from astropy import units as u
#| code-summary: Calculate the whistler mode diffusion coefficient.
def D_whistler(B_int: u.T**2, B_eq: u.T, L_m, ω_m, enerygy: u.eV):
    """Calculate the whistler mode diffusion coefficient."""
    from astropy.constants import c, m_e

    Ω_ce0 = gyrofrequency(B_eq, "e-")
    # trough density model
    Ω_ce0_over_Ω_pe0 = 1 / L_m
    # relativistic factor
    gamma = 1 + enerygy / (m_e * c**2)
    # whistler mode waves have a typical wave-normal angle width Δ𝜃CH ≈ 30°–45°
    dtheta_ch = 30 * u.deg
    B_ratio = (B_int / B_eq**2).to(u.dimensionless_unscaled)
    Ω_ratio = (Ω_ce0 / (ω_m)).to(u.dimensionless_unscaled)
    result = (
        B_ratio
        * Ω_ce0
        * Ω_ce0_over_Ω_pe0
        * (Ω_ratio) ** 0.5
        / np.tan(dtheta_ch)
        / gamma**2
        / 2
        / np.sqrt(3)
    ).to(u.rad / u.s)
    return result


def pd_df_D_whistler(df: pd.DataFrame):
    B_int = df["wave_intensity"].to_numpy() * u.pT**2
    B_eq = df["B_eq"].to_numpy() * u.nT
    L_m = df["Lm"].to_numpy()
    ω_m = (df["wave_freq"].to_numpy() * u.kHz).to(
        u.rad / u.s, equivalencies=[(u.cy / u.s, u.Hz)]
    )
    energy = df["energy"].to_numpy() * u.keV
    return D_whistler(B_int, B_eq, L_m, ω_m, energy)


pipeline_whistler = pdp.PdPipeline(
    [
        pdp.ApplyToRows(lambda row: ebins, "ebins"),
        pdp.df.explode("ebins"),
        pdp.ColRename({"ebins": "energy"}),
        pdp.ColByFrameFunc("D_whistler", pd_df_D_whistler),
    ]
)


# erg_pwe_D = pipeline_whistler(pwe_B_filtered)
# erg_pwe_D.describe()


def plot_diffusion_coefficient_energy_Lm_bins(df):
    """
    Plots the diffusion coefficient as a function of energy in bins of Lm.

    Parameters:
    df (pandas.DataFrame): DataFrame containing columns for bin values, energy values, and diffusion coefficients.

    Returns:
    None.
    """

    chart = (
        alt.Chart(df)
        .mark_line(point=alt.OverlayMarkDef(filled=False, fill="white"))
        .encode(
            x=alt.X(
                "energy",
                scale=alt.Scale(type="log"),
                title="Energy [keV]",
            ),
            y=alt.Y("mean(D_whistler)", scale=alt.Scale(type="log")),
            color=alt.Color(
                "Lm", bin=alt.Bin(nice=False,step=1), scale=alt.Scale(scheme="dark2")
            ),
        )
    )

    # Save the chart as an HTML file and display it
    # chart.save('diffusion_coefficient_energy_Lm_bins.html')
    return chart

# plot_diffusion_coefficient_energy_Lm_bins(erg_pwe_D) if PLOT else None


# %%
### Calculate combined diffusion coefficient

#| code-summary: Calculate the lifetime of electrons

from numpy import log, sin
from unyt import deg, rad

# alpha_0_max_EMIC = 50 * deg
# alpha_0_max_CH = 87.0 * deg

# tau_CH = lambda D_whistler: log(sin(alpha_0_max_CH) / sin(alpha_0_max_EMIC)) / (
#     4 * D_whistler
# )
# tau_EMIC = lambda D_aa, alpha_eq_LC: log(sin(alpha_0_max_EMIC) / sin(alpha_eq_LC)) / (
#     4 * D_aa
# )


# res_df = pd.concat(
#     [
#         erg_pwe_D.groupby("energy")["D_whistler"].mean(),
#         elf_df.groupby("energy")["D_EMIC", "alpha_eq_LC"].mean(),
#     ],
#     axis=1,
# )

# res_df["tau_CH"] = tau_CH(res_df["D_whistler"])
# res_df["tau_EMIC"] = tau_EMIC(res_df["D_EMIC"], res_df["alpha_eq_LC"])
# res_df["tau"] = res_df["tau_CH"] + res_df["tau_EMIC"]


#| code-summary: Format the time in a pretty way

def pretty_format_time(time:pd.Series):
    """Convert the times in a pretty format according to its mean value"""
    time = time.to_numpy() * un.second
    if time.mean() > 1 * un.day:
        return time.to(un.day)
    elif time.mean() > 1 * un.hour:
        return time.to(un.hour)
    elif time.mean() > 1 * un.minute:
        return time.to(un.minute)
    else:
        return time.to(un.second)

def pretty_format_time_df(raw_df:pd.DataFrame):
    """Convert the times in a pretty format according to its mean value"""
    df = raw_df.copy()
    for column in pdp.cq.StartsWith("tau")(df):
        df[column] = pretty_format_time(df[column])
    return df

# pretty_format_time_df(res_df)


def plot_lifetime(df, backend="hvplot"):

    xlabel = "Energy [keV]"
    ylabel = "Electron Lifetime [s]"
    df = df[["tau","tau_CH","tau_EMIC"]]
    match backend:
        case "hvplot":
            return df.hvplot(
                logx = True, xlabel = xlabel, xlim = (3e1, 4e3),
                logy = True, ylabel = ylabel,
            ) * df.hvplot.scatter()
        
        case "proplot":
            fig = pplt.figure()
            ax = fig.subplot()

            ax.plot(df, cycle='default')
            ax.legend( ncols=1, frame=False)

            ax.scatter(df, cycle='default')

            ax.format(
                xscale='log', xlabel=xlabel, xlim=(5e1, 3e3),
                yscale='log', ylabel=ylabel,
            )
            fig.savefig("../figures/lifetime.png")

# res_df.hvplot(y=["tau","tau_CH","tau_EMIC"],logy=True,title="Electron Lifetime")
# plot_lifetime(res_df, backend="proplot") if PLOT else None
# plot_lifetime(res_df, backend="hvplot") if PLOT else None

# %%
#| code-summary: read pose data from .sav file (processed first using IDL)
from scipy.io import readsav
import glob
import xarray as xr


def read_sav(file, xrrary=False):
    data = readsav(file)
    probe = data.probe.decode()

    def get_data(key):
        key = f"{probe}_{key}"
        dtype = data[key].y[0].dtype.name
        return data[key].y[0].astype(dtype)

    print(probe, data.keys())

    if xrrary:
        ds = xr.Dataset(
            data_vars={
                "mep_ele_flux_tel0": (
                    ("electron_energy", "time"),
                    get_data("mep_ele_flux_tel0"),
                ),
                "mep_ele_flux_tel90": (
                    ("electron_energy", "time"),
                    get_data("mep_ele_flux_tel90"),
                ),
                "mep_pro_flux_tel0": (
                    ("proton_energy", "time"),
                    get_data("mep_pro_flux_tel0"),
                ),
                "mep_pro_flux_tel90": (
                    ("proton_energy", "time"),
                    get_data("mep_pro_flux_tel90"),
                ),
            },
            coords={
                "time": pd.to_datetime(data[f"{probe}_l_igrf"].x[0], unit="s"),
                "l_igrf": ("time", get_data("l_igrf"), {"long_name": f"L-Shell"}),
                "mlt": ("time", get_data("mlt")),
            },
            attrs={
                "probe": probe
            }
        )
        # make mlt as a new coordinates
        return ds
    else:
        df = pd.DataFrame(
            {
                "time": pd.to_datetime(data[f"{probe}_l_igrf"][0][0], unit="s"),
                "l_igrf": data[f"{probe}_l_igrf"][0][1],
                "mlt": data[f"{probe}_mlt"][0][1],
                "mep_ele_flux_tel0_0": data[f"{probe}_mep_ele_flux_tel0"][0][1][0],
                "mep_ele_flux_tel90_0": data[f"{probe}_mep_ele_flux_tel90"][0][1][0],
                "mep_pro_flux_tel0_0": data[f"{probe}_mep_pro_flux_tel0"][0][1][0],
                "mep_pro_flux_tel90_0": data[f"{probe}_mep_pro_flux_tel90"][0][1][0],
            }
        ).set_index("time")

        df.attrs["name"] = probe
        # Calculate the ratio of parallel over perpendicular fluxes
        df["mep_pro_flux_ratio"] = (
            df["mep_pro_flux_tel0_0"] / df["mep_pro_flux_tel90_0"]
        )
        df["mep_ele_flux_ratio"] = (
            df["mep_ele_flux_tel0_0"] / df["mep_ele_flux_tel90_0"]
        )
        return df.replace([np.inf, -np.inf], np.nan).dropna()  # drop NaN value

# file_list = glob.glob("./data/*.sav")
# dss = list(map(lambda x: read_sav(x,xrrary=True), file_list))
# ds = dss[0]


def mep_pro_plot(ds, time_slice, x="lgrf", backend = "proplot", save=False):
    
    start_time = pd.to_datetime(time_slice[0]).strftime("%H:%M")
    end_time = pd.to_datetime(time_slice[1]).strftime("%H:%M")

    match backend:
        case "hvplot":
            time_slice = slice(*time_slice)

            plot01 = (
                ds["mep_pro_flux_tel0"]
                .sel(time=time_slice, proton_energy=0)
                .hvplot(x=x, logy=True)
            )

            plot02 = (
                ds["mep_pro_flux_tel90"]
                .sel(time=time_slice, proton_energy=0)
                .hvplot(x=x, logy=True)
            )

            plot03 = (
                ds["mep_pro_flux_tel0"]
                .sel(time=time_slice, proton_energy=1)
                .hvplot(x=x, logy=True)
            )

            plot04 = (
                ds["mep_pro_flux_tel90"]
                .sel(time=time_slice, proton_energy=1)
                .hvplot(x=x, logy=True)
            )

            return (plot01 * plot02 + plot03 * plot04).cols(1)
        
        case "proplot":
            # create figure with width/height ratio 2:1
            fig = pplt.figure(space=0, refaspect=2,refwidth=5)
            axs = fig.subplots(nrows=2, ncols=1)

            temp_ds = ds.sel(time=slice(*time_slice), proton_energy=1)
            cycle = pplt.Cycle(['black', 'red'])
            axs[0].plot(x="l_igrf",y="mep_pro_flux_tel0",data=temp_ds,
                        cycle=cycle)
            axs[0].plot(x="l_igrf",y="mep_pro_flux_tel90",data=temp_ds,
                       cycle=cycle)

            temp_ds = ds.sel(time=slice(*time_slice), proton_energy=0)
            axs[1].plot(x="l_igrf",y="mep_pro_flux_tel0",data=temp_ds,label=r"$\parallel$", cycle=cycle)
            axs[1].plot(x="l_igrf",y="mep_pro_flux_tel90",data=temp_ds,
                        label=r"$\perp$", cycle=cycle)

            axs.format(
                leftlabels=('80 - 240 keV', '30-80 keV'),
                bottomlabels=(f'{ds.probe}',)
            )
            axs[0].text(
                0.72, 0.92, f'Time : {start_time} - {end_time}',
                transform="axes", size='large'
            )
            axs[0].text(
                0.82, 0.82, f'<MLT> = {int(temp_ds.mlt.mean())}',
                transform="axes", size='large'
            )
            axs[1].legend(frame=False,ncol=1)
            # y axis name to energy range
            axs[0].set_ylabel(r"MEPED Proton FLux [$\#/cm^2/s/sr/keV$]")
            axs[0].set_yscale("log")
            axs[1].set_yscale("log")
            
            if save:
                fname = f"../figures/{ds.probe}_mep_pro_flux_{start_time} - {end_time}"
                fig.save(fname+'.png')
                fig.save(fname+'.svg')
            
            return fig


# time_slices = [
#     ["2021-04-17 01:15:00", "2021-04-17 01:20:00"],
#     ["2021-04-17 02:57:00", "2021-04-17 03:03:00"],
#     ["2021-04-17 04:38:00", "2021-04-17 04:44:00"],
# ]

# [mep_pro_plot(dss[0],time_slice) for time_slice in time_slices]


# time_slices = [
#     ["2021-04-17 01:47:00", "2021-04-17 01:52:00"],
#     ["2021-04-17 03:30:00", "2021-04-17 03:36:00"],
#     ["2021-04-17 05:12:00", "2021-04-17 05:18:00"],
# ]

# [mep_pro_plot(dss[4], time_slice) for time_slice in time_slices]


# %%
def format_energy(list):
    return [f'{erg:0.0f} keV' for erg in list]

def calculate_avg_diff_fluxes(trange, probe='16'):
    files = pyspedas.goes.mpsh(trange=trange, probe = probe, downloadonly=True)
    ds = xr.open_mfdataset(files, combine='nested', concat_dim='record_number')
    avg_proton_diff_chans = ds['DiffProtonEffectiveEnergy'].mean(dim=['record_number','telescopes'])
    avg_electron_diff_chans = ds['DiffElectronEffectiveEnergy'].mean(dim=['record_number','telescopes'])
    return avg_proton_diff_chans, avg_electron_diff_chans

def format_flux_data(flux, flux_type, avg_diff_chans, probe):
    """
    Format the flux data.

    Parameters:
    flux : Dataset
        The flux dataset to format.
    flux_type : str
        The type of the flux ('proton' or 'electron').
    avg_diff_chans : DataArray
        The average differential channels.
    probe : str
        The probe used.

    Returns:
    flux_formatted : Dataset
        The formatted flux dataset.
    """
    flux_formatted = flux.rename({
        'v_dim': 'telescopes',
        'v2_dim': f'{flux_type}_diff_chans',
    }).assign_coords(
        {f'{flux_type}_diff_chans': format_energy(avg_diff_chans.values)}
    )
    flux_formatted.attrs['units'] = r'$1/cm^2/s/sr/keV$'
    flux_formatted.attrs['long_name'] = f'G{probe} {flux_type.capitalize()} Flux'
    return flux_formatted

def goes_part_prod(trange, probe='16', flux_type='all'):
    """
    Data products from the GOES satellite.
    """
    # Time-averaged electron fluxes in several differential channels between 50 and 4,000 keV
    # Time-averaged proton fluxes in several differential channels between 80 and 10,000 keV

    avg_proton_diff_chans, avg_electron_diff_chans = calculate_avg_diff_fluxes(trange, probe)

    pyspedas.goes.mpsh(trange=trange, probe = probe, time_clip=True, prefix=f'g{probe}_')

    flux_data = []
    if flux_type in ['proton', 'all']:
        goes_pflux = get_data(f'g{probe}_AvgDiffProtonFlux', xarray=True)
        goes_pflux = format_flux_data(goes_pflux, 'proton', avg_proton_diff_chans, probe)
        flux_data.append(goes_pflux)

    if flux_type in ['electron', 'all']:
        goes_eflux = get_data(f'g{probe}_AvgDiffElectronFlux', xarray=True)
        goes_eflux = format_flux_data(goes_eflux, 'electron', avg_electron_diff_chans, probe)
        flux_data.append(goes_eflux)

    return tuple(flux_data)

def geos_mag_prod(trange, probe='16'):
    pyspedas.goes.mag(trange=trange, probe=probe, time_clip=True, prefix=f'g{probe}_')

    goes_mag = get_data(f'g{probe}_b_gsm', xarray=True)
    goes_mag.attrs['units'] = 'nT'
    goes_mag.attrs['long_name'] = f'G{probe} B'

    goes_mag = goes_mag.assign_coords(
        v_dim= [r"$B_x$",r"$B_y$",r"$B_z$"]
        )
    return goes_mag
