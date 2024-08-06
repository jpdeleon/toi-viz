#!/usr/bin/env python
import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import requests
import astropy.units as u
from astropy.coordinates import Distance
from pathlib import Path
import itertools
from astropy.coordinates import SkyCoord

DATA_PATH = '../data/'

def get_tois(
    clobber=False,
    outdir=DATA_PATH,
    verbose=True,
    remove_FP=True,
    remove_known_planets=False,
    add_FPP=False,
):
    """Download TOI list from TESS Alert/TOI Release.

    Parameters
    ----------
    clobber : bool
        re-download table and save as csv file
    outdir : str
        download directory location
    verbose : bool
        print texts

    Returns
    -------
    d : pandas.DataFrame
        TOI table as dataframe
    """
    dl_link = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
    fp = os.path.join(outdir, "TOIs.csv")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not os.path.exists(fp) or clobber:
        d = pd.read_csv(dl_link)  # , dtype={'RA': float, 'Dec': float})
        msg = f"Downloading {dl_link}\n"
        if add_FPP:
            fp2 = os.path.join(outdir, "Giacalone2020/tab4.txt")
            classified = ascii.read(fp2).to_pandas()
            fp3 = os.path.join(outdir, "Giacalone2020/tab5.txt")
            unclassified = ascii.read(fp3).to_pandas()
            fpp = pd.concat(
                [
                    classified[["TOI", "FPP-2m", "FPP-30m"]],
                    unclassified[["TOI", "FPP"]],
                ],
                sort=True,
            )
            d = pd.merge(d, fpp, how="outer").drop_duplicates()
        d.to_csv(fp, index=False)
        msg += f"Saved: {fp}\n"
    else:
        d = pd.read_csv(fp).drop_duplicates()
        msg = f"Loaded: {fp}\n"

    # remove False Positives
    if remove_FP:
        idx = (d["TFOPWG Disposition"] != "FP") & (d["TFOPWG Disposition"] != "FA") & ~(d.Comments.str.contains("SB2").astype(bool))
        d = d[idx]
        msg += "TOIs with TFPWG disposition==`False Positive` & `False Alarm` are removed.\n"
        msg += "TOIs with `SB2` in Comments are removed.\n"
    if remove_known_planets:
        idx = (d["TFOPWG Disposition"] != "KP") & (d["TFOPWG Disposition"] != "VP") & (d["TFOPWG Disposition"] != "CP")
        d = d[idx]
        msg += "TOIs with TFPWG disposition==`Known Planet` & `Validated Planet` are removed.\n"
    
    if verbose:
        print(msg)
    return d.sort_values("TOI")
    
def plot_xyz_uvw(
    df,
    target_gaiaid=None,
    match_id=True,
    df_target=None,
    target_label=None,
    target_color="w",
    target_size=50,
    color="teff_val",
    marker="o",
    verbose=True,
    figsize=(12, 8),
    cmap="viridis",
):
    """
    Plot 3D position in galactocentric (xyz) frame
    and proper motion with radial velocity in galactic cartesian velocities
    (UVW) frame

    Parameters
    ----------
    df : pandas.DataFrame
        contains ra, dec, parallax, pmra, pmdec, rv columns
    target_gaiaid : int
        target gaia DR2 id
    df_target : pandas.Series
        target's gaia parameters

    Note: U is positive towards the direction of the Galactic center (GC);
    V is positive for a star with the same rotational direction as the Sun going around the galaxy,
    with 0 at the same rotation as sources at the Sun’s distance,
    and W positive towards the north Galactic pole

    U,V,W can be converted to Local Standard of Rest (LSR) by subtracting V = 238 km/s,
    the adopted rotation velocity at the position of the Sun from Marchetti et al. (2018).

    See also https://arxiv.org/pdf/1707.00697.pdf which estimates Sun's
    (U,V,W) = (9.03, 255.26, 7.001)

    See also https://arxiv.org/pdf/1804.10607.pdf for modeling Gaia DR2 in 6D
    """
    assert len(df) > 0, "df is empty"
    fig, axs = pl.subplots(2, 3, figsize=figsize, constrained_layout=True)
    ax = axs.flatten()

    errmsg = f"color={color} not in {df.columns}"
    assert color in df.columns, errmsg
    if not np.all(df.columns.isin("X Y Z U V W".split())):
        df = get_transformed_coord(df, frame="galactocentric", verbose=verbose)
    if df_target is not None:
        df_target = get_transformed_coord(
            pd.DataFrame(df_target).T, frame="galactocentric"
        )
    cname = df.Name.unique()[0]
    kwargs = {"marker": r"$\star$",
              "c": target_color,
              "ec": "k",
              "s": target_size,
              "label": target_label,
              "zorder": 100
             }
    n = 0
    for (i, j) in itertools.combinations(["X", "Y", "Z"], r=2):
        if target_gaiaid is not None:
            idx = df.GaiaDR3.astype(int).isin([target_gaiaid])
            if match_id:
                errmsg = f"{cname} does not contain the target gaia id [{target_gaiaid}]"
                assert sum(idx) > 0, errmsg
                ax[n].scatter(
                    df.loc[idx, i],
                    df.loc[idx, j],
                    **kwargs                
                )
            else:
                assert df_target is not None, "provide df_target"
                ax[n].scatter(
                    df_target[i],
                    df_target[j],
                    **kwargs 
                )
        c = df[color] if color is not None else None
        cbar = ax[n].scatter(df[i], df[j], c=c, marker=marker, cmap=cmap)
        ax[n].set_xlabel(i + " [pc]")
        ax[n].set_ylabel(j + " [pc]")
        n += 1

    n = 3
    for (i, j) in itertools.combinations(["U", "V", "W"], r=2):
        if target_gaiaid is not None:
            idx = df.GaiaDR3.astype(int).isin([target_gaiaid])
            if match_id:
                errmsg = f"Given cluster does not contain the target gaia id [{target_gaiaid}]"
                assert sum(idx) > 0, errmsg
                ax[n].scatter(
                    df.loc[idx, i],
                    df.loc[idx, j],
                    **kwargs 
                )
            else:
                ax[n].scatter(
                    df_target[i],
                    df_target[j],
                    **kwargs 
                )
        # _ = df.plot.scatter(x=i, y=j, c=color, marker=marker, ax=ax[n], cmap=cmap)
        c = df[color] if color is not None else None
        cbar = ax[n].scatter(df[i], df[j], c=c, marker=marker, cmap=cmap)
        if (color is not None) and (n == 5):
            # show last colorbar only only
            fig.colorbar(cbar, ax=ax[n], label=color)
        ax[n].set_xlabel(i + " [km/s]")
        ax[n].set_ylabel(j + " [km/s]")
        n += 1
    if target_label:
        ax[0].legend()
    fig.suptitle(cname)
    return fig

def get_transformed_coord(df, frame="galactocentric", verbose=True):
    """
    Parameters
    ----------
    df : pandas.DataFrame
        catalog with complete kinematics parameters
    frame : str
        frame conversion

    Returns
    -------
    df : pandas.DataFrame
        catalog with transformed coordinates appended in columns

    Note
    ----
    Assumes galactic center distance distance of 8.1 kpc based on the GRAVITY
    collaboration, and a solar height of z_sun=0 pc.
    See also:
    http://learn.astropy.org/rst-tutorials/gaia-galactic-orbits.html?highlight=filtertutorials
    """
    assert len(df) > 0, "df is empty"
    if np.any(df["Plx"] < 0):
        # retain non-negative parallaxes including nan
        df = df[(df["Plx"] >= 0) | (df["parallax"].isnull())]
        if verbose:
            print("Some parallaxes are negative!")
            print("These are removed for the meantime.")
            print("For proper treatment, see:")
            print("https://arxiv.org/pdf/1804.09366.pdf\n")
    errmsg = f"RV is not in {df.columns}"
    assert df.columns.isin(["RV"]).any(), errmsg
    df2 = df.copy()
    icrs = SkyCoord(
        ra=df2["RA_ICRS"].values * u.deg,
        dec=df2["DE_ICRS"].values * u.deg,
        distance=Distance(parallax=df2["Plx"].values * u.mas),
        radial_velocity=df2["RV"].values * u.km / u.s,
        pm_ra_cosdec=df2["pmRA"].values * u.mas / u.yr,
        pm_dec=df2["pmDE"].values * u.mas / u.yr,
        frame="fk5",
        # equinox="J2015.5",
    )
    # transform to galactocentric frame
    if frame == "galactocentric":
        # xyz = icrs.transform_to(
        #     Galactocentric(z_sun=0 * u.pc, galcen_distance=8.1 * u.kpc)
        # )
        xyz = icrs.galactocentric
        df2["X"] = xyz.x.copy()
        df2["Y"] = xyz.y.copy()
        df2["Z"] = xyz.z.copy()
        df2["U"] = xyz.v_x.copy()
        df2["V"] = xyz.v_y.copy()
        df2["W"] = xyz.v_z.copy()

    elif frame == "galactic":
        # transform to galactic frame
        # gal = icrs.transform_to("galactic")
        gal = icrs.galactic
        df2["gal_l"] = gal.l.deg.copy()
        df2["gal_b"] = gal.b.deg.copy()
        df2["gal_pm_b"] = gal.pm_b.copy()
        df2["gal_pm_l_cosb"] = gal.pm_l_cosb.copy()
    else:
        raise ValueError(f"frame={frame} is unavailable")
    return df2


def get_quadrant(lat,lon):
    """
    quadrant in Perren+2023
    https://ucc.ar/database/
    """
    if -90 <= lat < 0:
        if 0 <= lon < 90:
            return "Q1N"
        elif 90 <= lon < 180:
            return "Q2N"
        elif 180 <= lon < 270:
            return "Q3N"
        elif 270 <= lon < 360:
            return "Q4N"
    elif 0 <= lat < 90:
        if 0 <= lon < 90:
            return "Q1P"
        elif 90 <= lon < 180:
            return "Q2P"
        elif 180 <= lon < 270:
            return "Q3P"
        elif 270 <= lon < 360:
            return "Q4P"
    else:
        return "Invalid latitude or longitude values"

