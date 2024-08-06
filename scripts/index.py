import os
import sys
from urllib.parse import urlencode
from pathlib import Path
from pprint import pprint

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import Angle, Distance
from astropy.table import Table
from astroquery.vizier import Vizier
import altair as alt

outpath = '../index.html'
DATA_PATH = '../data/'
EXOFOP_URL = 'https://exofop.ipac.caltech.edu/tess/target.php?'

def get_tois(
    clobber=False,
    outdir=DATA_PATH,
    verbose=False,
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
    else:
        d = pd.read_csv(fp).drop_duplicates()
        msg = f"Loaded: {fp}\n"
    assert len(d) > 1000, f"{fp} likely has been overwritten!"

    # remove False Positives
    if remove_FP:
        d = d[d["TFOPWG Disposition"] != "FP"]
        msg += "TOIs with TFPWG disposition==FP are removed.\n"
    if remove_known_planets:
        planet_keys = [
            "HD",
            "GJ",
            "LHS",
            "XO",
            "Pi Men" "WASP",
            "SWASP",
            "HAT",
            "HATS",
            "KELT",
            "TrES",
            "QATAR",
            "CoRoT",
            "K2",  # , "EPIC"
            "Kepler",  # "KOI"
        ]
        keys = []
        for key in planet_keys:
            idx = ~np.array(
                d["Comments"].str.contains(key).tolist(), dtype=bool
            )
            d = d[idx]
            if idx.sum() > 0:
                keys.append(key)
        msg += f"{keys} planets are removed.\n"
    msg += f"Saved: {fp}\n"
    if verbose:
        print(msg)
    return d.sort_values("TOI")

VIZIER_KEYS_CLUSTER_CATALOG = {
    #regularly updated in zenodo
    #"Perren": "https://zenodo.org/records/8250524", 
    # 3794 open clusters parameters (Hao+, 2021)
    "Hao2022": "J/A+A/660/A4",
    # c.f. Evolution of the local spiral structure of the Milky Way revealedby open clusters
    "Hao2021": "J/A+A/652/A102",
    # 1656 new star clusters found in Gaia EDR3
    # "He2022c": "https://ui.adsabs.harvard.edu/abs/2022arXiv220908504H/abstract",
    # 886 Clusters within 1.2 kpc of the Sun
    "He2022b": "J/ApJS/262/7",
    # 541 new open cluster candidates
    "He2022a": "J/ApJS/260/8",
    # 628 new open clusters found with OCfinder
    "CastroGinard2022": "J/A+A/661/A118",
    # 570 new open clusters in the Galactic disc
    "CastroGinard2020": "J/A+A/635/A45",
    # 1481 clusters and their members
    "CantatGaudin2020": "J/A+A/633/A99",
    # open clusters in the Galactic anticenter
    "CastroGinard2019": "J/A+A/627/A35",
    #
    "CantatGaudin2018": "J/A+A/618/A93",
    # HRD of Gaia DR2
}

CATALOG_LIST = [key for key in VIZIER_KEYS_CLUSTER_CATALOG.keys()]

class CatalogDownloader:
    """download tables from vizier

    Attributes
    ----------
    tables : astroquery.utils.TableList
        collection of astropy.table.Table downloaded from vizier
    """

    def __init__(
        self, catalog_name, data_loc=DATA_PATH, verbose=True, clobber=False
    ):
        self.catalog_name = catalog_name
        self.catalog_dict = VIZIER_KEYS_CLUSTER_CATALOG
        self.verbose = verbose
        self.clobber = clobber
        if not Path(data_loc).exists():
            Path(data_loc).mkdir()
        self.data_loc = Path(data_loc, self.catalog_name)
        self.tables = None

    def get_tables_from_vizier(self, row_limit=50, save=False, clobber=None):
        """row_limit-1 to download all rows"""
        clobber = self.clobber if clobber is None else clobber
        if row_limit == -1:
            msg = "Downloading all tables in "
        else:
            msg = f"Downloading the first {row_limit} rows of each table "
        msg += f"{self.catalog_dict[self.catalog_name]} from vizier."
        if self.verbose:
            print(msg)
        # set row limit
        Vizier.ROW_LIMIT = row_limit

        tables = Vizier.get_catalogs(self.catalog_dict[self.catalog_name])
        errmsg = "No data returned from Vizier."
        assert tables is not None, errmsg
        self.tables = tables

        if self.verbose:
            pprint({k: tables[k]._meta["description"] for k in tables.keys()})

        if save:
            self.save_tables(clobber=clobber)
        return tables

    def save_tables(self, clobber=None):
        errmsg = "No tables to save"
        assert self.tables is not None, errmsg
        clobber = self.clobber if clobber is None else clobber

        if not self.data_loc.exists():
            self.data_loc.mkdir()

        for n, table in enumerate(self.tables):
            fp = Path(self.data_loc, f"{self.catalog_name}_tab{n}.txt")
            if not fp.exists() or clobber:
                table.write(fp, format="ascii")
                if self.verbose:
                    print(f"Saved: {fp}")
            else:
                print("Set clobber=True to overwrite.")

    def get_vizier_url(self, catalog_name=None):
        if catalog_name is None:
            catalog_name = self.catalog_name
        base_url = "https://vizier.u-strasbg.fr/viz-bin/VizieR?-source="
        vizier_key = self.catalog_dict[catalog_name]
        return base_url + vizier_key

    def __repr__(self):
        """Override to print a readable string representation of class
        """
        included_args = ["catalog_name", "cluster_name"]
        args = []
        for key in self.__dict__.keys():
            val = self.__dict__.get(key)
            if key in included_args:
                if val is not None:
                    args.append(f"{key}={val}")
        args = ", ".join(args)
        return f"{type(self).__name__}({args})"

class ClusterCatalog(CatalogDownloader):
    def __init__(
        self,
        catalog_name="CantatGaudin2020",
        verbose=True,
        clobber=False,
        data_loc=DATA_PATH,
    ):
        super().__init__(
            catalog_name=catalog_name,
            data_loc=data_loc,
            verbose=verbose,
            clobber=clobber,
        )
        """Initialize the catalog

        Attributes
        ----------
        data_loc : str
            data directory
        all_members: pd.DataFrame
            list of all members in catalog
        all_clusters : pd.DataFrame
            list of all clusters in catalog

        Note:
        setting `all_members` as a method (as opposed to attribute)
        seems not
        """
        self.catalog_list = CATALOG_LIST
        self.all_clusters = None  # self.query_catalog(return_members=False)
        self.all_members = None

        if self.data_loc.exists():  # & len(files)<2:
            if self.clobber:
                _ = self.get_tables_from_vizier(
                    row_limit=-1, save=True, clobber=self.clobber
                )
        else:
            _ = self.get_tables_from_vizier(
                row_limit=-1, save=True, clobber=self.clobber
            )
        if self.verbose:
                print("Vizier URL:", self.get_vizier_url())

    def query_catalog(self, name=None, return_members=False, **kwargs):
        """Query catalogs

        Parameters
        ----------
        name : str
            catalog name; see `self.catalog_list`
        return_members : bool
            return parameters for all members instead of the default

        Returns
        -------
        df : pandas.DataFrame
            dataframe parameters of the cluster or its individual members
        Note:
        1. See self.vizier_url() for details
        2. Use the following:
        if np.any(df["parallax"] < 0):
            df = df[(df["parallax"] >= 0) | (df["parallax"].isnull())]
            if verbose:
                print("Some parallaxes are negative!")
                print("These are removed for the meantime.")
                print("For proper treatment, see:")
                print("https://arxiv.org/pdf/1804.09366.pdf\n")

        FIXME: self.all_clusters and self.all_members are repeated each if else block
        """
        self.catalog_name = name if name is not None else self.catalog_name
        if self.verbose:
            print(f"Using {self.catalog_name} catalog.")

        if self.catalog_name == "Hao2022":
            if return_members:
                df_mem = self.get_members_Hao2022()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_Hao2022()
                self.all_clusters = df
                return df
        elif self.catalog_name == "He2022a":
            if return_members:
                df_mem = self.get_members_He2022a()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_He2022a()
                self.all_clusters = df
                return df
        elif self.catalog_name == "He2022b":
            if return_members:
                df_mem = self.get_members_He2022b()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_He2022b()
                self.all_clusters = df
                return df
        elif self.catalog_name == "CastroGinard2022":
            if return_members:
                df_mem = self.get_members_CastroGinard2022()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_CastroGinard2022()
                self.all_clusters = df
                return df
        elif self.catalog_name == "Bouma2019":
            if return_members:
                df_mem = self.get_members_Bouma2019()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_Bouma2019()
                self.all_clusters = df
                return df
        elif self.catalog_name == "CantatGaudin2020":
            if return_members:
                df_mem = self.get_members_CantatGaudin2020()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_CantatGaudin2020()
                self.all_clusters = df
                return df
        elif self.catalog_name == "CastroGinard2020":
            if return_members:
                df_mem = self.get_members_CastroGinard2020()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_CastroGinard2020()
                self.all_clusters = df
                return df
        elif self.catalog_name == "CastroGinard2019":
            if return_members:
                df_mem = self.get_members_CastroGinard2019()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_CastroGinard2019()
                self.all_clusters = df
                return df
        elif self.catalog_name == "CantatGaudin2018":
            if return_members:
                df_mem = self.get_members_CantatGaudin2018()
                self.all_members = df_mem
                return df_mem
            else:
                df = self.get_clusters_CantatGaudin2018()
                self.all_clusters = df
                return df

    def get_clusters_CantatGaudin2020(self):
        """Cantat-Gaudin et al. 2020:
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        tab = Table.read(fp, format="ascii")
        df = tab.to_pandas()
        df = _decode_n_drop(df, ["SimbadName", "dmode_01", "dmode-01"])
        df = df.rename(
            columns={
                "RA_ICRS": "raJ2015",
                "DE_ICRS": "decJ2015",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "Plx": "parallax",
                "dmode": "distance",
                "pmRA": "pmra",
                "pmDE": "pmdec",
                "N": "Nstars",
            }
        )
        return df

    def get_members_CantatGaudin2020(self):
        """Cantat-Gaudin et al. 2020:
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        tab = Table.read(fp, format="ascii")
        df = tab.to_pandas()
        df = df.applymap(
            lambda x: x.decode("ascii") if isinstance(x, bytes) else x
        )
        df = df.rename(
            columns={
                "RA_ICRS": "raJ2015",
                "DE_ICRS": "decJ2015",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "Source": "source_id",  # 'RV':'radial_velocity',
                "pmRA": "pmra",
                "pmDE": "pmdec",
                "Plx": "parallax",
                "Gmag": "phot_g_mean_mag",
                "BP-RP": "bp_rp",
            }
        )
        return df

    def get_clusters_Hao2022(self):
        """Hao+2022: Gaia EDR3 new Galactic open clusters
        https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract

        tab1: (c)Mean parameters for the reported open clusters (704 rows)

        c.f. Hao+2021: Evolution of the local spiral structure of the Milky Way revealed by open clusters
        Parameters for 3794 clusters based on the Gaia EDR3 (3794 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                "age": "log10_age",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "e_plx": "e_parallax",
                "pmRA": "pmra",
                "e_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "e_pmDE": "e_pmdec",
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_members_Hao2022(self):
        """Hao+2022:
        https://ui.adsabs.harvard.edu/abs/2022A%26A...660A...4H/abstract

        tab: (c)Members for the reported open clusters (19425 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                # "GaiaEDR3": "source_id_edr3",
                "GaiaEDR3": "source_id",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "e_plx": "e_parallax",
                "pmRA": "pmra",
                "e_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "e_pmDE": "e_pmdec",
                "RV": "radial_velocity",
                "e_RV": "e_radial_velocity",
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_clusters_He2022b(self):
        """He+2022b: 886 Clusters within 1.2 kpc of the Sun
        https://ui.adsabs.harvard.edu/abs/2022ApJS..262....7H/abstract

        tab1: Parameters for the 886 objects (886 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                "logAge": "log10_age",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "s_plx": "e_parallax",
                "pmRA": "pmra",
                "e_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "e_pmDE": "e_pmdec",
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_members_He2022b(self):
        """He+2022b: 886 Clusters within 1.2 kpc of the Sun
        https://ui.adsabs.harvard.edu/abs/2022yCat..22620007H/abstract

        tab1: Gaia EDR3 parameters of the member stars (134192 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                # "GaiaEDR3": "source_id_edr3",
                "GaiaEDR3": "source_id",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "e_plx": "e_parallax",
                "pmRA": "pmra",
                "e_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "e_pmDE": "e_pmdec",
                # "RV": "radial_velocity",
                # "e_RV": "e_radial_velocity"
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_clusters_He2022a(self):
        """He+2022: New Open-cluster Candidates Found in the Galactic Disk Using Gaia DR2/EDR3 Data
        https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract


        tab1: Parameters of median astrometric values and isochrone fits for 541 new open cluster candidates (541 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                "logAge": "log10_age",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "Plx": "parallax",
                "s_Plx": "e_parallax",
                "pmRA": "pmra",
                "s_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "s_pmDE": "e_pmdec",
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_members_He2022a(self):
        """He+2022: New Open-cluster Candidates Found in the Galactic Disk Using Gaia DR2/EDR3 Data
        https://ui.adsabs.harvard.edu/abs/2022ApJS..260....8H/abstract

        tab2: Member stars for 541 new OCCs (66468 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                # "GaiaEDR3": "source_id_edr3",
                "GaiaEDR3": "source_id",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "Plx": "parallax",
                "e_Plx": "e_parallax",
                "pmRA": "pmra",
                "e_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "e_pmDE": "e_pmdec",
                "RV": "radial_velocity",
                "e_RV": "e_radial_velocity",
            }
        ).reset_index(drop=True)
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_clusters_CastroGinard2022(self):
        """Castro-Ginard et al. 2022: Hunting for open clusters in Gaia EDR3: 628 new open clusters found with OCfinder
        https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract

        tab1: Mean parameters for the reported UBC clusters (628 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                "RA_ICRS": "raJ2015",
                "s_RA_ICRS": "e_raJ2015",
                "DE_ICRS": "decJ2015",
                "s_DE_ICRS": "e_decJ2015",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "s_plx": "e_parallax",
                "pmRA": "pmra",
                "s_pmRA": "e_pmra",
                "pmDE": "pmdec",
                "s_pmDE": "e_pmdec",
                "RV": "radial_velocity",
                "s_RV": "e_radial_velocity",
            }
        )
        # add distance
        df["distance"] = Distance(parallax=df.parallax.values * u.mas).pc
        return df

    def get_members_CastroGinard2022(self):
        """Castro-Ginard et al. 2022:
        https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.118C/abstract

        tab2: Members for the reported UBC clusters (25466 rows)
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        df = _read_clean_df(fp)
        df = df.rename(
            columns={
                # "GaiaEDR3": "source_id_edr3",
                "GaiaEDR3": "source_id",
                "RA_ICRS": "raJ2015",
                "DE_ICRS": "decJ2015",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "plx": "parallax",
                "pmRA": "pmra",
                "pmDE": "pmdec",
                "Source": "source_id",
            }
        )
        return df

    def get_clusters_CantatGaudin2018(self):
        """Cantat-Gaudin et al. 2018:
        http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab0.txt")
        tab = Table.read(fp, format="ascii")
        df = tab.to_pandas()
        df = _decode_n_drop(df, ["SimbadName", "dmode_01", "dmode-01"])
        df = df.rename(
            columns={
                "RAJ2000": "ra",
                "DEJ2000": "dec",
                "dmode": "distance",
                "pmRA": "pmra",
                "pmDE": "pmdec",
                "plx": "parallax",
            }
        )
        return df

    def get_members_CantatGaudin2018(self):
        """Cantat-Gaudin et al. 2018:
        http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93
        """
        fp = Path(self.data_loc, f"{self.catalog_name}_tab1.txt")
        tab = Table.read(fp, format="ascii")
        df = tab.to_pandas()
        df = df.applymap(
            lambda x: x.decode("ascii") if isinstance(x, bytes) else x
        )
        df = df.rename(
            columns={
                "RA_ICRS": "raJ2015",
                "DE_ICRS": "decJ2015",
                "_RA.icrs": "ra",
                "_DE.icrs": "dec",
                "Source": "source_id",
                "o_Gmag": "phot_g_n_obs",
                "pmRA": "pmra",
                "pmDE": "pmdec",
                "plx": "parallax",
                # 'Gmag': 'phot_g_mean_mag',
                "BP-RP": "bp_rp",
            }
        )
        return df

def _read_clean_df(fp):
    tab = Table.read(fp, format="ascii")
    df = tab.to_pandas()
    df = df.applymap(
        lambda x: x.decode("ascii") if isinstance(x, bytes) else x
    )
    return df

def _decode_n_drop(df, columns):
    """
    columns : list of columns to drop
    """
    # bytes to str
    df = df.applymap(
        lambda x: x.decode("ascii") if isinstance(x, bytes) else x
    )
    # remove columns
    df = df.drop(columns, axis=1)
    return df

def make_exofop_url(ticid):
    return EXOFOP_URL + urlencode({'id': ticid})

def plot_catalog(df, width=800, height=400):
    chart = (
        alt.Chart(df)
        .mark_point(color="red", clip=True)
        .encode(
            x='ra:Q',
            y='dec:Q',
            # x=alt.X(
            #     "ra:Q",
            #     axis=alt.Axis(title="RA"),
            #     scale=alt.Scale(domain=(0, 360)),
            # ),
            # y=alt.Y(
            #     "dec:Q",
            #     axis=alt.Axis(title="Dec"),
            #     scale=alt.Scale(domain=(-90, 90)),
            # ),
            tooltip=[
                "Cluster:N",
                "distance:Q",
                "parallax:Q",
                "pmra:Q",
                "pmdec:Q",
                "Nstars:Q",
            ],
        )
        .properties(width=width, height=height)
        .interactive()
        )
    return chart

def plot_tois(df, width=800, height=400):
    chart = (
        alt.Chart(df, title="TOI")
        # .transform_calculate(
        #     url = EXOFOP_URL+"id="+alt.datum.TIC_ID
        # )
        .mark_point(color="black", clip=True)
        .encode(
            x='RA:Q',
            y='Dec:Q',
            # x = alt.X(
            #     "RA:Q",
            #     axis=alt.Axis(title="RA"),
            #     scale=alt.Scale(domain=(0, 360)),
            # ),
            # y = alt.Y(
            #     "Dec:Q",
            #     axis=alt.Axis(title="Dec"),
            #     scale=alt.Scale(domain=(-90, 90)),
            # ),
            href='url',
            tooltip = [
                "TOI:Q",
                "TIC ID:Q",
                "url:N",
                "Stellar Distance (pc):Q",
                "PM RA (mas/yr):Q",
                "PM Dec (mas/yr):Q",
            ],
        )
        .properties(width=width, height=height)
        .interactive()
    )
    return chart

def plot_clusters(df, width=800, height=400):
    chart = (
        alt.Chart(df)
        .mark_circle(clip=True)
        .encode(
            # x="ra:Q",
            # y="dec:Q",
            x = alt.X(
                "ra:Q",
                axis=alt.Axis(title="RA"),
                scale=alt.Scale(domain=(0, 360)),
            ),
            y = alt.Y(
                "dec:Q",
                axis=alt.Axis(title="Dec"),
                scale=alt.Scale(domain=(-90, 90)),
            ),
            color="Cluster:N",
            tooltip=[
                "source_id:O",
                "parallax:Q",
                "pmra:Q",
                "pmdec:Q",
                "phot_g_mean_mag:Q",
            ],
        )
        .properties(width=width, height=height)
        .interactive()
    )
    return chart

def plot_interactive(
        catalog_name="CantatGaudin2020",
        clobber_toi=False,
        min_parallax=1.5,
        thin=10,
        width=800,
        height=400,
    ):
    """show altair plots of TOI and clusters

    Parameters
    ----------
    plx_cut : float
        parallax cut in mas; default=2 mas < 100pc
    thin : integer
        thinning factor to use ony every nth cluster member
    """
    try:
        import altair as alt
    except ModuleNotFoundError:
        print("pip install altair")

    if sys.argv[-1].endswith("json"):
        print("import altair; altair.renderers.enable('notebook')")

    cc = ClusterCatalog(verbose=False)
    df0 = cc.query_catalog(catalog_name=catalog_name, return_members=False)
    df2 = cc.query_catalog(catalog_name=catalog_name, return_members=True)
    idx = df0.parallax >= min_parallax
    df0 = df0.loc[idx]
    df0["distance"] = Distance(parallax=df0["parallax"].values * u.mas).pc
    # plot catalog
    chart0 = plot_catalog(df0, width=width, height=height)

    # get TOI list
    toi = get_tois(verbose=False, clobber=clobber_toi)
    toi["TIC_ID"] = toi["TIC ID"]
    toi["RA"] = Angle(toi["RA"].values, unit="hourangle").deg
    toi["Dec"] = Angle(toi["Dec"].values, unit="deg").deg
    toi["url"] = toi['TIC ID'].apply(make_exofop_url)
    # plot TOI
    chart1 = plot_tois(toi, width=width, height=height)

    # plot cluster members
    idx = df2.parallax >= min_parallax
    df2 = df2.loc[idx]
    # skip other members
    if thin<10:
        alt.data_transformers.disable_max_rows()
    df2 = df2.iloc[::thin, :]
    chart2 = plot_clusters(df2, width=width, height=height)

    return chart2 + chart1 + chart0


if __name__=="__main__":
    import panel as pn

    # title = '## TOI visualization'

    chart = plot_interactive(
                            catalog_name='He2022b',
                            clobber_toi=True,
                            thin=1,
                            min_parallax=1.5,
                            width=1600,
                            height=900,
                            )
    # to allow new browser tab pop-up
    chart['usermeta'] = {
            "embedOptions": {
            'loader': {'target': '_blank'}
            }
        }
    chart.save(outpath)

    # pn.Row(
    #     pn.Column(title, ticker, window)
    # )
    # panel.save(outpath, embed=True)
    print("Saved: ", outpath)
