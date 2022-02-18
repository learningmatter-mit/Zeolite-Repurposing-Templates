import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


LIMS1 = [[50, 400], [-19, 0]]
LIMS2 = [[3.5, 14], [3, 10.5]]

GRIDSPEC_KWS = {
    "width_ratios": (0.45, 0.45, 0.02),
    "hspace": 0.2,
}

CMAP = "inferno_r"
SCATTER_KWS = {
    "cmap": CMAP,
    "linewidths": 0.5,
    "edgecolors": "#505050",
    "s": 60,
}
ANNOT_SCATTER_KWS = {
    **SCATTER_KWS,
    "linewidths": 1.5,
    "edgecolors": "k",
    "s": 200,
}
HEXBIN_KWS = {
    "mincnt": 1,
    "gridsize": 10,
    "reduce_C_function": np.mean,
    "edgecolors": "w",
    "cmap": CMAP,
    "linewidths": (0.0,),
}

FIGSIZE = (6.3, 2.3)


def get_literature_markers(in_literature, default="o"):
    if in_literature == 1.0:
        return "^"
    return default


def mscatter(x, y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers

    ax = ax or plt.gca()
    sc = ax.scatter(x, y, **kw)

    if (m is not None) and (len(m) == len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc


class ShapePlotter:
    def __init__(
        self,
        df,
        figsize=FIGSIZE,
        color_option="Binding (SiO2)",
        gridspec_kws=GRIDSPEC_KWS,
        scatter_kws=SCATTER_KWS,
        hexbin_kws=HEXBIN_KWS,
        annot_scatter_kws=ANNOT_SCATTER_KWS,
        norm=None,
    ):
        self.df = df
        self.figsize = figsize
        self.color_opt = color_option
        self.gridspec_kws = gridspec_kws
        self.scatter_kws = scatter_kws
        self.annot_scatter_kws = annot_scatter_kws
        self.hexbin_kws = hexbin_kws
        self.norm = norm

    def get_figure(self):
        fig, ax_fig = plt.subplots(
            1,
            3,
            figsize=self.figsize,
            gridspec_kw=self.gridspec_kws,
            constrained_layout=True,
        )
        return fig, ax_fig

    def get_color(self, values, cmap):
        norm = plt.Normalize(values.min(), values.max())
        cmap_ = cm.get_cmap(cmap)
        return cmap_(norm(values))

    def get_subset(self, zeolite, lims1, lims2):
        xlim1, ylim1 = lims1
        xlim2, ylim2 = lims2

        subdf = self.df.loc[
            (self.df["Zeolite"] == zeolite)
            & (self.df["Volume (Angstrom3)"] > xlim1[0])
            & (self.df["Volume (Angstrom3)"] < xlim1[1])
            & (self.df["Binding (SiO2)"] > ylim1[0])
            & (self.df["Binding (SiO2)"] < ylim1[1])
        ].sort_values("Binding (SiO2)", ascending=False)
        subdf["color"] = subdf[self.color_opt]
        subdf["markers"] = (
            subdf["In literature?"].apply(get_literature_markers).values.tolist()
        )

        return subdf

    def get_annotation_df(self, subdf, osdas):
        smiles = list(osdas.keys())

        df_ = subdf.copy()
        df_["color"] = self.get_color(
            subdf["color"].values, self.scatter_kws.get("cmap", "inferno_r")
        ).tolist()

        annotdf = df_.loc[subdf["SMILES"].isin(smiles)]

        annotdf["markers"] = (
            annotdf["In literature?"]
            .apply(lambda x: get_literature_markers(x, "s"))
            .values.tolist()
        )
        annotdf["labels"] = [osdas[smi] for smi in annotdf["SMILES"]]

        return annotdf

    def annotate(self, ax, annotdf, xvar, yvar):
        scat = mscatter(
            annotdf[xvar],
            annotdf[yvar],
            ax=ax,
            c=annotdf["color"],
            m=annotdf["markers"],
            **self.annot_scatter_kws,
        )

        for _, row in annotdf.iterrows():
            ax.annotate(
                row["labels"],
                (row[xvar], row[yvar]),
                zorder=3,
                ha="center",
                va="center",
            )

        return scat

    def plot_volume(self, ax, subdf, lims, osdas):
        x = "Volume (Angstrom3)"
        y = "Binding (SiO2)"

        scat = mscatter(
            subdf[x],
            subdf[y],
            ax=ax,
            c=subdf["color"],
            m=subdf["markers"],
            norm=self.norm,
            **self.scatter_kws,
        )

        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_xlim(lims[0])

        if osdas is not None:
            annotdf = self.get_annotation_df(subdf, osdas)
            self.annotate(ax, annotdf, x, y)

        return scat

    def plot_shape(self, ax, subdf, lims, osdas):
        x = "Axis 1 (Angstrom)"
        y = "Axis 2 (Angstrom)"

        LIM_TOL = 0.4
        xlim = lims[0]
        ylim = lims[1]
        newlimx = [xlim[0] + LIM_TOL, xlim[1] - LIM_TOL]
        newlimy = [ylim[0] + LIM_TOL, ylim[1] - LIM_TOL]

        hbin = ax.hexbin(
            subdf[x],
            subdf[y],
            C=subdf[self.color_opt],
            extent=(newlimx + newlimy),
            **self.hexbin_kws,
            norm=self.norm,
        )
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_xlim(lims[0])
        ax.set_ylim(lims[1])

        if osdas is not None:
            annotdf = self.get_annotation_df(subdf, osdas)
            self.annotate(ax, annotdf, x, y)

        return hbin

    def plot_cbar(self, ax, fig, scat):
        cbar = fig.colorbar(scat, cax=ax)
        ax.set_ylabel(self.color_opt)

        return cbar

    def plot(self, zeolite, lims1, lims2, osdas=None):
        fig, ax_fig = self.get_figure()
        subdf = self.get_subset(zeolite, lims1, lims2)

        scat = self.plot_volume(ax_fig[0], subdf, lims1, osdas)
        hbin = self.plot_shape(ax_fig[1], subdf, lims2, osdas)
        cbar = self.plot_cbar(ax_fig[2], fig, scat)

        return fig, ax_fig
