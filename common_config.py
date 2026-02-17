import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica"]
plt.rcParams["font.size"] = 6 
plt.rcParams["svg.fonttype"] = "none"

rsc_lightgreen = "#6EC207"
rsc_lightblue = "#56B4E9"
rsc_orange = "#D55E00"
rsc_grey = "#959595"


def figure_journal_basic():
    sc.set_figure_params(
        dpi=100,
        dpi_save=300,
        figsize=(4, 4),
        frameon=False,
        facecolor=None,
        transparent=True,
        vector_friendly=True,
    )


def figure_jupyter():
    sc.set_figure_params(
        dpi=80,
        dpi_save=150,
        figsize=(8, 8),
        frameon=False,
        facecolor="white",
        transparent=False,
        vector_friendly=True,
    )
