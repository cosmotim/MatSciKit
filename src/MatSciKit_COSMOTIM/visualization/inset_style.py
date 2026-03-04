"""
Inset-style figure export.

Formats matplotlib figures for inset panels in journal publications.

Translated from plotExportInsetStyle.m
"""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from typing import Optional, Union


def export_inset_figure(fig: Optional[Figure] = None,
                        filename: Union[str, Path] = 'figure_inset',
                        format: str = 'tiff',
                        dpi: int = 600,
                        width: float = 2.0,
                        height: float = 1.5) -> None:
    """
    Export a matplotlib figure with inset-style formatting.

    Applies compact formatting suitable for inset panels in scientific
    publications, with smaller fonts and tighter layout than journal style.

    Parameters
    ----------
    fig : Figure, optional
        Figure handle. Defaults to current figure if None.
    filename : str or Path
        Output filename without extension.
    format : str
        Output format ('tiff', 'png', 'pdf', 'eps'). Default is 'tiff'.
    dpi : int
        Resolution in dots per inch. Default is 600.
    width : float
        Figure width in inches. Default is 2.0.
    height : float
        Figure height in inches. Default is 1.5.
    """
    if fig is None:
        fig = plt.gcf()

    ax = fig.gca() if fig.axes else None

    if ax:
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)

        ax.tick_params(
            which='both',
            direction='in',
            length=3,
            width=1.0,
            labelsize=9,
            top=True,
            right=True
        )

        ax.minorticks_off()

        for item in ([ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(9)
            item.set_fontname('Arial')

        if ax.get_title():
            ax.title.set_fontsize(9)
            ax.title.set_fontname('Arial')

        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_fontsize(8)
                text.set_fontname('Arial')

    fig.patch.set_facecolor('white')
    fig.set_size_inches(width, height)
    fig.tight_layout(pad=0.2)

    filename = Path(filename)
    if not filename.suffix:
        filename = filename.with_suffix(f'.{format}')

    fig.savefig(
        filename,
        dpi=dpi,
        format=format,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none',
        pad_inches=0.03
    )

    print(f'Figure exported as {filename}')
