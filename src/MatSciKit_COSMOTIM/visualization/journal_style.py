"""
Journal-Style Figure Export Module

This module provides functions to format matplotlib figures for journal publication
and export them with appropriate settings.
"""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from typing import Optional, Union


def export_journal_figure(fig: Optional[Figure] = None, 
                          filename: Union[str, Path] = 'figure',
                          format: str = 'tiff',
                          dpi: int = 600,
                          width: float = 4.0,
                          height: float = 3.0) -> None:
    """
    Export a matplotlib figure with journal-style formatting.
    
    Applies professional formatting suitable for scientific publication including:
    - Arial font at 11pt for labels and 10pt for legend
    - 1.5pt line widths for axis spines
    - Inward-facing ticks on all sides
    - Proper figure sizing and high-resolution export
    
    Args:
        fig: Figure handle (defaults to current figure if None)
        filename: Output filename without extension
        format: Output format ('tiff', 'png', 'pdf', 'eps')
        dpi: Resolution in dots per inch (default: 600)
        width: Figure width in inches (default: 4.0)
        height: Figure height in inches (default: 3.0)
    """
    if fig is None:
        fig = plt.gcf()
    
    # Get axes handle
    ax = fig.gca() if fig.axes else None
    
    if ax:
        # Set axis spine properties
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        
        # Set tick parameters
        ax.tick_params(
            which='both',
            direction='in',
            length=4,
            width=1.5,
            labelsize=11,
            top=True,
            right=True
        )
        
        # Disable minor ticks
        ax.minorticks_off()
        
        # Set font properties for axis labels and tick labels
        for item in ([ax.xaxis.label, ax.yaxis.label] + 
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(11)
            item.set_fontname('Arial')
        
        # Set title font if exists
        if ax.get_title():
            ax.title.set_fontsize(11)
            ax.title.set_fontname('Arial')
        
        # Set legend font if exists
        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_fontsize(10)
                text.set_fontname('Arial')
    
    # Set figure properties
    fig.patch.set_facecolor('white')
    fig.set_size_inches(width, height)
    
    # Adjust layout to prevent label clipping
    fig.tight_layout(pad=0.3)
    
    # Handle filename and extension
    filename = Path(filename)
    if not filename.suffix:
        filename = filename.with_suffix(f'.{format}')
    
    # Export the figure
    fig.savefig(
        filename,
        dpi=dpi,
        format=format,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none',
        pad_inches=0.05
    )
    
    print(f'Figure exported as {filename}')


