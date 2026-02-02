#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from dta.gui.matplotlib import SiteSelectorManager, SiteSelector


def example_usage():
    """
    Create a demo figure and attach one `SiteSelector` to each of two Axes.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    axes : tuple[matplotlib.axes.Axes, matplotlib.axes.Axes]
        The two polar Axes in the figure.
    manager : SiteSelectorManager
        The manager responsible for routing events to selectors.

    Notes
    -----
    This function is intended for manual exploration and debugging in an
    interactive environment. It creates two polar pcolormesh plots and registers
    a separate selector for each Axes.
    """
    theta_edges = np.linspace(0.0, 2.0 * np.pi, 49)
    r_edges = np.linspace(0.0, 1.0, 21)

    theta, r = np.meshgrid(theta_edges, r_edges)
    data = np.random.rand(20, 48)

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        subplot_kw={"projection": "polar"},
        figsize=(10, 5),
        constrained_layout=True,
    )

    pcm1 = ax1.pcolormesh(
        theta,
        r,
        data,
        shading="auto",
        cmap="viridis",
    )
    ax1.set_title("Polar Data A")
    fig.colorbar(pcm1, ax=ax1, pad=0.1)

    pcm2 = ax2.pcolormesh(
        theta,
        r,
        data,
        shading="auto",
        cmap="plasma",
    )
    ax2.set_title("Polar Data B")
    fig.colorbar(pcm2, ax=ax2, pad=0.1)

    plotting_kwargs = {
        "color": "red",
        "lw": 2.0,
        "zorder": 20,
    }

    selector_a = SiteSelector(
        ax1,
        theta_edges=theta_edges,
        r_edges=r_edges,
        plot_kwargs=plotting_kwargs,
    )

    plotting_kwargs["color"] = "cyan"

    selector_b = SiteSelector(
        ax2,
        theta_edges=theta_edges,
        r_edges=r_edges,
        plot_kwargs=plotting_kwargs,
    )

    manager = SiteSelectorManager(fig)
    manager.register(selector_a, active=True)
    manager.register(selector_b, active=True)

    return fig, (ax1, ax2), manager
