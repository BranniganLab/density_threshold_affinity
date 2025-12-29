import numpy as np
import matplotlib.pyplot as plt

#TODO: define theta_edges, r_edges, state variables

# ----- Helper functions -----
def clear_artists(artists):
    for a in artists:
        a.remove()
    artists.clear()

def bins_in_region(t0, t1, r0, r1):
    r0, r1 = sorted([r0, r1])
    t0 = t0 % (2*np.pi)
    t1 = t1 % (2*np.pi)

    dt = (t1 - t0) % (2*np.pi)
    if dt <= np.pi:
        def theta_in_bin(t_low, t_high):
            return ((t_low - t0) % (2*np.pi) < dt) or ((t_high - t0) % (2*np.pi) < dt) or ((t0 - t_low) % (2*np.pi) < (t_high - t_low))
    else:
        dt = (t0 - t1) % (2*np.pi)
        def theta_in_bin(t_low, t_high):
            return ((t_low - t1) % (2*np.pi) < dt) or ((t_high - t1) % (2*np.pi) < dt) or ((t1 - t_low) % (2*np.pi) < (t_high - t_low))

    bins = []
    for ti in range(len(theta_edges)-1):
        t_low, t_high = theta_edges[ti], theta_edges[ti+1]
        if theta_in_bin(t_low, t_high):
            for ri in range(len(r_edges)-1):
                r_low, r_high = r_edges[ri], r_edges[ri+1]
                if r_high < r0 or r_low > r1:
                    continue
                bins.append((ti, ri))
    return bins

def draw_outer_edges_circular(mask_bins, color='red', lw=2):
    if not mask_bins:
        return []

    n_r = len(r_edges) - 1
    n_t = len(theta_edges) - 1
    mask = np.zeros((n_r, n_t), dtype=bool)
    for ti, ri in mask_bins:
        mask[ri, ti] = True

    artists = []
    for ri in range(n_r):
        for ti in range(n_t):
            if not mask[ri, ti]:
                continue

            top = ri == n_r-1 or not mask[ri+1, ti]
            bottom = ri == 0 or not mask[ri-1, ti]

            left_ti = (ti - 1) % n_t
            right_ti = (ti + 1) % n_t
            left = not mask[ri, left_ti]
            right = not mask[ri, right_ti]

            t0, t1 = theta_edges[ti], theta_edges[(ti+1)%n_t]
            r0, r1 = r_edges[ri], r_edges[ri+1]

            if top:
                artists.append(ax.plot([t0, t1], [r1, r1], color=color, lw=lw, zorder=10)[0])
            if bottom:
                artists.append(ax.plot([t0, t1], [r0, r0], color=color, lw=lw, zorder=10)[0])
            if left:
                artists.append(ax.plot([t0, t0], [r0, r1], color=color, lw=lw, zorder=10)[0])
            if right:
                artists.append(ax.plot([t1, t1], [r0, r1], color=color, lw=lw, zorder=10)[0])
    return artists

# ----- Event handlers -----
def on_press(event):
    global drag_start, mode, selected_bins
    if event.inaxes != ax or event.button != 1:
        return
    drag_start = (event.xdata, event.ydata)

    if event.key == "shift":
        mode = "add"
    elif event.key == "control":
        mode = "subtract"
    else:
        mode = "replace"
        selected_bins.clear()  # Clear old selection immediately

def on_motion(event):
    if drag_start is None or event.inaxes != ax:
        return
    if event.xdata is None or event.ydata is None:
        return

    clear_artists(hover_artists)
    dx = event.xdata - drag_start[0]
    dy = event.ydata - drag_start[1]
    threshold = 1e-8

    if abs(dx) < threshold and abs(dy) < threshold:
        ti = np.searchsorted(theta_edges, event.xdata % (2*np.pi), side="right") - 1
        ri = np.searchsorted(r_edges, event.ydata, side="right") - 1
        bins = [(ti, ri)] if 0 <= ti < len(theta_edges)-1 and 0 <= ri < len(r_edges)-1 else []
    else:
        bins = bins_in_region(drag_start[0], event.xdata, drag_start[1], event.ydata)

    # Immediate hover display
    if mode == "replace":
        temp_mask = set(bins)
    elif mode == "add":
        temp_mask = selected_bins.union(bins)
    elif mode == "subtract":
        temp_mask = selected_bins.difference(bins)
    else:
        temp_mask = selected_bins.union(bins)

    hover_artists.extend(draw_outer_edges_circular(temp_mask, color='orange', lw=1.5))
    fig.canvas.draw_idle()

def on_release(event):
    global drag_start
    if drag_start is None or event.xdata is None or event.ydata is None:
        drag_start = None
        return

    dx = event.xdata - drag_start[0]
    dy = event.ydata - drag_start[1]
    threshold = 1e-8

    if abs(dx) < threshold and abs(dy) < threshold:
        ti = np.searchsorted(theta_edges, event.xdata % (2*np.pi), side="right") - 1
        ri = np.searchsorted(r_edges, event.ydata, side="right") - 1
        bins = [(ti, ri)] if 0 <= ti < len(theta_edges)-1 and 0 <= ri < len(r_edges)-1 else []
    else:
        bins = bins_in_region(drag_start[0], event.xdata, drag_start[1], event.ydata)

    if mode == "replace":
        selected_bins.update(bins)  # Already cleared on press
    elif mode == "add":
        selected_bins.update(bins)
    elif mode == "subtract":
        selected_bins.difference_update(bins)

    clear_artists(selected_artists)
    selected_artists.extend(draw_outer_edges_circular(selected_bins, color='red', lw=2))
    clear_artists(hover_artists)

    drag_start = None
    fig.canvas.draw_idle()