# Python Package Architecture

This folder contains the core Python implementation for interactive polar-bin selection using Matplotlib.  
The code is designed for **interactive environments** (including Jupyter notebooks with widget backends) and emphasizes **correctness, debuggability, and separation of concerns**.

---

## Overview

The system allows a user to:

- Interactively select polar bins using click-and-drag gestures
- Add to, subtract from, or replace selections using modifier keys
- Drag outside an Axes and still commit the selection correctly
- Work with multiple Axes in the same Figure without state leakage

To support this reliably, the code is organized using a **Model–View–Controller (MVC)** architecture.

---

## Model–View–Controller (MVC)

### Roles

- **Model**  
  Owns domain state and pure geometry rules.  
  Has no knowledge of Matplotlib events or rendering.

- **View**  
  Responsible only for drawing on Matplotlib Axes.  
  Does not encode selection semantics.

- **Controller**  
  Translates user input (mouse + modifiers) into model updates and view updates.  
  Owns all gesture and interaction state.

---

## MVC Class Diagram

```mermaid
classDiagram

class BinSelectionModel:::model
class PolarBinGrid:::model

class PolarBinRenderer:::view
class BinEdge:::view

class SiteSelector:::controller
class SiteSelectorManager:::controller
class SelectorDragState:::controller
class SelectorDrawState:::controller
class SelectionOperation:::controller

<<Model>> BinSelectionModel
<<Model>> PolarBinGrid

<<View>> PolarBinRenderer
<<View>> BinEdge

<<Controller>> SiteSelector
<<Controller>> SiteSelectorManager
<<Controller>> SelectorDragState
<<Controller>> SelectorDrawState
<<Controller>> SelectionOperation

class BinSelectionModel {
    -_bins: Set[(int,int)]
    +snapshot()
    +bins()
    +set(bins)
    +add(bins)
    +remove(bins)
    +clear()
}

class PolarBinGrid {
    +map_coord_to_bin_idx(r, theta)
    +bins_in_region(r0, t0, r1, t1)
    +exposed_edges(bins)
    +bin_in_theta_arc(theta_start, theta_end, bin_start, bin_end)
}

class PolarBinRenderer {
    +draw_edges(edges, plot_kwargs)
}

class SiteSelector {
    +on_press(event)
    +on_motion(event)
    +on_release(event)
}

class SiteSelectorManager {
    +register(selector)
    +set_active(selector)
}

SiteSelector --> PolarBinGrid
SiteSelector --> BinSelectionModel
SiteSelector --> PolarBinRenderer
SiteSelector *-- SelectorDragState
SiteSelector *-- SelectorDrawState
SiteSelector --> SelectionOperation

SiteSelectorManager o-- SiteSelector

classDef model fill:#E8F5E9,stroke:#2E7D32,color:#1B5E20;
classDef view fill:#E3F2FD,stroke:#1565C0,color:#0D47A1;
classDef controller fill:#FFF3E0,stroke:#EF6C00,color:#E65100;
