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

%% =========================
%% MVC Stereotypes
%% =========================

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


%% =========================
%% Model
%% =========================
class BinSelectionModel {
    -_bins: Set<(int,int)>
    +snapshot(): FrozenSet<(int,int)>
    +bins(): Set<(int,int)>
    +set(bins: Set<(int,int)>)
    +add(bins: Set<(int,int)>)
    +remove(bins: Set<(int,int)>)
    +clear()
}

class PolarBinGrid {
    -theta_edges: ndarray
    -r_edges: ndarray
    -n_t: int
    -n_r: int
    +map_coord_to_bin_idx(r: float, theta: float): (int, int) | None
    +bins_in_region(r0: float, t0: float, r1: float, t1: float): Iterable<(int,int)>
    +exposed_edges(bins: Set<(int,int)>): List<BinEdge>
    +bin_in_theta_arc(theta_start: float, theta_end: float, bin_start: float, bin_end: float): bool
}

%% =========================
%% View
%% =========================
class BinEdge {
    +theta_endpoints: (float, float)
    +r_endpoints: (float, float)
}

class PolarBinRenderer {
    -ax: Axes
    +plot_kwargs: dict
    +draw_edges(edges: List<BinEdge>, plot_kwargs: dict | None): List<Artist>
}

%% =========================
%% Controller
%% =========================
class SelectionOperation {
    <<enumeration>>
    REPLACE
    ADD
    SUBTRACT
}

class SelectorDragState {
    <<dataclass>>
    +drag_start: (float,float) | None
    +last_theta: float | None
    +last_preview_bins: Set<(int,int)> | None
    +mods: FrozenSet<string>
}

class SelectorDrawState {
    <<dataclass>>
    +selected_artists: List
    +hover_artists: List
}

class SiteSelector {
    +ax: Axes
    +grid: PolarBinGrid
    +renderer: PolarBinRenderer
    +model: BinSelectionModel
    +draw_tracker: SelectorDrawState
    +drag_tracker: SelectorDragState
    +operation: SelectionOperation
    +on_press(event)
    +on_motion(event)
    +on_release(event)
}

class SiteSelectorManager {
    -fig: Figure
    -_drag_owner: SiteSelector | None
    +register(selector, active=False)
    +set_active(selector)
}

%% =========================
%% Relationships
%% =========================

PolarBinGrid --> BinEdge : produces
PolarBinRenderer --> BinEdge : renders

SiteSelector --> PolarBinGrid : queries geometry
SiteSelector --> BinSelectionModel : updates selection
SiteSelector --> PolarBinRenderer : renders outlines
SiteSelector *-- SelectorDragState : owns
SiteSelector *-- SelectorDrawState : owns
SiteSelector --> SelectionOperation : mode

SiteSelectorManager o-- SiteSelector : routes events

%% =========================
%% Styling (MVC color coding)
%% =========================
classDef model fill:#E8F5E9,stroke:#2E7D32,stroke-width:1px,color:#1B5E20;
classDef view fill:#E3F2FD,stroke:#1565C0,stroke-width:1px,color:#0D47A1;
classDef controller fill:#FFF3E0,stroke:#EF6C00,stroke-width:1px,color:#E65100;
```

## Sequence Diagram Depicting User Actions and Code Response
```mermaid
sequenceDiagram
    autonumber
    actor User
    participant Canvas as Matplotlib Figure Canvas
    participant Mgr as SiteSelectorManager
    participant Sel as SiteSelector
    participant Model as BinSelectionModel

    note over User,Canvas: Mouse press inside Axes
    User->>Canvas: button_press_event
    Canvas->>Mgr: event(inaxes=Axes)
    Mgr->>Mgr: set drag owner\nlatch modifier keys
    Mgr->>Sel: on_press(event)
    Sel->>Sel: initialize drag state\ncompute initial preview
    Sel->>Canvas: draw_idle()

    note over User,Canvas: Mouse drag
    User->>Canvas: motion_notify_event
    Canvas->>Mgr: event
    Mgr->>Sel: on_motion(event)
    Sel->>Sel: update preview\n(unless cursor left Axes)
    Sel->>Canvas: draw_idle()

    note over User,Canvas: Mouse release (anywhere)
    User->>Canvas: button_release_event
    Canvas->>Mgr: event
    Mgr->>Sel: on_release(event)
    Sel->>Model: commit last preview selection
    Sel->>Sel: clear hover state\nreset drag state
    Sel->>Canvas: draw_idle()
    Mgr->>Mgr: clear drag owner

```
| If you want to change…                     | Edit this class                                |
| ------------------------------------------ | ---------------------------------------------- |
| Which bins are selected by a drag          | `PolarBinGrid`                                 |
| How angular wraparound works               | `PolarBinGrid.bin_in_theta_arc`                |
| How selection state is stored              | `BinSelectionModel`                            |
| How outlines are drawn                     | `PolarBinRenderer`                             |
| Selection semantics (replace/add/subtract) | `SiteSelector`                                 |
| Modifier-key behavior                      | `SiteSelectorManager._mods_from_mouse_event`   |
| Drag ownership / event routing             | `SiteSelectorManager`                          |
| Hover vs committed rendering               | `SiteSelector`                                 |
| Undo/redo behavior                         | Override `SiteSelector.on_selection_committed` |
