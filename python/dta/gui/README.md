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

class PolarBinGrid:::model
class BinEdge:::model
class BinAddress:::model

class SelectionRenderer:::view

class SiteSelector:::controller
class SiteSelectorManager:::controller
class SelectorDragState:::controller
class SelectorDrawState:::controller
class SelectorOperations:::controller
class BinSelection:::controller

<<Model>> PolarBinGrid
<<Model>> BinEdge
<<Model>> BinAddress

<<View>> SelectionRenderer

<<Controller>> SiteSelector
<<Controller>> SiteSelectorManager
<<Controller>> SelectorDragState
<<Controller>> SelectorDrawState
<<Controller>> SelectorOperations
<<Controller>> BinSelection


%% =========================
%% Model
%% =========================
namespace DTA.bin_logic.selection {
class BinSelection {
    -bins : Set~BinAddress~
    +set(bins: Set~BinAddress~)
    +add(bins: Set~BinAddress~)
    +remove(bins: Set~BinAddress~)
    +clear()
    +snapshot() FrozenSet~BinAddress~
    +get_bins() Set~BinAddress~
}
}

namespace DTA.bin_logic.polar_grid {
class PolarBinGrid {
    -r_edges: ndarray
    -theta_edges: ndarray
    -n_r: int
    -n_t: int
    +map_coord_to_bin_idx(r: float, theta: float) BinAddress | None
    +bins_in_region(r0: float, t0: float, r1: float, t1: float) Set~BinAddress~
    +exposed_edges(bins: Set~BinAddress) List~BinEdge~
    +bin_in_theta_arc(theta_start: float, theta_end: float, bin_start: float, bin_end: float) bool
    -_edge_geometry(bin: BinAddress, side: str) BinEdge
}
}

namespace DTA.bin_logic.utils {
class BinEdge {
    +r_endpoints: float, float
    +theta_endpoints: float, float
}

class BinAddress {
    +r_index: Int
    +theta_index: Int
}
}

%% =========================
%% View
%% =========================
namespace DTA.gui.renderers {
class SelectionRenderer {
    -ax: Axes
    +plot_kwargs : dict
    +draw_edges(edges: List~BinEdge~, plot_kwargs: dict | None) List~Artist~
    +shade_interior_region()*
}
}

%% =========================
%% Controller
%% =========================
namespace DTA.gui.selector_state {
class SelectorOperations {
    <<enumeration>>
    REPLACE
    ADD
    SUBTRACT
}

class SelectorDragState {
    <<dataclass>>
    +drag_start: float,float | None
    +last_theta: float | None
    +last_preview_bins: Set~BinAddress~ | None
    +mods: FrozenSet~string~
}

class SelectorDrawState {
    <<dataclass>>
    +selected_artists: List~Artist~
    +hover_artists: List~Artist~
}
}

namespace DTA.gui.selectors {
class SiteSelector {
    +ax: Axes
    +model: BinSelection
    +grid: PolarBinGrid
    +renderer: SelectionRenderer
    +drag_tracker: SelectorDragState
    +draw_tracker: SelectorDrawState
    +operation: SelectorOperations
    +on_activate()
    +on_deactivate()
    +on_press(event)
    +on_motion(event)
    +on_release(event)
    -bins_from_drag(start: tuple, end: tuple) Set~BinAddress~ | BinAddress | None
    -apply_preview(bins: Set~BinAddress~) Set~BinAddress~
    -apply_commit(bins: Set~BinAddress~)
    -commit_preview_selection(preview_bins: Set~BinAddress~)
    -draw_hover(bins: Set~BinAddress~)
    -draw_committed()
    -clear_artists(artists: List~Artist~)
    -on_selection_committed(before: Set~BinAddress~, after: Set~BinAddress~)*
}

class SiteSelectorManager {
    -fig: Figure
    -drag_owner: SiteSelector | None
    +register(selector: SiteSelector, active: Bool)
    +set_active(selector: SiteSelector)
    -mods_from_mouse_event(event: MouseEvent) Set~str~
    -dispatch(method: str) callable 
}
}
%% =========================
%% Relationships
%% =========================

SiteSelector o-- PolarBinGrid : defines polar lattice logic
SiteSelector o-- BinSelection : holds & updates bin selection
SiteSelector o-- SelectionRenderer : draws things
SiteSelector *-- SelectorDragState : keeps track of mouse drags and key-press modifiers
SiteSelector *-- SelectorDrawState : keeps track of what is currently drawn
SiteSelector ..> SelectorOperations : enumerates allowed modes

SiteSelectorManager o-- SiteSelector : coordinates behavior for one Axes

%% =========================
%% Styling (MVC color coding)
%% =========================
classDef model fill:#E8F5E9,stroke:#2E7D32,stroke-width:1px,color:#1B5E20;
classDef view fill:#E3F2FD,stroke:#1565C0,stroke-width:1px,color:#0D47A1;
classDef controller fill:#FFF3E0,stroke:#EF6C00,stroke-width:1px,color:#E65100;
```

## Sequence Diagram Depicting User Actions and Code Response
```mermaid
%% Self-messages require a trailing <br/> to avoid label/loop overlap (Mermaid limitation)
%%{init: {
  "sequence": {
    "messageMargin": 150,
    "rightAngles": true,
    "wrap": false,
    "messageAlign": "center"
  }
}}%%
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
    Mgr->>Mgr: set drag owner<br/>latch modifier keys<br/> 
    Mgr->>Sel: on_press(event)
    Sel->>Sel: initialize drag state<br/>compute initial preview<br/> 
    Sel->>Canvas: draw_idle()

    note over User,Canvas: Mouse drag
    User->>Canvas: motion_notify_event
    Canvas->>Mgr: event
    Mgr->>Sel: on_motion(event)
    Sel->>Sel: update preview<br/>(unless cursor left Axes)<br/> 
    Sel->>Canvas: draw_idle()

    note over User,Canvas: Mouse release (anywhere)
    User->>Canvas: button_release_event
    Canvas->>Mgr: event
    Mgr->>Sel: on_release(event)
    Sel->>Model: commit last preview selection
    Sel->>Sel: clear hover state<br/>reset drag state<br/> 
    Sel->>Canvas: draw_idle()
    Mgr->>Mgr: clear drag owner<br/> 

```
## Where to change what
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
