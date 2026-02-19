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
class BinAddress:::model
class Coordinate:::model
class BinEdge:::model

class SelectionRenderer:::view

class SiteSelector:::controller
class SiteSelectorManager:::controller
class SelectorDragState:::controller
class SelectorOperations:::controller
class BinSelection:::controller

<<Model>> PolarBinGrid
<<Model>> BinAddress
<<Model>> Coordinate
<<Model>> BinEdge

<<View>> SelectionRenderer

<<Controller>> SiteSelector
<<Controller>> SiteSelectorManager
<<Controller>> SelectorDragState
<<Controller>> SelectorOperations
<<Controller>> BinSelection


%% =========================
%% Model
%% =========================

namespace DTA.bin_logic.polar_grid {
class PolarBinGrid {
    -r_edges: ndarray
    -theta_edges: ndarray
    -n_r: int
    -n_t: int
    +map_coord_to_bin_idx(coord: Coordinate) BinAddress | None
    +bins_in_region(start: Coordinate, end: Coordinate) Set~BinAddress~
    +exposed_edges(bins: Set~BinAddress) List~BinEdge~
    +bin_in_theta_arc(theta_start: float, theta_end: float, bin_start: float, bin_end: float) bool
    -_edge_geometry(bin: BinAddress, side: str) BinEdge
}
}

namespace DTA.bin_logic.utils {
class NamedTuple {}
class BinAddress {
    +r_index: Int
    +theta_index: Int
}
class Coordinate {
    +r_coord: float
    +theta_coord: float
}
class BinEdge {
    +endpoint1: Coordinate
    +endpoint2: Coordinate
}
}

%% =========================
%% View
%% =========================
namespace DTA.gui.renderers {
class SelectionRenderer {
    -ax: Axes
    +plot_kwargs : dict
    +selected_artists: List~Artist~
    +hover_artists: List~Artist~
    +draw_edges(edges: List~BinEdge~, plot_kwargs: dict | None) List~Artist~
    +shade_interior_region()*
}
}

%% =========================
%% Controller
%% =========================

namespace DTA.gui.selector_state {
class BinSelection {
    -bins : Set~BinAddress~
    -selection_history*
    +update(bins: Set~BinAddress~)
    +get_bins() Set~BinAddress~
    +clear()
    +undo()*
    +redo()*
}

class SelectorOperations {
    <<enumeration>>
    REPLACE
    ADD
    SUBTRACT
}

class SelectorDragState {
    <<dataclass>>
    +drag_start: Coordinate | None
    +last_theta: float | None
    +current_preview_bins: Set~BinAddress~ | None
    +mods: FrozenSet~string~
}
}

namespace DTA.gui.selectors {
class SiteSelector {
    +selection: BinSelection
    +grid: PolarBinGrid
    +renderer: SelectionRenderer
    +drag_tracker: SelectorDragState
    +operation: SelectorOperations
    +on_activate()
    +on_press(event)
    +on_motion(event)
    +on_release(event)
    +on_deactivate()
    -calculate_preview_bins(bins: Set~BinAddress~) Set~BinAddress~
    -draw_preview(bins: Set~BinAddress~)
    -draw_selection()
    -clear_artists(artists: List~Artist~)
    -save_to_selection_history(before: Set~BinAddress~, after: Set~BinAddress~)*
}

class SiteSelectorManager {
    +fig: Figure
    -selectors: dict
    -active: dict
    -drag_owner: SiteSelector | None
    -cids: List~int~
    +register(selector: SiteSelector, active: bool)
    +set_active(selector: SiteSelector)
    -on_press_event(event: MouseEvent)
    -on_drag_event(event: MouseEvent)
    -on_release_event(event: MouseEvent)
    -mods_from_mouse_event(event: MouseEvent) Set~str~
}
}
%% =========================
%% Relationships
%% =========================

SiteSelector *-- PolarBinGrid : defines polar lattice logic
SiteSelector *-- SelectionRenderer : draws things & keeps track of what's drawn
SiteSelector *-- BinSelection : holds & updates bin selection
SiteSelector *-- SelectorDragState : keeps track of mouse drags and key-press modifiers
SiteSelector ..> SelectorOperations : enumerates allowed modes

SiteSelectorManager o-- SiteSelector : coordinates behavior for one Axes

NamedTuple <|-- BinAddress : is a
NamedTuple <|-- Coordinate : is a
NamedTuple <|-- BinEdge : is a

%% =========================
%% Styling (MVC color coding)
%% =========================
classDef model fill:#E8F5E9,stroke:#2E7D32,stroke-width:1px,color:#1B5E20;
classDef view fill:#E3F2FD,stroke:#1565C0,stroke-width:1px,color:#0D47A1;
classDef controller fill:#FFF3E0,stroke:#EF6C00,stroke-width:1px,color:#E65100;
```

```mermaid
classDiagram
  direction TB
  namespace Legend {
    class Owner
    class Dog
    class SocialSecurityNumber
    class Handedness {
      RIGHTHANDED
      LEFTHANDED
      AMBIDEXTROUS
      NO_HANDS
    }}
    Owner o-- Dog : Owner "has a" Dog - they could have multiple dogs
    Owner *-- SocialSecurityNumber : Owner "has a" SSN - can only have one
    Owner ..> Handedness : Owner "can be one of these enumerations"
  
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
    participant Model as BinSelection

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
```mermaid
flowchart TD
    A[Event received by SiteSelectorManager] --> B{Event type?}

    B -->|button_press_event| C[SiteSelectorManager._on_press_event]
    B -->|motion_notify_event| D[SiteSelectorManager._on_motion_event]
    B -->|button_release_event| E[SiteSelectorManager._on_release_event]

    %% PRESS FLOW
    C --> C1{Is MouseEvent inside an Axes?}
    C1 -->|No| C0
    C1 -->|Yes| C2{Is there an active SiteSelector assigned to that Axes?}
    C2 --> |No| C0
    C2 --> |Yes| C3[self._drag_owner = active SiteSelector]
    C3 --> C4[Check if shift/ctrl pressed; save to self._drag_owner.drag_tracker.mods ]
    C4 --> C5[Call self._drag_owner's on_press method] --> C6{Was anything new drawn?}
    C6 --> |No| C0
    C6 --> |Yes| C7[self.fig.canvas.draw_idle]
    C7 --> C0[Return]

    %% MOTION FLOW
    D --> D1{Is there a self._drag_owner?}
    D1 -->|No| D0
    D1 -->|Yes| D2[Call self._drag_owner's on_motion method]
    D2 --> D6{Was anything new drawn?}
    D6 --> |No| D0
    D6 --> |Yes| D7[self.fig.canvas.draw_idle]
    D7 --> D0[Return]

    %% RELEASE FLOW
    E --> E1{Is there a self._drag_owner?}
    E1 --> |No|E0[Return]
    E1 --> |Yes|E2[Call self._drag_owner's on_release method]
    E2 --> E3{Was anything new drawn?}
    E3 --> |No|E6
    E3 --> |Yes|E4[self.fig.canvas.draw_idle]
    E4 -->E6[Reset self._drag_owner]
    E6 --> E0
```
## Where to change what
| If you want to change…                     | Edit this class                                |
| ------------------------------------------ | ---------------------------------------------- |
| Which bins are selected by a drag          | `PolarBinGrid`                                 |
| How angular wraparound works               | `PolarBinGrid.bin_in_theta_arc`                |
| How selection state is stored              | `BinSelection`                                 |
| How outlines are drawn                     | `PolarBinRenderer`                             |
| Selection semantics (replace/add/subtract) | `SiteSelector`                                 |
| Modifier-key behavior                      | `SiteSelectorManager._mods_from_mouse_event`   |
| Drag ownership / event routing             | `SiteSelectorManager`                          |
| Hover vs committed rendering               | `SiteSelector`                                 |
| Undo/redo behavior                         | Override `SiteSelector.on_selection_committed` |
