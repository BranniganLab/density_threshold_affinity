# Architecture and Design Rationale

This document explains **why** the interactive polar-bin selection system is designed the way it is.

It complements `README.md`:

- **README.md** answers: *“How do I use this?”*
- **ARCHITECTURE.md** answers: *“Why does this work the way it does?”*

If you are refactoring, extending, or debugging interaction behavior, **read this file first**.

---

## Problem Context

Interactive selection in Matplotlib is deceptively difficult due to:

- Inconsistent event semantics across backends (desktop vs Jupyter widgets)
- Modifier-key state that may not be reliably reported
- Mouse gestures that cross Axes boundaries
- The need to keep hover previews and final commits consistent

Early implementations tended to:
- Leak state across gestures
- Commit incorrect selections on mouse release
- Break when dragging outside an Axes
- Behave differently depending on backend

The current architecture exists to **eliminate those classes of bugs**.

---

## Architectural Goals

The system is designed to satisfy the following goals:

1. **Gesture correctness**
   - A selection gesture is atomic: press → drag → release
   - Final commit must equal the last hover preview

2. **Backend robustness**
   - Must work in Jupyter widget backends (e.g. `ipympl`)
   - Must not depend on reliable modifier-key reporting mid-gesture

3. **Multiple Axes safety**
   - Dragging in one Axes must not affect another
   - Switching Axes must not leak gesture state

4. **Separation of concerns**
   - Geometry, rendering, and interaction logic must not be interwoven
   - Each responsibility must have a single owner

---

## Model–View–Controller (MVC)

The system uses a strict MVC architecture.

### Model

**Responsibilities**
- Own domain state
- Define geometry rules
- Be deterministic and testable

**Classes**
- `PolarBinGrid`
- `BinSelectionModel`

**Explicit non-responsibilities**
- No Matplotlib event handling
- No rendering logic

---

### View

**Responsibilities**
- Draw visual representations on Matplotlib Axes
- Accept geometry and style inputs

**Classes**
- `PolarBinRenderer`
- `BinEdge`

**Explicit non-responsibilities**
- No selection semantics
- No gesture logic

---

### Controller

**Responsibilities**
- Translate user input into model updates
- Manage gesture lifecycle and transient state
- Enforce interaction invariants

**Classes**
- `SiteSelector`
- `SiteSelectorManager`
- `SelectorDragState`
- `SelectorDrawState`
- `SelectionOperation`

**Explicit non-responsibilities**
- No geometry rules
- No persistent domain state outside the model

---

## Key Architectural Decisions

### 1. Gesture ownership is latched at mouse press

**Decision**
- The selector receiving the mouse press becomes the *drag owner*
- All subsequent motion and release events are routed to that selector

**Why**
- Prevents cross-Axes interference
- Makes drag behavior predictable
- Avoids reliance on `event.inaxes` during motion

**Consequence**
- Hover preview freezes if the cursor leaves the Axes
- Release outside the Axes still commits correctly

---

### 2. Modifier keys are latched per gesture

**Decision**
- Modifier state (Shift / Control) is captured at mouse press
- Modifier state does not change mid-drag

**Why**
- Jupyter widget backends do not reliably report modifier changes
- Live modifier switching caused non-deterministic behavior

**Consequence**
- Users must release and re-press to change modes
- Behavior is consistent across environments

---

### 3. Preview selection is computed separately from commit

**Decision**
- `_apply_preview` computes a *final preview selection*
- Commit uses the preview result verbatim

**Why**
- Prevents recomputation mismatches on release
- Eliminates “delete everything” and inversion bugs
- Makes the invariant explicit:

  > Final commit equals last hover preview

**Consequence**
- Preview state must be stored explicitly
- Commit logic is simpler and safer

---

### 4. Manager never mutates selection state

**Decision**
- `SiteSelectorManager` only routes events
- All selection logic lives in `SiteSelector`

**Why**
- Prevents split-brain state
- Keeps selection semantics in one place
- Makes testing and reasoning tractable

---

## Invariants (Must Always Hold)

These invariants are relied upon throughout the codebase:

1. **One drag → one owner**
   - `_drag_owner` is set on press and cleared on release

2. **One gesture → one mode**
   - `SelectionOperation` does not change mid-gesture

3. **Preview == Commit**
   - The committed selection must exactly match the last hover preview

4. **No cross-Axes leakage**
   - A selector never reacts to motion from another Axes

5. **Model purity**
   - Models never depend on Matplotlib state

If any of these invariants are violated, interaction bugs will reappear.

---

## Known Tradeoffs

These behaviors are intentional and documented:

- **No live modifier switching**
  - Required for backend robustness

- **Hover freezes outside Axes**
  - Prevents accidental selection distortion

- **Release commits even outside Axes**
  - Matches user intent for “select up to edge”

These are not bugs; they are design choices.

---

## What Not to Do

The following changes will almost certainly reintroduce bugs:

- Computing selection deltas on release instead of committing previews
- Allowing modifier changes mid-drag
- Letting the manager modify selection state
- Mixing geometry logic into controllers
- Drawing directly from the model

If a proposed change violates an invariant, it should be rejected.

---

## Final Note

This architecture is the result of **iterative debugging of real interaction failures**.  
It may appear conservative, but each constraint exists to prevent a known class of bugs.

If you change something and think:

> “This restriction seems unnecessary…”

Check this file first — there is usually a reason.
