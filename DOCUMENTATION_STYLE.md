# Documentation style guide

This guide applies to documentation contributed by people and AI agents in
this repository. It covers Python docstrings, the top-level README, user
guides, configuration references, notebooks, and explanatory comments.

The goal is documentation that is accurate, useful at the point of need, and
easy to maintain. Documentation should explain the behavior a user or
maintainer needs to understand; it should not merely restate the implementation.

## Core principles

- **Be accurate.** Describe behavior that is present in the code or supported
  by a runnable example. Do not guess about defaults, units, supported inputs,
  or failure modes.
- **Be specific.** Name the relevant file, parameter, object, output, unit, or
  condition. Replace vague phrases such as “works as expected” with an
  observable result.
- **Be concise.** Include the smallest explanation that lets the reader use,
  debug, or extend the feature. Remove historical detail that does not affect
  current behavior.
- **Prefer user outcomes.** Explain what a setting changes, when to use it,
  and what the reader should expect—not just how the implementation is
  arranged.
- **Keep documentation close to its audience.** Put API behavior in
  docstrings, workflow instructions in user guides, and repository-wide rules
  in this file or the README.
- **Use the project’s vocabulary.** Match existing names such as `Rmax`,
  `atomsels`, `leaflet_sorting_algorithm`, Site, and SymmetricSite. Do not
  silently rename code concepts in prose.

## Include

Documentation should include information that changes what the reader does or
expects, including:

- the purpose and scope of a module, class, function, command, or workflow;
- required inputs, accepted types or formats, units, and meaningful defaults;
- the important relationship between arguments, especially when lists must
  have matching lengths or when one option changes the meaning of another;
- return values, generated files, mutations, plots, and other side effects;
- errors and validation rules a user can reasonably encounter, including how
  to correct them;
- assumptions, limitations, supported cases, and behavior that is intentionally
  not supported;
- non-obvious rationale, scientific assumptions, and invariants that explain
  why an approach was chosen or what must remain true when modifying it;
- short examples for non-obvious workflows, command syntax, configuration, or
  API usage;
- references to the source of a scientific method, external format, or
  domain-specific definition when that reference is necessary to interpret the
  result.

## Do not include

- claims that have not been checked against the implementation or an executable
  example;
- line-by-line descriptions of obvious code, private implementation details,
  or internal names that provide no help to a user or maintainer;
- duplicated documentation that will drift from a canonical explanation;
- obsolete behavior, commented-out alternatives, speculative future plans, or
  changelog entries in a reference section;
- unexplained acronyms, jargon, magic numbers, or units;
- promises about performance, compatibility, reproducibility, or scientific
  interpretation unless they are established and relevant;
- output that exposes secrets, local absolute paths, personal data, or
  machine-specific setup;
- instructions that depend on an unmentioned framework, package, or tool;
- filler such as “This function does what its name says” or “See above.”

## Python docstrings

Use a docstring for all public and private modules, classes, functions,
methods, and properties, including those with simple behavior.

Start with a one-sentence imperative or descriptive summary. Add detail only
when it helps the reader make a decision or use the API correctly. Follow the
NumPy-style sections already used in the Python package:

This excerpt documents the existing `dta.density.calculate_density` function;
the function body is omitted to focus on the docstring.

```python
def calculate_density(avg_counts: np.ndarray, grid: PolarBinGrid) -> np.ndarray:
    """Calculate the average bead density in each bin.

    Parameters
    ----------
    avg_counts : numpy.ndarray
        Two-dimensional array of average bead counts, indexed by radial bin
        and then angular bin. Its shape must match the grid's bin counts.
    grid : PolarBinGrid
        Grid supplying the bin counts and areas.

    Returns
    -------
    numpy.ndarray
        Average bead counts divided by bin areas, with the same shape as
        ``avg_counts``. Units are beads per unit area of the grid.

    Raises
    ------
    TypeError
        If ``avg_counts`` is not a NumPy array or ``grid`` is not a
        ``PolarBinGrid``.
    ValueError
        If ``avg_counts`` is not two-dimensional or its shape does not match
        ``(grid.r.n_bins, grid.theta.n_bins)``.
    """
```

Use these sections when applicable:

- `Parameters`: when this section is present, every parameter other than
  implicit `self` or `cls`, with its type and applicable units, default, and
  constraints;
- `Returns`: the returned value, its type, shape or structure, and units;
- `Raises`: exceptions callers can reasonably trigger and their causes;
- `Attributes`: public state exposed by a class;
- `Side Effects`: mutations, files, plots, logging, callbacks, or changes to
  external objects;
- `Notes`: constraints or relationships that do not fit naturally elsewhere;
- `Examples`: only examples that are short, correct, and useful.

Do not document an argument as optional unless the signature or implementation
actually accepts its omission. Accepting `None` is a separate property: state
whether `None` is accepted and what it means. A required argument can accept
`None`, and an argument with a default need not accept `None`.

Use double backticks for code identifiers in docstrings and single backticks
in Markdown prose. Use the exact type and unit terminology used by the code.
For a trivial method, a one-line docstring is preferable to a long template
with empty sections.

## User guides and README files

Organize task-oriented documentation around what the reader is trying to do:

1. State the purpose and prerequisites.
2. Give the shortest complete procedure, with commands and paths in code
   formatting.
3. Explain the settings or choices that require judgment.
4. Show the expected result and the most likely failure or validation check.
5. Link to deeper reference material rather than repeating it.

Use headings that describe the subject or task. Use numbered steps for a
sequence and bullets for independent options. Keep examples copyable: use
explicit placeholders such as `/absolute/path/to/config.tcl`, and say what the
reader must replace. Prefer repository-relative paths in general guidance and
avoid embedding a contributor’s local path.

Configuration references should follow the existing pattern in `README.md`:
describe each setting in terms of its effect, valid values or units, default
behavior, and interactions with related settings. When a rule applies only to
one algorithm or mode, say so directly beside that rule.

When documenting analysis and plotting workflows, distinguish the step that
computes data from the step that visualizes it.

The agent workflow for applying this guide lives in
[the documentation skill](.codex/skills/documentation-workflow/SKILL.md).
