---
name: documentation-workflow
description: Create or update software documentation, including docstrings, READMEs, user guides, and configuration references, while verifying claims against the repository.
metadata:
  short-description: Create and review software documentation
---

# Documentation workflow

Use this skill when creating or editing documentation that explains software
behavior to users, maintainers, or future agents. It applies to Python
docstrings, Markdown user guides, READMEs, configuration references, and
documentation comments. It does not replace specialized artifact skills for
Word documents, PDFs, presentations, or spreadsheets.

## Workflow

1. Identify the audience, task, and documentation location. Search the
   repository for equivalent functionality, terminology, examples, and
   existing guidance before writing.
2. Before editing, read the repository's
   [documentation style guide](../../../DOCUMENTATION_STYLE.md) and any
   applicable contributor guidance. Use the guide for content, terminology,
   formatting, and docstring conventions.
3. Read the relevant implementation, tests, configuration examples, and nearby
   documentation. Verify behavior, defaults, types, units, outputs, side
   effects, limitations, and failure modes from evidence in the repository.
   If implementation, tests, and documented intent disagree, report the
   discrepancy. Do not silently describe a suspected bug as the intended
   contract or change runtime behavior during a documentation-only task.
4. Apply the style guide to the requested documentation. Preserve surrounding
   structure unless the task calls for reorganizing it. When behavior changes
   are part of the task, update the affected documentation to match.
5. Keep the change scoped. Preserve existing behavior and surrounding
   structure, and do not invent APIs, defaults, citations, supported cases, or
   scientific claims. If a necessary fact cannot be verified, flag the gap
   instead of guessing.
6. After editing, reread the documentation style guide and review the changed
   documentation against its applicable guidelines at least once. Check both
   what should be included and what should be omitted, along with terminology,
   formatting, and docstring conventions. Correct any mismatches and recheck
   the corrections before proceeding.
7. Select validation based on the change using the guidance below, then inspect
   the final diff for unrelated edits.
8. Report the documentation changes, validation performed, and unresolved
   issues or uncertainty.

## Validation

- For prose edits, inspect Markdown structure, code fences, changed links, and
  whitespace. A wording change alone does not require running the full test
  suite.
- For API documentation, compare signatures, defaults, types, units, array
  shapes, side effects, and exception conditions with the relevant code. Use
  targeted tests when behavior is uncertain.
- For commands and examples, execute them when practical within the authorized
  scope. If execution needs unavailable data, dependencies, or external
  services, inspect the example and report that it was not executed.
- Run relevant existing documentation checks when available. Do not introduce
  a new tool or dependency solely to validate a documentation edit without
  asking.

In the final report, distinguish checks that were run from examples or claims
that were only inspected, and identify unresolved discrepancies.

## Completion standard

A documentation change is complete when a reader can follow it without relying
on undocumented assumptions, its claims are supported by current repository
evidence or authoritative references, the review against the style guide is
complete, and the final report states what was checked and what remains
unresolved.
