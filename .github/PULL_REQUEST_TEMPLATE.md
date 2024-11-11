## Description
Provide a brief description of the PR's purpose here.

## Usage Changes
Document any changes to the parameters

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Pre-Review checklist (PR maker)
- [ ] Approval tests pass
- [ ] Acceptance tests pass
    - [ ] All three modules run without errors 
    - [ ] ModCoor and CoorSpec agree (no subsampling). Based on plots.
    - [ ] New acceptance tests are documented 
    - [ ] Running `do_run.sh params.conf` passes the acceptance tests
    - [ ] ModCoor r^2 values converge to 1 with sufficient sampling (1000 samples)
- [ ] Linting isn't worse than it was (linting checks pass)
- [ ] Pytests pass
- [ ] New functions are documented
- [ ] New configuration parameters are documented (Usage changes above)

## Review checklist (Reviewer)
- [ ] Approval tests pass
- [ ] Acceptance tests pass
    - [ ] All three modules run without errors 
    - [ ] ModCoor and CoorSpec agree (no subsampling). Based on plots.
    - [ ] New acceptance tests are documented 
    - [ ] Running `do_run.sh params.conf` passes the acceptance tests
    - [ ] ModCoor r^2 values converge to 1 with sufficient sampling (1000 samples)
- [ ] Linting scores haven't gotten worse (linting checks pass)
- [ ] Pytests pass
- [ ] New functions are documented 
- [ ] New functions/variables are named appropriately
- [ ] No missed "low-hanging fruit" that would substantially aid readability.
- [ ] Any "high-hanging" or "rotten" fruit is added to the issues list.
- [ ] I understand what the changes are doing and how
- [ ] I understand the motivation for this PR (the PR itself is appropriately documented)

## Status
- [ ] Ready for review