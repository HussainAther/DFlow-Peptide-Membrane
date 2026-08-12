# Repository scope

This repository contains source code, small configurations, tests, and project documentation needed to run and develop DFlow.

The cleaned distribution intentionally excludes generated simulation output, molecular-dynamics trajectories, figures, manuscripts, reference papers, operating-system metadata, caches, backup files, and embedded Git history. Those materials should be stored separately from the source repository.

`dflow_reversible.py` is the current reference runner. The modules under `src/` and the commands under `cli/` remain work in progress and are retained for continued consolidation.

The removed `src/space_init.py` was an incomplete, syntactically invalid fragment and was not imported elsewhere in the repository.
