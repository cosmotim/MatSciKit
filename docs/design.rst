Design Document
===============

For the full architecture design document, see
`DESIGN.md <https://github.com/cosmotim/MatSciKit/blob/main/docs/DESIGN.md>`_
on GitHub.

Key Design Principles
---------------------

1. **Separate I/O from analysis from plotting** — analysis functions never import matplotlib
2. **Analysis functions return data** — they compute and return, never plot
3. **Physical constants in one place** — ``constants.py``, not redefined per function
4. **Material parameters from crystallography** — CIF reader (future) or manual input
5. **No hardcoded paths** — all file paths are user-provided arguments
6. **Porosity is sample-specific** — not a property of the Material class
7. **Pipelines are independently runnable** — Pipeline 2 accepts θ_D/v_s as arguments
