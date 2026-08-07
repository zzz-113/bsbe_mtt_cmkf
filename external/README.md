# Third-party code

Files in this folder were written by others and are redistributed unmodified so
that the pipeline runs without additional downloads. They are **not** covered by
this repository's LICENSE; each retains the licence of its origin.

## `gsw_sound_speed.m`

Sound speed from the TEOS-10 75-term specific-volume expression.

- Source: [GSW Oceanographic Toolbox for MATLAB](https://www.teos-10.org/software.htm), version 3.06.12
- Authors: Paul Barker and Trevor McDougall
- Licence: MIT
- Reference: IOC, SCOR and IAPSO (2010), *The international thermodynamic
  equation of seawater – 2010*. Intergovernmental Oceanographic Commission,
  Manuals and Guides No. 56, UNESCO.

Used by `src/soundSpeedFromCTD.m`.

## `viridis.m`

The `viridis` colormap as a MATLAB lookup table.

- Source: colormap data from [matplotlib](https://matplotlib.org/), MATLAB port
  widely circulated on the MATLAB File Exchange
- Colormap authors: Nathaniel J. Smith, Stéfan van der Walt and Eric Firing
- Licence: the colormap data is released into the public domain (CC0) by its
  authors

Used by `src/plotTrackingOverview.m` for the echogram in panel (a).
