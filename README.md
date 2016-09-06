# OE-Indices
Package to calculate the uncertainty of convection indices using Monte Carlo sampling.

These two scripts interact with the SHARPpy libraries to generate
convection indices with errorbars from OE-generated thermodynamic profiles.

The branch named "samples" will actually output the entire distribution of convection indices.

Currently these scripts do not support real-time use.  They will not append
to the files created.  This will be done before PECAN, preferibly refined during the
testing.

They also require the implementation of an IPython cluster to support parallelization.
Parallelization significantly speeds up the Monte Carlo convection index computations.
