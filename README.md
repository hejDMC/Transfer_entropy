# Transfert entropy

A set of functions to compute the transfert entropy between all pairs of single units recorded simultaneously. Take as an input a NWB file (HDF5). Save a matrix in csv and plot the obtained transfert entropy values as a heatmap. Builds on Julia's [CausalityTools.jl](https://github.com/JuliaDynamics/CausalityTools.jl) package.

You can obatain normalized transfert entropy values (te) by computing shuffled te values (circular shifts) for the same single units.

If you know the spatial location of your units. You can rearrange the units according to brain regions.