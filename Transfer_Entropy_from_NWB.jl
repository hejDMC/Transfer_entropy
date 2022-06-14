
using Plots
using DelimitedFiles
using HDF5
using CausalityTools
include("simple_psth.jl")

# plot backend
plotlyjs()         # or gr() plotlyjs() pyplot()

function Transfert_entropy_for_NWB(filename)
# open data
#filename = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/152417_20191023-probe0.nwb"
nwb = h5open(filename, "r")

# Get units spike times, Load jagged arrays
unit_times_data = read(nwb["units/spike_times"]);
unit_times_idx = read(nwb["units/spike_times_index"]);
pushfirst!(unit_times_idx,1);
unit_ids = read(nwb["units/id"]);
spk_times = [unit_times_data[unit_times_idx[i]:unit_times_idx[i+1]] for i in 1:length(unit_ids)];

# list comprehension simple psth
rec_dur = maximum(maximum(spk_times))
p = [simple_psth(spk_times[i],0,[0 rec_dur],0.1) for i in 1:length(unit_ids)]

# Transfert entropy
# Binning-based estimator, discrete time aprroach
est = VisitationFrequency(RectangularBinning(4))
# list comprehension
#te = [transferentropy(p[i], p[j], est) for i in 1:length(unit_ids), j in 1:length(unit_ids)] #length(unit_ids)
# for loop with multiple threads
# preallocate]
te = zeros(length(unit_ids),length(unit_ids)) 
Threads.@threads for i = 1:length(unit_ids)
    for j = 1:length(unit_ids)
        te[i,j] = transferentropy(p[i], p[j], est)
    end
end 

# Save csv
writedlm("/Users/pielem/Desktop/te_152417_20191023-probe0.csv", te)

# plot heatmap
heatmap(te./0.1, clims=(0, 0.02), c = :thermal)

end



