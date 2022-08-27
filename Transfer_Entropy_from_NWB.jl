using Pkg
Pkg.activate("Transfer_entropy")
Pkg.status()

using Plots
using DelimitedFiles
using HDF5
using CausalityTools
using NaNStatistics

# plot backend
plotlyjs()         # or gr() plotlyjs() pyplot()

# open data
filename = "/Volumes/T7/LeMerre_dataset/mPFC/Aversion/273858_20200925-probe0.nwb"
nwb = h5open(filename, "r")

# Get units spike times, Load jagged arrays
unit_times_data = read(nwb["units/spike_times"]);
unit_times_idx = read(nwb["units/spike_times_index"]);
pushfirst!(unit_times_idx,1);
unit_ids = read(nwb["units/id"]);
spk_times = [unit_times_data[unit_times_idx[i]:unit_times_idx[i+1]] for i in 1:length(unit_ids)];

# PSTH function
function simple_psth(spk,eve,win,binSize)
    output = [histcountindices(spk, eve[i]+win[1]:binSize:eve[i]+win[2]) for i in 1:length(eve)];
    spk_cnt = output[1][1]
end

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
writedlm("/Users/pierre/Desktop/te_273858_20200925-probe0.csv", te)

# plot heatmap
heatmap(te./0.1, clims=(0, 0.02), c = :thermal)