using Pkg
Pkg.activate("Transfer_entropy")
Pkg.status()

using Plots
using DelimitedFiles
using HDF5
using CausalityTools 
using NaNStatistics
using Glob 
 
# plot backend
plotlyjs()         # or gr() plotlyjs() pyplot()

# get NWB files
fileloc = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/"
filelist = glob("*.nwb",fileloc)

#get PSTH function
include("simple_psth.jl")

# Saving folder
Destfolder = "/Volumes/labs/dmclab/Pierre/Transfert_Entropy/"

# Main loop through files
for f in filelist

    # create output file name
    path_elements = split(f,['/','.'])
    csvname = "te_"*path_elements[9]*".csv"
    # if csv file already exists do nothing
    if isfile(Destfolder*csvname) == 0

    # open data
    #filename = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/152417_20191023-probe0.nwb"
    nwb = h5open(f, "r")

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
    # for loop with multiple threads
    # preallocate]
    te = zeros(2,2) #length(unit_ids),length(unit_ids)
    Threads.@threads for i = 1: 2#length(unit_ids)
        for j = 1:2#length(unit_ids)
            te[i,j] = transferentropy(p[i], p[j], est)
        end
    end 

    # Save csv
    writedlm(Destfolder*csvname, te)

    # plot heatmap
    heatmap(te./0.1, clims=(0, 0.02), c = :thermal)

    end
end

