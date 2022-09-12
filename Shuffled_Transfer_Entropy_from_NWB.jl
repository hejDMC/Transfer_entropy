using Pkg
Pkg.activate("Transfer_entropy")
Pkg.status()

using Plots
using DelimitedFiles
using HDF5
using CausalityTools 
using NaNStatistics
using Glob
#using ProgressMeter
 
# plot backend
plotly()         # or gr() plotlyjs() pyplot()

# get NWB files
fileloc = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Passive/"
filelist = glob("*.nwb",fileloc)

#get PSTH function
include("simple_psth.jl")

# Saving folder
Destfolder = "/Volumes/labs/dmclab/Pierre/Transfert_Entropy/Passive/Shuffled/"

num_of_iterations = 30

f = filelist[3]
println("Working on "*f)

# for loop with multiple threads
#k = 1
Threads.@threads for k in 1:num_of_iterations
    
    # create output file name
    path_elements = split(f,['/','.'])
    csvname = "shuf_te_"*path_elements[9]*"_iteration"*string(k)*".csv"
    
    # if csv file already exists do nothing
    if isfile(Destfolder*csvname) == 0

        # open data
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
        te = zeros(length(unit_ids),length(unit_ids)) 
        for i = 1: length(unit_ids)
            for j = 1:length(unit_ids)         
                te[i,j] = transferentropy(circshift(p[i][1],rand(1:length(p))), p[j][1], est)
            end
        end 

        # Save csv
        writedlm(Destfolder*csvname, te)

    end
    
end

println("Done")
