using Pkg
#Pkg.add("Plots")
#Pkg.add("PlotlyJS") # or  "PyPlot","PlotlyJS"
#Pkg.add("HDF5")
using DelimitedFiles
using Plots
using HDF5
using Glob
using Statistics
using NaNStatistics

# plot backend
gr()         # or gr() plotlyjs() pyplot()

# Parameters
Task = "Aversion"
src1 = "/Volumes/labs/dmclab/Pierre/Transfert_Entropy/"*Task*"/"
filelist1 = glob("*.csv",src1)
src2 = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/"*Task*"/"
filelist2 = glob("*.nwb",src2)

# Spatial Binning
bins4histo = -5:0.1:0 #in mm

# preallocate
spa_TE = zeros(length(filelist1), length(bins4histo),length(bins4histo))

for f in 1:length(filelist1)
    #f = 2
    # Read csv
    te = readdlm(filelist1[f])
    # Read NWB
    nwb = h5open(filelist2[f], "r")

    # get ede labels
    ede = read(nwb["general/extracellular_ephys/electrodes/location"]);
    main_ch = read(nwb["units/electrodes"]).+ 1;

    # get DV coordinates
    dv = read(nwb["general/extracellular_ephys/electrodes/DV"]);
    dv_ch = dv[main_ch] 

    # sort matrix by DV
    new_ii = sortperm(dv_ch);
    te = te[new_ii, reverse(new_ii)] # sorted by dv
    xlab = reverse(ede[main_ch[new_ii]])
    ylab = ede[main_ch[new_ii]]

    heatmap(te, 
        clims=(0, 0.01),
        c = :thermal, 
        xticks=(1:10:size(te,1),xlab[1:10:size(te,1)]),
        xrotation=90,
        xlabel= "Input unit region", 
        yticks=(1:10:size(te,1),
        ylab[1:10:size(te,1)]),
        ylabel= "output unit region",
        size = (750, 700)
    )

        
    histc = histcountindices(dv_ch, bins4histo)
    bin_ids = 1:length(bins4histo) 
    spa_te = [mean(te[findall(x->x==b1, histc[2]),findall(x->x==b2, histc[2])]) for b1 in bin_ids, b2 in bin_ids]
    replace!(spa_te, NaN=>0)

    xl = ["-5.0" ,"-4.5", "-4.0", "-3.5","-3.0","-2.5","-2.0","-1.5","-1.0","-0.5","0"]
    yl = reverse(xl)

    spa_TE[f,:,:] = spa_te

end

T = mean(spa_TE, dims=1)
T = dropdims(T, dims = (findall(size(T) .== 1)...,))

    heatmap(T,
        clims=(0, 0.001),
        c = :thermal, 
        xticks=(0:10:length(bins4histo),yl),
        xrotation=90,
        xlabel= "Depth of input unit (mm)", 
        yticks=(0:10:length(bins4histo),xl),
        ylabel= "Depth of output unit (mm)",
        size = (750, 700)
    )

