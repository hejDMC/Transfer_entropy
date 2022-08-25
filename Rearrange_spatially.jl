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
Mouse = "PL031"
Task = "Passive"

# Read csv
file2open = "/Volumes/labs/dmclab/Pierre/Transfert_Entropy/"*Task*"/te_"*Mouse*"_20190524-probe0.csv"
te = readdlm(file2open)

# Read shuffled folder
src = "/Volumes/labs/dmclab/Pierre/Transfert_Entropy/"*Task*"/Shuffled/"
filelist = glob("*PL031*.csv",src)
te_s = [readdlm(f) for f in filelist]
# get mean shuffled te matrix
m_te_s = mean(te_s)

# Normalized te matrix
n_te = te - m_te_s

# Read NWB
filename = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Passive/PL031_20190524-probe0.nwb"
nwb = h5open(filename, "r")

# get ede labels
ede = read(nwb["general/extracellular_ephys/electrodes/location"]);
main_ch = read(nwb["units/electrodes"]);

# get DV coordinates
dv = read(nwb["general/extracellular_ephys/electrodes/DV"]);
dv_ch = dv[main_ch] 

# sort matrix by DV
new_ii = sortperm(dv_ch);
n_te = n_te[new_ii, reverse(new_ii)] # sorted by dv
xlab = reverse(ede[main_ch[new_ii]])
ylab = ede[main_ch[new_ii]]

heatmap(n_te, 
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

# Spatial Binning
bins4histo = -5:0.1:0 #in mm    
histc = histcountindices(dv_ch, bins4histo)
bin_ids = 1:length(bins4histo) 
spa_te = [mean(n_te[findall(x->x==b1, histc[2]),findall(x->x==b2, histc[2])]) for b1 in bin_ids, b2 in bin_ids]
replace!(spa_te, NaN=>0)

xl = ["-5.0" ,"-4.5", "-4.0", "-3.5","-3.0","-2.5","-2.0","-1.5","-1.0","-0.5","0"]
yl = reverse(xl)

heatmap(spa_te, 
    clims=(0, 0.001),
    c = :thermal, 
    xticks=(0:5:51,yl),
    xrotation=90,
    xlabel= "Depth of input unit (mm)", 
    yticks=(0:5:51,xl),
    ylabel= "Depth of output unit (mm)",
    size = (750, 700)
    )

