using Pkg
Pkg.activate("Transfer_entropy")
Pkg.status()


using Glob
include("Transfer_Entropy_from_NWB.jl")

# plot backend
plotlyjs()         # or gr() plotlyjs() pyplot()

fileloc = "/Volumes/labs/dmclab/Pierre/NPX_Database/mPFC/Aversion/"
filelist = glob("*.nwb",fileloc)

for f in filelist
    Transfert_entropy_for_NWB(f)
end