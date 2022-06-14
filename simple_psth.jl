using NaNStatistics

# PSTH function
function simple_psth(spk,eve,win,binSize)
    output = [histcountindices(spk, eve[i]+win[1]:binSize:eve[i]+win[2]) for i in 1:length(eve)];
    spk_cnt = output[1][1]
end