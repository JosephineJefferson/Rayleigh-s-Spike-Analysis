#Julia script destined to read Matlab files including
#spike train data and a sine wave showing the stimulus.

#A key part of my LA internship. Will make Rayleigh Z plots
#aiming to ellucidate whether there is phase locking in
#vestibular afferents of animals treated with low doses of
#gentamycin.

#Basically making something to read matlab files about spikes
#and make a sensible Rayleigh Z plot.

#Author: Josephine Jefferson
#Start date: 14 Jan 2020

using MAT
using Plots
gr()


file = matopen("C:/Users/josep/OneDrive/Desktop/RayleighSpikePlotter/Data/CA87_06_38_02152006_unknown_mgp.mat") #change file here
spikeTimes = read(file, "spt")
time = read(file, "t") #time in seconds
dt=time[2]-time[1]
stimCurve = read(file,"u") #points along curve

stimulusPlot = plot(time[1:5463], stimCurve[1:5463], legend=false,
title="Stimulus Waveform")
display(stimulusPlot)

#Next up: identifying frequencies, and isolating cycles.
#Can I do calculus on this curve??
#Probably not... is just a set of points

#Code identifying frequency

possFreq = Float64[0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]

locOfPeakRAW=findmax(stimCurve)
cart = (locOfPeak[2])
startIndex = cart[2] #these 3 lines for extracting index of highest
count = 0
i = startIndex

while stimCurve[i]>=stimCurve[i+1] || stimCurve[i-1]>=stimCurve[i+3]
    global count +=1
    global i +=1
end

cycleTime=count*dt*2
println(cycleTime)
approxHz=1/cycleTime
println(approxHz)

#and then subtract approxHz from each off possFreq
#square answers. Whichever possFreq produces lowest result
#must be intended stimulus frequency


#Code chopping out middle bit of stimulus, certain number of
#high amplitude cycles

#Make same chop for all 3 relevant arrays

#Code saying how many cycles we have in chopped bit

#Code using frequency to identify cycle length in s

#function transforming time in s to angle in radians based on frequency
#OR be more accurate and base transformation on stimCurve at
#each timepoint

#for every angle, while angle>2pi, angle=angle-2pi

# struct Spike_Count:
#     location :: Float64
#     frequency :: Int
# end
