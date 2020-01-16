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
using StatsBase
gr()

file = matopen("C:/Users/josep/OneDrive/Desktop/RayleighSpikePlotter/Data/CA87_06_38_02152006_unknown_mgp.mat") #change file here
spikeTimes = read(file, "spt") #spike times in s
time = read(file, "t") #time in seconds
dt=time[2]-time[1] #timestep, derived from time array
stimCurve = read(file,"u") #points along curve
choppedStimCurve= Float64[]
choppedTime=Float64[]
choppedSpikeTimes=zeros(length(spikeTimes)) #to put spikes which occur within window of interest

#plots stimulus
function plotStimulusWaveform()
    stimulusPlot = plot(time[:], stimCurve[:], legend=false,
    title="Stimulus Waveform")
    display(stimulusPlot)
end

#plots trimmed stimulus
function plotChoppedStimulusWaveform()
    choppedStimulusPlot = plot(choppedTime[:], choppedStimCurve[:],
    legend=false,title="Trimmed Stimulus Waveform")
    display(choppedStimulusPlot)
end

#Code identifying frequency
function findFrequency()
    possFreq = Float64[0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]
    #all of the frequencies which were tested
    locOfPeakRAW=findmax(stimCurve)
    cart = (locOfPeakRAW[2])
    startIndex = cart[2] #these 3 lines for extracting index of highest
    count = 0
    i = startIndex

    while stimCurve[i]>=stimCurve[i+1] || stimCurve[i-1]>=stimCurve[i+3]
        count +=1
        i +=1
    end

    cycleTime=count*dt*2
    #println(cycleTime)
    approxHz=1/cycleTime
    #println(approxHz)

    #and then subtract approxHz from each off possFreq
    #square answers. Whichever possFreq produces lowest result
    #must be intended stimulus frequency
end

#trims data to only include region of high stimulus amplitude
function chopData()
    positveThresh=percentile(stimCurve[:], 98)
    negativeThresh=percentile(stimCurve[:], 2)
    i=1
    while stimCurve[i]<=positveThresh && stimCurve[i]>=negativeThresh
        i+=1
    end
    j=length(stimCurve[:])
    while stimCurve[j]<=positveThresh && stimCurve[j]>=negativeThresh
        j-=1
    end
    global choppedStimCurve=stimCurve[i:j]
    global choppedTime=time[i:j]
    choppedStart=choppedTime[1]
    choppedEnd=choppedTime[length(choppedTime)]
    a=1
    for k in spikeTimes
        if k>=choppedStart && k<=choppedEnd
            global choppedSpikeTimes[a]=k
            a+=1
        end
    end
    global choppedSpikeTimes=choppedSpikeTimes[1:a-1]
    # println("start")
    # println(length(choppedStimCurve))
    # println(length(choppedTime))
    # println(choppedStart)
    # println(choppedEnd)
    # println(length(choppedSpikeTimes))
end

plotStimulusWaveform()
#findFrequency()
chopData()
plotChoppedStimulusWaveform()

#Code saying how many cycles we have in chopped bit

#Code using frequency to identify cycle length in s

#function transforming time in s to angle in radians based on frequency
#OR be more accurate and base transformation on stimCurve at
#each timepoint <- YES do this


#for every angle, while angle>2pi, angle=angle-2pi

# struct Spike_Count:
#     location :: Float64
#     frequency :: Int
# end
