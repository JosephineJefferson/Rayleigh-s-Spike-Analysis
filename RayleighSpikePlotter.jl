#Julia script destined to read Matlab files including
#spike train data and a sine wave showing the stimulus.

#A part of my LA internship. Will make Rayleigh Z plots
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
using DSP
#using LsqFit
#using ApproxFun
gr()

file = matopen("C:/Users/josep/OneDrive/Desktop/RayleighSpikePlotter/Data/CA87_15_89_02152006_unknown_mgp.mat") #change file here
spikeTimes = read(file, "spt") #spike times in s
time = read(file, "t") #time in seconds
dt=time[2]-time[1] #timestep, derived from time array
stimCurve = read(file,"u") #points along curve
chopping=true
choppedStimCurve= Float64[]
choppedTime=Float64[]
choppedSpikeTimes=zeros(length(spikeTimes)) #to put spikes which occur within window of interest
freq=0
#Parameters of equation
A=0
k=0
b=0
#array to store y values of fitted line
simData=zeros(length(time))
#array to store converted spike times
sptAsRadians=zeros(length(spikeTimes))

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

function plotStimRelativeToAngle(chopped=chopping)
    if chopped == false
        t = time[:]
    else
        t = choppedTime[:]
    end
    radiansg=zeros(length(t))
    co=1
    while co<=length(radiansg)
        radiansg[co] = timeToRadians(t[co])
        co +=1
    end
    display(plot(radiansg[:],simData[:]))
end

#trims data to only include region of high stimulus amplitude
function chopData()
    positveThresh=percentile(stimCurve[:], 99)
    negativeThresh=percentile(stimCurve[:], 1)
    i=1
    while stimCurve[i]<positveThresh && stimCurve[i]>negativeThresh
        i+=1
    end
    j=length(stimCurve)
    while stimCurve[j]<positveThresh && stimCurve[j]>negativeThresh
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
    global simData = zeros(length(choppedTime))
    global sptAsRadians = zeros(length(choppedSpikeTimes))
end

#finds frequency of high amplitude portion of data
function mikesPeriodogram(show=false, chopped=chopping)
    if chopped == false
        s = stimCurve[:]
    else
        s = choppedStimCurve[:]
    end
    pg = periodogram(s[:], fs = 1/dt)         # signal power as a fcn of freq
    i15Hz = minimum(findall(x -> x>15.0, pg.freq)) # index of 15Hz
    f = pg.freq[findmax(pg.power)[2]]     # frequency of peak power
    if show==true
        display(plot(pg.freq[1:i15Hz], pg.power[1:i15Hz],
            title="power", xlabel = "frequency (Hz)")) # plot up to 15Hz
        display(plot(pg.freq[1:i15Hz], pg.power[1:i15Hz], yaxis=:log,
            title="log power", xlabel="frequency (Hz) (nb power at 0 = signal mean )"))
        println("Frequency = ", f, " (", round(f; digits=1), ")")
    end
    return f
end

#finds the parameters of the fitted sine wave
function findParams(chopped=chopping)
    if chopped == false
        t = time[:]
        s = stimCurve[:]
    else
        t = choppedTime[:]
        s = choppedStimCurve[:]
    end
    global A = findmax(s[:])[1]
    global k = freq*2*pi
    index=length(t)÷2
    while s[index]>=A || s[index]<=-A
        index+=20
    end
    global b = asin(s[index]/A)/-k + t[index]
end

#called by other functions to make calculations based on
#fitted curve
function yBasedOnTime(t1::Float64)
    answer = A*sin(k*(t1-b))
    return answer
end

#plots the fit curve over the stimulus wave
function plotFitOnStim(chopped=chopping)
    if chopped == false
        t = time[:]
        s = stimCurve[:]
    else
        t = choppedTime[:]
        s = choppedStimCurve[:]
    end
    c=1
    while c<=length(simData)
        global simData[c]= yBasedOnTime(t[c])
        c += 1
    end
    p1=plot(t[:], s[:], label="stimulus")
    plot!(t[:], simData[:], label="fit")
    display(p1)
end

#finds y intercepts of fitted curve and returns an array with
#the indexes of these intercepts
#called by timeToRadians
function findFittedIntercepts(chopped=chopping)
    if chopped == false
        t = time[:]
    else
        t = choppedTime[:]
    end
    u=simData
    index=1
    interceptNo=0
    interceptsInd=zeros(length(t))
    while index<length(u)
        if (u[index]>=0 && u[index+1]<=0) || (u[index]<=0 && u[index+1]>=0)
            interceptNo+=1
            interceptsInd[interceptNo]=index
            index+=20
        end
        index+=1
    end
    interceptsInd=interceptsInd[1:interceptNo]
    # interceptvals=zeros(length(interceptsInd))
    # c=1
    # while c>=length(interceptsInd)
    #     interceptvals[c]=u[interceptsInd[c]]
    #     c+=1
    # end
    # println(interceptvals)
    return interceptsInd
end

#takes a time and converts it to radians based on phase length
#and a shift. Shift calculated based on time and gradient-sign
#of first y intercept
function timeToRadians(t1::Float64, chopped=chopping)
    if chopped == false
        t = time[:]
    else
        t = choppedTime[:]
    end
    u=simData[:]
    phLen = 1/freq
    intcpt1 = Int(findFittedIntercepts()[1])
    mult=2pi/phLen
    shift=0
    if u[intcpt1]<u[intcpt1+1]
        shift=-mult*t[intcpt1]
    else
        shift=pi-mult*t[intcpt1]
    end
    answer = t1*mult + shift
    while answer>=2pi
        answer=answer-2pi
    end
    while answer<0
        answer+=2pi
    end
    return answer
end

#converts all (relevant) spiketimes to angles in radians
function convertSpikeTimes(chopped=chopping)
    if chopped == false
        spt = spikeTimes[:]
    else
        spt = choppedSpikeTimes[:]
    end
    count=1
    while count<=length(spt)
        sptAsRadians[count]=timeToRadians(spt[count])
        count+=1
    end
end

#now we do stats!

#main code, calling functions

if chopping
    chopData()
end
freq = mikesPeriodogram()
findParams()
plotFitOnStim()
convertSpikeTimes()
