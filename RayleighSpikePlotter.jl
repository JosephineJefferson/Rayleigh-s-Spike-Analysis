#Julia script reading Matlab files including
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
gr()

file = matopen("C:/Users/josep/OneDrive/Desktop/RayleighSpikePlotter/Data/CA88_08_33_02222006_3.2_mgp.mat") #change file here
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
simStim=zeros(length(time))
#array to store converted spike times
sptAsRadians=zeros(length(spikeTimes))
#values related to mean vector
meanVsize=0
meanVangle=0

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
    display(plot(radiansg[:],simStim[:]))
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
    global simStim = zeros(length(choppedTime))
    global sptAsRadians = zeros(length(choppedSpikeTimes))
end

#finds frequency of high amplitude portion of data
function mikesPeriodogram( show=false, correct=true, chopped=chopping)
    if chopped == false
        s = stimCurve[:]
    else
        s = choppedStimCurve[:]
    end
    pg = periodogram(s[:], fs = 1/dt)         # signal power as a fcn of freq
    i15Hz = minimum(findall(x -> x>15.0, pg.freq)) # index of 15Hz
    f = pg.freq[findmax(pg.power)[2]]     # frequency of peak power
    if show
        display(plot(pg.freq[1:i15Hz], pg.power[1:i15Hz],
            title="power", xlabel = "frequency (Hz)")) # plot up to 15Hz
        display(plot(pg.freq[1:i15Hz], pg.power[1:i15Hz], yaxis=:log,
            title="log power", xlabel="frequency (Hz) (nb power at 0 = signal mean )"))
        println("Frequency = ", f, " (", round(f; digits=1), ")")
    end
    freq=f
    if correct
        poss = [0.2,0.4,0.8,1.6,3.2,6.4,12.8]
        smallDiff=20.0
        count=1
        while abs(f-poss[count])<smallDiff
            smallDiff=abs(f-poss[count])
            freq=poss[count]
            count+=1
        end
        if show
            println("Corrected Frequency: ", freq)
        end
    end
    return freq
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
    #b calculation is unreliable
    global b = asin(s[index]/A)/-k + t[index]
    #global b = acsc(A/s[index])/-k + t[index]
end

#called by other functions to make calculations based on
#fitted curve
function yBasedOnTime(t1::Float64)
    answer = A*sin(k*(t1-b))
    return answer
end

#plots the fit curve over the stimulus wave
function mapFitOnStim(show=false, chopped=chopping)
    if chopped == false
        t = time[:]
        s = stimCurve[:]
    else
        t = choppedTime[:]
        s = choppedStimCurve[:]
    end
    c=1
    while c<=length(simStim)
        global simStim[c]= yBasedOnTime(t[c])
        c += 1
    end
    if show
        p1=plot(t[:], s[:], label="stimulus")
        plot!(t[:], simStim[:], label="fit")
        display(p1)
    end
end

#finds x intercepts of fitted curve and returns an array with
#the indexes of these intercepts
#called by timeToRadians
function findFittedIntercepts(chopped=chopping)
    if chopped == false
        t = time[:]
    else
        t = choppedTime[:]
    end
    u=simStim
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
    u=simStim[:]
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
    while answer>=2*pi
        answer=answer-2*pi
    end
    while answer<0
        answer+=2*pi
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

#now we do stats! badly, apparently
function meanVectorCalc(print=false)
    n = length(sptAsRadians) #-1
    sumCos=0.0
    for z in sptAsRadians[:]#[1:length(sptAsRadians)-1]
        sumCos+=cos(z)
    end
    X=sumCos/n

    sumSin=0.0
    for z in sptAsRadians[:]#[1:length(sptAsRadians)-1]
        sumSin+=sin(z)
    end
    Y=sumSin/n

    r=sqrt((X^2) + (Y^2))
    global meanVsize=r
    angle=atan(Y,X)
    global meanVangle = angle
    if print
        println("Vector strength: ", r)
        if angle<0
            angle+=2pi
        end
        println("Vector angle (rad): ", angle)
        println("Vector angle (deg): ", rad2deg(angle))
    end
end

#function displaying polarplot incl mean vector
function displayPolar()
    scale=1
    if meanVsize<0.01
        scale=0.02
    elseif meanVsize<0.05
        scale=0.1
    elseif meanVsize<0.1
        scale=0.2
    end
    x=zeros(length(sptAsRadians))
    count=1
    while count<=length(x)
        x[count]=scale+scale/10
        count+=1
    end

    sizeArray = 0:meanVsize:meanVsize#*20
    angleArray = zeros(length(sizeArray))
    count=1
    while count <= length(angleArray)
        angleArray[count]=meanVangle
        count+=1
    end

    polarplot=scatter(sptAsRadians[:], x, proj=:polar, label="Spikes")
    polarplot=plot!(angleArray[:], sizeArray[:], proj=:polar,
        yaxis=[0,scale+scale/5], linewidth=3, label="Mean Vector")
    display(polarplot)
end

#main code, calling functions

if chopping
    chopData()
end
freq = mikesPeriodogram(true)
findParams()
mapFitOnStim(true)
#convertSpikeTimes()
#plotStimRelativeToAngle()
#meanVectorCalc(true)
#displayPolar()
