using DataFrames
using CSV
using LaTeXStrings
using HypothesisTests
import PyPlot
const plt = PyPlot
plt.pygui(true)
const path = "./results/resultsBinVar/No cones for separation/"
const pathPP = "./performanceProfile/"
const fields = [:Quality, 
                :Number_Of_Cycles, 
                :Max_Number_Of_Cycles, 
                :Total_number_of_nd_points, 
                :Total_number_of_nd_points_found_by_GM,
                :Time,
                :Nb_of_feasible_points_found,
                :Nb_of_maxTime_reached,
                :Nb_of_maxTrials_reached
                ]
const pathRatio = "./RatioBeforeAfter/"


const fieldsAndCorrespondingOrderFunction = Dict{Symbol,Function}( # used in performanceProfile
                                                                    :Quality => maximum,
                                                                    :Number_Of_Cycles => minimum, 
                                                                    :Max_Number_Of_Cycles => minimum, 
                                                                    :Total_number_of_nd_points => maximum, 
                                                                    :Total_number_of_nd_points_found_by_GM => maximum,
                                                                    :Time => minimum,
                                                                    :Nb_of_feasible_points_found => maximum,
                                                                    :Nb_of_maxTime_reached => minimum,
                                                                    :Nb_of_maxTrials_reached => minimum
                                                                 )
#=
If characteristicBis is null only the values for the given characteristic is given, 
the ratio of characteristic/characteristicBis is ploted otherwise.
=#
function plotQualities(filename::String, characteristic::Symbol, characteristicBis=nothing, refFile::String=path*"resultRef.csv")
    # parsing 
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refFile))
    # getting some informations
    names = df[!,:Instance]
    characteristicCurrent = (characteristicBis==nothing ? df[!,characteristic] : (df[!,characteristic] ./ df[!,characteristicBis]))
    characteristicRef = (characteristicBis==nothing ? refDf[!,characteristic] : (refDf[!,characteristic] ./ refDf[!,characteristicBis]))
    #println(characteristicCurrent)
    #println(characteristicRef)
    #println(names)
    
    xPos = collect(1:length(names))

    width = 0.4
    fig, ax1 = plt.subplots()
    #ax2 = ax1.twinx() # second axes
    fig.set_dpi(300)
    fig.set_size_inches(20.,13.)
    ax1.tick_params(labelsize=6)
    plt.bar(xPos .- width/2,characteristicCurrent,width=width/2,label="New method",color="blue",)
    plt.bar(xPos .+ width/2,characteristicRef,width=width/2,label="State-of-the-art results",color="green")

    # Setting the y_axis with a personalised density of ticks
    nbpoints = 20
    max_y = max(maximum(characteristicRef),maximum(characteristicCurrent))
    yticks = range(start=0.,step=max_y/nbpoints,stop=max_y)
    plt.yticks(yticks, labels=[string(v) for v in yticks])
    
    title::String = (characteristicBis==nothing ? string(characteristic) : string(characteristic)*"/"*string(characteristicBis))

    plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)))#,fontsize=5)
    plt.xlabel("Instance")
    plt.ylabel(title)
    plt.title("Comparison between GM ref and "*filename[1:end-3]*"-"*title)
    meanRef = sum(characteristicRef)/length(characteristicRef)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanRef,meanRef],color="red")
    meanCurrent = sum(characteristicCurrent)/length(characteristicCurrent)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color="orange")
    #ax2.set_yticks([meanCurrent,meanRef],["Means for SOTA","Mean new method"])
    #plt.yticks([meanRef],["Mean-New Method"],color="orange")
    #plt.yticks([meanCurrent],["Mean_State of the Art"],color="red")
    plt.legend()
    #plt.show()
    #println("Saving path : ", filename[1:end-5]*"/"*string(characteristic)*".png")
    plt.savefig(filename[1:end-4]*"/"*title*".png")
    plt.close(fig)
end

function plotQualitiesScatterLine(characteristic::Symbol, characteristicBis=nothing, refFile::String=path*"0.0-0.0.csv")
    # parsing 
    refDf = DataFrame(CSV.File(refFile))
    # getting some informations
    
    
    characteristicRef = (characteristicBis==nothing ? refDf[!,characteristic] : (refDf[!,characteristic] ./ refDf[!,characteristicBis]))
    #println(characteristicCurrent)
    #println(characteristicRef)
    #println(names)
    names = refDf[!,:Instance]
    xPos = collect(1:length(names))

    width = 0.4
    fig, ax1 = plt.subplots(figsize=(20,13),dpi=100)
    #ax2 = ax1.twinx() # second axes
    #fig.set_dpi(300)
    #fig.set_size_inches(20.,13.)
    ax1.tick_params(labelsize=6)
    
    plt.plot(xPos ,characteristicRef,color="green")
    plt.scatter(xPos ,characteristicRef,label="State-of-the-art results",color="green")

    # Setting the y_axis with a personalised density of ticks
    nbpoints = 20
    
    maxY = 0
    title::String = (characteristicBis==nothing ? string(characteristic) : string(characteristic)*"/"*string(characteristicBis))
    plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)))#,fontsize=5)
    plt.xlabel("Instance")
    plt.ylabel(title)
    plt.title("Comparison for different ratio of fixed binary variables "*title)
    meanRef = sum(characteristicRef)/length(characteristicRef)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanRef,meanRef],color="green",linestyle="dashed",label="Mean - ref State-of-the-art")

    colorIdx = 1
    colorSet = ["seagreen","maroon","chocolate","darkviolet","lawngreen","royalblue"]

    for filename in [file for file in readdir(path) if (file[end-3:end]==".csv") && (file!="0.0-0.0.csv")]
        df = DataFrame(CSV.File(path*filename))
        characteristicCurrent = (characteristicBis==nothing ? df[!,characteristic] : (df[!,characteristic] ./ df[!,characteristicBis]))
        maxy = max(maximum(characteristicRef),maximum(characteristicCurrent),maxY)
        plt.plot(xPos ,characteristicCurrent,color=colorSet[colorIdx])
        plt.scatter(xPos ,characteristicCurrent,label=filename[1:end-4],color=colorSet[colorIdx])
        meanCurrent = sum(characteristicCurrent)/length(characteristicCurrent)
        #stdDev = 100*(characteristicCurrent .- characteristicRef ./ characteristicRef) 
        pval = pvalue(HypothesisTests.SignedRankTest(characteristicCurrent,characteristicRef))
        plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color=colorSet[colorIdx],linestyle="dashed",label="Mean "*filename[1:end-4])
        println(string(characteristic)*"---"*filename[1:end-4]*" - "*string(pval))
        colorIdx += 1
    end 
    if maxY>0
        yticks = range(start=0.,step=maxY/nbpoints,stop=maxY)
        plt.yticks(yticks, labels=[string(v) for v in yticks])
    end
    #ax2.set_yticks([meanCurrent,meanRef],["Means for SOTA","Mean new method"])
    #plt.yticks([meanRef],["Mean-New Method"],color="orange")
    #plt.yticks([meanCurrent],["Mean_State of the Art"],color="red")
    plt.legend()
    #plt.show()
    #println("Saving path : ", filename[1:end-5]*"/"*string(characteristic)*".png")
    plt.savefig(path*title*"scatline.png")
    plt.close(fig)
end
#=
function plotQualitiesScatterLineParam(filename::String, characteristic::Symbol, characteristicBis=nothing, refFile::String=path*"resultRef.csv")
    # parsing 
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refFile))
    # getting some informations
    names = df[!,:Instance]
    characteristicCurrent = (characteristicBis==nothing ? df[!,characteristic] : (df[!,characteristic] ./ df[!,characteristicBis]))
    characteristicRef = (characteristicBis==nothing ? refDf[!,characteristic] : (refDf[!,characteristic] ./ refDf[!,characteristicBis]))
    #println(characteristicCurrent)
    #println(characteristicRef)
    #println(names)
    
    xPos = collect(1:length(names))

    width = 0.4
    fig, ax1 = plt.subplots()
    #ax2 = ax1.twinx() # second axes
    fig.set_dpi(300)
    fig.set_size_inches(20.,13.)
    ax1.tick_params(labelsize=6)
    plt.plot(xPos ,characteristicCurrent,color="blue",)
    plt.plot(xPos ,characteristicRef,color="green")
    plt.scatter(xPos ,characteristicCurrent,label="New method",color="blue")
    plt.scatter(xPos ,characteristicRef,label="State-of-the-art results",color="green")

    # Setting the y_axis with a personalised density of ticks
    nbpoints = 20
    max_y = max(maximum(characteristicRef),maximum(characteristicCurrent))
    yticks = range(start=0.,step=max_y/nbpoints,stop=max_y)
    plt.yticks(yticks, labels=[string(v) for v in yticks])
    
    title::String = (characteristicBis==nothing ? string(characteristic) : string(characteristic)*"/"*string(characteristicBis))

    plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)))#,fontsize=5)
    plt.xlabel("Instance")
    plt.ylabel(title)
    plt.title("Comparison between GM ref and "*filename[1:end-3]*"-"*title)
    meanRef = sum(characteristicRef)/length(characteristicRef)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanRef,meanRef],color="green",linestyle="dashed",label="Mean - ref State-of-the-art")
    meanCurrent = sum(characteristicCurrent)/length(characteristicCurrent)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color="blue",linestyle="dashed",label="Mean - new method")
    #ax2.set_yticks([meanCurrent,meanRef],["Means for SOTA","Mean new method"])
    #plt.yticks([meanRef],["Mean-New Method"],color="orange")
    #plt.yticks([meanCurrent],["Mean_State of the Art"],color="red")
    plt.legend()
    #plt.show()
    #println("Saving path : ", filename[1:end-5]*"/"*string(characteristic)*".png")
    plt.savefig(filename[1:end-4]*"/"*title*"scatline.png")
    plt.close(fig)
end=#

function plotPerformances()
    csvf::Vector{String} = [file for file in readdir(path) if file[end-3:end]==".csv"]

    println("Files :", csvf)
    for file in csvf
        try
            mkdir(path*file[1:end-4])
        catch; Base.IOError # the directory already exists
            println("WARNING: The directory already exists -> Figures inside the directory will still be updated") 
        end

        for characteristic in fields
                #plotQualitiesScatterLine(path*file, characteristic, characteristicBis=:Time) # plotQualities(path*file, characteristic, characteristicBis=:Time) before
            try
                plotQualitiesScatterLine(path*file, characteristic)
            catch ArgumentError # compatibility with the former version of GMBenchmark which includes less measures
                println("WARNING: retro-Compatibility mode activated")
                nothing
            end
        end
    end
end

function plotRatioMIP()
    directories = readdir(pathRatio)
    for dir in directories
        files = [file for file in readdir(pathRatio*dir) if (file[end-3:end]==".csv")] # to avoid existing .png files
        println(files)
        fig, ax = plt.subplots()
        for file in files
            nbgen = length(files)
            println(pathRatio*dir*"/"*file)
            data = DataFrame(CSV.File(pathRatio*dir*"/"*file))[!,:x1]
            #map!(x-> x==NaN ? 0. : x, data, data) # MODIFY THIS
            if !isempty(data)
                println("Data :", data)
                plt.plot(1:length(data),data,label=file[1:end-4])
                plt.xticks(1:length(data),labels=[string(i) for i=1:length(data)])
                plt.title("Instance "*dir[end-5:end-4]*"-Nb variables in R before / Nb of variables in R after")
            end
        end
        plt.legend()
        plt.show()
        plt.savefig(pathRatio*dir*"/"*dir[end-5:end-4]*"vg")
        plt.close(fig)
    end
end

#= From: https://link.springer.com/article/10.1007/s101070100263
Each subarray of data is a flavour of Gravity Machine. The lengths of names and data must be same.
The data must represents the same metric to be compared with each other.
=#
function performanceProfile(data::Vector{Vector{Float64}}, names::Vector{String}, title::String, operator::Function)
    
    ρs(τ::Float64, P::Vector{Float64}) = length([true for rps in P if rps <= τ])/length(P) # function as defined in the article
    
    @assert (length(data)==length(names)) "The number of names does not match the number of solvers"

    ns = length(data) # number of solvers
    nbInst = length(data[1]) # number of instances
    bestValues = [operator(dataForOnesolver) for dataForOnesolver in data]
    ratios::Vector{Vector{Float64}} = [collect(map(x-> (bestValues[s]==0. ? 0.0 : x/bestValues[s]),data[s])) for s in 1:ns]
    τrange = sort(unique(Iterators.flatten(ratios))) # all ratio can be captured
    rMax::Float64 = maximum([maximum(setOfRatios) for setOfRatios in ratios])
    println("rMax =", rMax)
    #τrange = range(start=0,stop=rMax,length=length(data[1])) # TODO: Can we be better?
    yaxis::Vector{Vector{Float64}} = [[ρs(τ, setOfRatios) for τ in τrange] for setOfRatios in ratios]

    fig, ax = plt.subplots()
    fig.set_dpi(300)
    fig.set_size_inches(20.,13.)
    plt.title(title)

    for s in eachindex(ratios)
        plt.plot(τrange,yaxis[s],label=names[s])
    end

    plt.legend()
    plt.show()
    plt.savefig(pathPP*title*".png")
    plt.close(fig)
end

function performanceProfileFromRaw()
    try 
        mkdir(pathPP)
    catch e;
        if e.code == -17 # The directory already exists
            println("WARNING: The figures directory"*pathPP*"already exists. Existing figures will be erased.")
        else
            println("Error not handled")
            rethrow()
        end
    end
    
    files::Vector{String} = [file for file in readdir("./results/";join=true) if file[end-3:end]==".csv"] # CSV files
    names::Vector{String} = [file[1:end-4] for file in readdir("./results/") if file[end-3:end]==".csv"]
    mixedData::Vector{DataFrame} = [DataFrame(CSV.File(file)) for file in files]
    #println(mixedData)

    for metric in fields
        try 
            data::Vector{Vector{Float64}} = [convert(Vector{Float64}, df[!,metric]) for df in mixedData]
            title = string(metric)
            performanceProfile(data,names,title, fieldsAndCorrespondingOrderFunction[metric])
        catch ArgumentError;
                println("[WARNING] The following metric has been ignored : "*string(metric))
        end
    end
end

function main()
    #plotPerformances()
    #plotRatioMIP()
    #performanceProfileFromRaw()
    for characteristic in fields
        plotQualitiesScatterLine(characteristic)
    end
end

main()