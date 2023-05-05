using DataFrames
using CSV
using LaTeXStrings
using HypothesisTests
import PyPlot
const plt = PyPlot
plt.pygui(true)
const path = "./results/resultsBinVar/21Aprils/"
const refPath = "./results/Ref/"
const pathPP = "./performanceProfile/"
const pathrefExact = "./results/timeExactSolving/purifiedData21April.csv" # "Instances","Times","Quality","Number of (weakly) efficient points","Number of variables","Number of constraints"
const refDfExact = DataFrame(CSV.File(pathrefExact))
global orderByVar = sortperm(refDfExact[!,:Instances],by=x->refDfExact[refDfExact.Instances .== x,"Number of variables"][1])
global orderByCons = sortperm(refDfExact[!,:Instances],by=x->refDfExact[refDfExact.Instances .== x,"Number of constraints"][1])
const fields = [:Quality, 
                :Number_Of_Cycles, 
                :Max_Number_Of_Cycles, 
                :Total_number_of_nd_points, 
                :Total_number_of_nd_points_found_by_GM,
                :Time,
                :Nb_of_feasible_points_found,
                :Nb_of_maxTime_reached,
                :Nb_of_maxTrials_reached,
                :Avg_Ratio_Non_Bin_Var,
                :Supported_Sol
                ]
const pathRatio = "./RatioBeforeAfter/"

const subDir = ["MILPProjection/","NonMILPProjection/"]

const fieldsAndCorrespondingOrderFunction = Dict{Symbol,Function}( # used in performanceProfile
                                                                    :Quality => maximum,
                                                                    :Number_Of_Cycles => minimum, 
                                                                    :Max_Number_Of_Cycles => minimum, 
                                                                    :Total_number_of_nd_points => maximum, 
                                                                    :Total_number_of_nd_points_found_by_GM => maximum,
                                                                    :Time => minimum,
                                                                    :Nb_of_feasible_points_found => maximum,
                                                                    :Nb_of_maxTime_reached => minimum,
                                                                    :Nb_of_maxTrials_reached => minimum,
                                                                    :Avg_Ratio_Non_Bin_Var => minimum
                                                                 )

const characteristicToName = Dict{Symbol,String}(
                :Quality => "Quality (%)", 
                :Number_Of_Cycles => "Number of cycles", 
                :Max_Number_Of_Cycles => "Max number of cycles", 
                :Total_number_of_nd_points => "Exact number of non-dominated points", 
                :Total_number_of_nd_points_found_by_GM => "Number of non-dominated points found by GM",
                :Time => "Time (s)",
                :Nb_of_feasible_points_found => "Number of feasible points found by GM",
                :Nb_of_maxTime_reached => L"Number of $maxTime$ limit is reached",
                :Nb_of_maxTrials_reached => L"Number of $maxTrial$ limit reached",
                :Avg_Ratio_Non_Bin_Var => "Average ratio of fractional variables (%)",
                :Supported_Sol => "Number of feasible solutions found by GM"
)

const characteristicToFilename = Dict{Symbol,String}(
                :Total_number_of_nd_points => "TotalNbNdPoints",
                :Avg_Ratio_Non_Bin_Var => "AverageNonBinVarPurified05mai"
)

#=
If characteristicBis is null only the values for the given characteristic is given, 
the ratio of characteristic/characteristicBis is ploted otherwise.
=#
function plotQualities(filename::String, characteristic::Symbol, characteristicBis=nothing)
    # parsing 
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refPath*"0.0-0.0.csv"))
    charRefDf = DataFrame(CSV.file(refPath*"AverageNonBinVarPurified05mai.csv"))
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

function plotQualitiesScatterLine(characteristic::Symbol, order::String="nothing", characteristicBis=nothing)
    FONTSIZELABEL = 16
    FONTSIZETICK = 13
    # parsing and getting info
    refDf = DataFrame(CSV.File(refPath*"0.0-0.0.csv"))
    refColor = "black"
    uniqueColor = "royalblue"
    #df = DataFrame()
    name = ""
    titlefile = ""
    if order == "cons"
        name = " sorted by ascending number of constraints"
        titlefile = "bycons"
        refDf[:,"Number of constraints"] = refDfExact[:,"Number of constraints"]
        refDf = sort(refDf,"Number of constraints")
    elseif order == "var"
        name = " sorted by ascending number of variables"
        titlefile = "byvar"
        refDf[:,"Number of variables"] = refDfExact[:,"Number of variables"]
        refDf = sort(refDf,"Number of variables")
    end

    for subdir in subDir
        characteristicRef = (characteristicBis==nothing ? refDf[!,characteristic] : (refDf[!,characteristic] ./ refDf[!,characteristicBis]))
        #characteristicCurrent = (characteristicBis==nothing ? df[!,characteristic] : (df[!,characteristic] ./ df[!,characteristicBis]))

        names = refDf[!,:Instance]
        xPos = collect(1:length(names))

        width = 0.4
        fig, ax1 = plt.subplots(figsize=(20,13),dpi=100)

        #ax2 = ax1.twinx() # second axes
        #fig.set_dpi(300)
        #fig.set_size_inches(20.,13.)
        ax1.tick_params(labelsize=6)
        ax1.set_xlim(xmax=42)
        # Setting the y_axis with a personalised density of ticks
        nbpoints = 17
        
        maxY = 0
        title::String = (characteristicBis==nothing ? characteristicToName[characteristic] : characteristicToName[characteristic]*"/"*characteristicToName[characteristic])
        plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)),fontsize=FONTSIZETICK)#,fontsize=5)
        plt.xlabel("Instances"*name,fontsize=FONTSIZELABEL,weight="bold")
        plt.ylabel(title,fontsize=FONTSIZELABEL,weight="bold")
        #plt.title("Comparison for different vale of"*L"$\alpha$"*title)

        secax = nothing # setting the secondary axis
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(xPos)
        if order == "cons"
            ax2.set_xticklabels(refDf[!,"Number of constraints"],rotation=45,fontsize=FONTSIZETICK)
            ax2.set_xlabel("Number of constraints",fontsize=FONTSIZELABEL,weight="bold")
        elseif order == "var"
            ax2.set_xticklabels(refDf[!,"Number of variables"],rotation=45,fontsize=FONTSIZETICK)
            ax2.set_xlabel("Number of variables",fontsize=FONTSIZELABEL,weight="bold")
        end

        colorIdx = 1 
        colorSet = ["crimson","darkgoldenrod","darkorange","darkviolet","green","royalblue"]
        alphas = collect(0.:0.25:1.)
        
        if (characteristic in [:Avg_Ratio_Non_Bin_Var,:Total_number_of_nd_points])
            println("Path :", refPath*characteristicToFilename[characteristic]*".csv")
            charRefDf = DataFrame(CSV.File(refPath*characteristicToFilename[characteristic]*".csv")) # Extracting the DataFrame from the current file
            if order == "cons"
                charRefDf[:,"Number of constraints"] = refDfExact[:,"Number of constraints"]
                charRefDf = sort(charRefDf,"Number of constraints")
            elseif order == "var"
                charRefDf[:,"Number of variables"] = refDfExact[:,"Number of variables"]
                charRefDf = sort(charRefDf,"Number of variables")
            end
            #println(charRefDf)
            characteristicCurrent = charRefDf[!,characteristic] # characteristic which should be considered as ratio
            maxY = maximum(characteristicCurrent)
            if characteristic == :Avg_Ratio_Non_Bin_Var
                meanRatioFeasibleGenerators = sum(charRefDf[!,:RatioFeasibleGenerators])/length(charRefDf[!,:RatioFeasibleGenerators])
                ax3 = ax1.twinx()
                ax3.set_ylabel("Feasible initial generators (%)", color = "red",fontsize=FONTSIZELABEL,weight="bold")
                ax3.bar(xPos,charRefDf[!,:RatioFeasibleGenerators],width=width,label=L"$\mu$",color="red")
                ax3.plot([xPos[1]-width,xPos[end]+width],[meanRatioFeasibleGenerators,meanRatioFeasibleGenerators],color="red",linestyle="dashed",label=L"Mean $\mu$")
            end
            plt.plot(xPos ,characteristicCurrent,color=uniqueColor) # CURVES with specific color to each characteristic
            plt.scatter(xPos ,characteristicCurrent,color=uniqueColor,label=(characteristic==:Avg_Ratio_Non_Bin_Var ? L"$\eta$" : L"$|Y_N|$")) # SCATTER with specific color to each characteristic
            meanCurrent = sum(characteristicCurrent)/length(characteristicCurrent) # Average value of the current characteristic over the instances
            #stdDev = 100*(characteristicCurrent .- characteristicRef ./ characteristicRef) 
            #pval = pvalue(HypothesisTests.SignedRankTest(characteristicCurrent,characteristicRef)) # Gives the p-value from wilcoxon test
            plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color=uniqueColor,linestyle="dashed",label="Mean "*(characteristic==:Avg_Ratio_Non_Bin_Var ? L"$\eta$" : L"$|Y_N|$"))#*characteristicToName[characteristic])
            
        else 
            # Gives a color to each parametized variant 
            # - - - 
            colorIdx += (subdir=="NonMILPProjection/" ? 1 : 0)  
            plt.plot(xPos ,characteristicRef,color=refColor)
            plt.scatter(xPos ,characteristicRef,label="GM ref",color=refColor)
            meanRef = sum(characteristicRef)/length(characteristicRef)
            plt.plot([xPos[1]-width,xPos[end]+width],[meanRef,meanRef],color=refColor,linestyle="dashed",label="Mean-GM ref")    
            for filename in [file for file in readdir(path*subdir) if (file[end-3:end]==".csv") && (file!="ref.csv")]
                df = DataFrame(CSV.File(path*subdir*filename)) # Extracting the DataFrame from the current file
                if order == "cons"
                    df[:,"Number of constraints"] = refDfExact[:,"Number of constraints"]
                    df = sort(df,"Number of constraints")
                elseif order == "var"
                    df[:,"Number of variables"] = refDfExact[:,"Number of variables"]
                    df = sort(df,"Number of variables")
                end
                characteristicCurrent = (characteristicBis==nothing ? df[!,characteristic] : (df[!,characteristic] ./ df[!,characteristicBis])) # characteristic which should be considered as ratio
                maxY = max(maximum(characteristicRef),maximum(characteristicCurrent),maxY)
                plt.plot(xPos ,characteristicCurrent,color=colorSet[colorIdx]) # CURVES with specific color to each characteristic
                plt.scatter(xPos ,characteristicCurrent,label=L"$\alpha=$"*string(alphas[colorIdx]),color=colorSet[colorIdx]) # SCATTER with specific color to each characteristic
                meanCurrent = sum(characteristicCurrent)/length(characteristicCurrent) # Average value of the current characteristic over the instances
                #stdDev = 100*(characteristicCurrent .- characteristicRef ./ characteristicRef) 
                pval = pvalue(HypothesisTests.SignedRankTest(characteristicCurrent,characteristicRef)) # Gives the p-value from wilcoxon test
                stattest_passed = pval < 0.05
                println("Statistical test passed? -->", stattest_passed)
                plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color=colorSet[colorIdx],linestyle="dashed",
                                                                                        label="Mean "*L"$\alpha=$"*string(alphas[colorIdx])*(stattest_passed ? "✓" : ""))
                println(characteristicToName[characteristic]*"---"*filename[1:end-4]*" - "*string(pval))
                colorIdx += 1
            end 
        end
        if maxY>0
            yticks = range(start=0.,step=maxY/nbpoints,stop=maxY)
            #plt.yticks(yticks, labels=[string(round(v;digits=2)) for v in yticks])
            if characteristic == :Avg_Ratio_Non_Bin_Var
                ax3.set_yticks(yticks, labels = string.(round.(range(start=0.,step=100/nbpoints,stop=100);digits=2)), color = "red")
                ax1.set_ylabel(title,fontsize=FONTSIZELABEL,weight="bold",color=uniqueColor)
                ax1.set_yticks(yticks, labels=[string(round(v;digits=2)) for v in yticks],fontsize=FONTSIZETICK, color=uniqueColor) 
            else
                ax1.set_yticks(yticks, labels=[string(round(v;digits=2)) for v in yticks],fontsize=FONTSIZETICK) 
            end
        end
        #ax2.set_yticks([meanCurrent,meanRef],["Means for SOTA","Mean new method"])
        #plt.yticks([meanRef],["Mean-New Method"],color="orange")
        #plt.yticks([meanCurrent],["Mean_State of the Art"],color="red")
        #plt.legend(bbox_to_anchor=(1.05,1),loc="center right",borderaxespad=0) # legend outside of the plot
        if characteristic == :Avg_Ratio_Non_Bin_Var
            plt.legend(bbox_to_anchor=(1.0, 1.1), loc="upper left")
        else
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")
        end
        #plt.show()
        #println("Saving path : ", filename[1:end-5]*"/"*characteristicsToName[characteristic]*".png")
        plt.subplots_adjust(left=0.065,right=0.91,top=0.91,bottom=0.06)
        plt.savefig(path*subdir*title*titlefile*"scatline.png")
        plt.close(fig)
    end
end

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
                plt.title("Instances "*dir[end-5:end-4]*"-Nb variables in R before / Nb of variables in R after")
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
        for order in ["cons","var"]
            plotQualitiesScatterLine(characteristic,order)
        end
    end
end

main()