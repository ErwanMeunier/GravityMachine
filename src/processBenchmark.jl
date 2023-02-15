using DataFrames
using CSV
using LaTeXStrings
import PyPlot
const plt = PyPlot
plt.pygui(true)
const path = "./results/"
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
#=
If characteristicBis is null only the values for the fiven characteristic is given, 
the ratio of characteristic/characteristicBis is ploted otherwise.
=#
function plotQualities(filename::String, characteristic::Symbol, refFile::String=path*"resultRef.csv", characteristicBis=nothing)
    # parsing 
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refFile))
    # getting some informations
    names = df[!,:Instance]
    characteristicCurrent = df[!,characteristic]
    characteristicRef = refDf[!,characteristic]
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
    

    plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)))#,fontsize=5)
    plt.xlabel("Instance")
    plt.ylabel(string(characteristic))
    plt.title("Comparison between GM ref and "*filename[1:end-3]*"-"*string(characteristic))
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
    plt.savefig(filename[1:end-4]*"/"*string(characteristic)*".png")
    plt.close(fig)
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
            try
                plotQualities(path*file, characteristic)
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

function main()
    plotRatioMIP()
end

main()