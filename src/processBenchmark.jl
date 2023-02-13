using DataFrames
using CSV
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

#=
output[!, :Instance]=filenames
    output[!, :Quality]=qualities
    output[!, :Number_Of_Cycles]=nbcyclestotalVect
    output[!, :Max_Number_Of_Cycles]=nbcyclesMaxVect
    output[!, :Total_number_of_nd_points]=tpns
    output[!, :Total_number_of_nd_points_found_by_GM]=fpns
    output[!, :Time]=times
    output[!, :Nb_of_feasible_points_found]=feasibleReached
    output[!, :Nb_of_maxTime_reached]=maxTimeReached
    output[!, :Nb_of_maxTrials_reached]=maxTrialsReached
=#

function plotQualities(filename::String, characteristic::Symbol, refFile::String=path*"resultRef.csv")
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
    plt.bar(xPos .+ width/2,characteristicRef,width=width/2,label="State-of-the-art results",color="green",)
    yticks = [10*i for i=0:Int(ceil(max(maximum(characteristicRef),maximum(characteristicCurrent))/10))]
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
    plt.savefig(filename[1:end-5]*"/"*string(characteristic)*".png")
    plt.close(fig)
end

function main()
    csvf::Vector{String} = [file for file in readdir(path) if file[end-3:end]==".csv"]

    println("Files :", csvf)

    for file in csvf
        try
            mkdir(path*file[1:end-5])
        catch; Base.IOError # the directory already exists
            nothing 
        end

        for characteristic in fields
            try
                plotQualities(path*file, characteristic)
            catch ArgumentError # compatibility with the former version of GMBenchmark which includes less measures
                nothing
            end
        end
    end
end

main()