using DataFrames
using CSV
import PyPlot
const plt = PyPlot
plt.pygui(true)

function plotQualities(filename::String,refFile::String="resultRef.csv")
    # parsing 
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refFile))
    # getting some informations
    names = df[!,:Instance]
    qualitiesCurrent = df[!,:Quality]
    qualitiesRef = refDf[!,:Quality]
    println(qualitiesCurrent)
    println(qualitiesRef)
    println(names)
    
    xPos = collect(1:length(names))

    
    width = 0.4
    fig, ax = plt.subplots()
    fig.set_dpi(300)
    fig.set_size_inches(15.,10.)
    ax.tick_params(labelsize=6)
    plt.bar(xPos .- width/2,qualitiesCurrent,width=width/2,label="New method",color="blue",)
    plt.bar(xPos .+ width/2,qualitiesRef,width=width/2,label="State-of-the-art results",color="green",)
    #ax.bar_label(rect1, padding=3)
    #ax.bar_label(rect2, padding=3)
    
    #ax.set_ylim(0,101)
    #ax.set_locator_params(axis="y", nbins=51)
    plt.yticks([10*i for i=0:10],[string(10*i) for i=0:10])
    plt.xticks(xPos, collect(map(x->x[end-5:end-4],names)))#,fontsize=5)
    plt.xlabel("Instance")
    plt.ylabel("Quality")
    plt.title("Comparison between GM ref and "*filename[1:end-4])
    meanRef = sum(qualitiesRef)/length(qualitiesRef)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanRef,meanRef],color="r")
    meanCurrent = sum(qualitiesCurrent)/length(qualitiesCurrent)
    plt.plot([xPos[1]-width,xPos[end]+width],[meanCurrent,meanCurrent],color="orange")
    #plt.yticks([meanRef],["Mean-New Method"],color="orange")
    #plt.yticks([meanCurrent],["Mean_State of the Art"],color="red")
    plt.legend()
    plt.show()
    plt.savefig(filename*".png")
    plt.close(fig)
end

csvf = [file for file in readdir() if file[end-3:end]==".csv"]

for file in csvf
    plotQualities(file)
end