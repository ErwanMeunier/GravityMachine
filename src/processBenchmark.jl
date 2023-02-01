using DataFrames
using CSV
import PyPlot
const plt = PyPlot
plt.pygui(true)

function plotQualities(filename::String,refFile::String="resultRef.csv")
    
    df = DataFrame(CSV.File(filename))
    refDf = DataFrame(CSV.File(refFile))
    names = df[!,:Instance]
    qualitiesCurrent = df[!,:Quality]
    qualitiesRef = refDf[!,:Quality]
    println(names)
    # Figure parameters
    fig, ax = plt.subplots()
    #fig.set_size_inches(18.5,10.5)
    fig.set_dpi(100)
    xPos = collect(length(names))
    width = 0.1
    # ---
    rect1 = ax.bar(xPos .- width/2, qualitiesCurrent, width, label=filename[1:end-5])#,tick_label=map(x->round(x),qualities))
    rect2 = ax.bar(xPos .+ width/2, qualitiesRef, width, label=refFile[1:end-5])
    
    ax.set_xlabel("Instance")
    ax.set_ylabel("Quality")
    ax.set_title("Comparison between GM ref and "*filename[1:end-5])
    ax.set_xticks(xPos, collect(map(x->x[end-5:end-4],names)))
    ax.legend()

    ax.bar_label(rect1, padding=3)
    ax.bar_label(rect2, padding=3)
    
    #ax.set_ylim(0,101)
    #ax.set_locator_params(axis="y", nbins=51)
    fig.tight_layout()
    plt.show()
    plt.savefig(filename*".png",)
    plt.close(fig)
end

csvf = [file for file in readdir() if file[end-3:end]==".csv"]

for file in csvf
    plotQualities(file)
end