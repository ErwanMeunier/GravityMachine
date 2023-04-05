include("GMmain.jl")
using DataFrames, CSV


const instanceDirectory = "../SPA/instances/"
const tailleSampling = 6
const maxTrial = 20
const maxTime = 20

function benchmarkGM()
    filenames::Vector{String} = readdir(instanceDirectory)
    for max_ratio_bv_pr in 1:0.25:1
        for max_ratio_bv_gi in (1:0.25:1)
            qualities = Vector{Float64}()
            nbcyclestotalVect = Vector{Int64}()
            nbcyclesMaxVect = Vector{Int64}()
            tpns = Vector{Int64}()
            fpns = Vector{Int64}()
            times = Vector{Float64}()
            feasibleReached = Vector{Int64}()
            maxTimeReached = Vector{Int64}()
            maxTrialsReached = Vector{Int64}()
            
            for instance in filenames
                    etime, quality, nbcyclestotal, nbcyclesMax, tpn, fpn, nbFeasible, nbMaxTime, nbMaxTrials = GM(instance[4:end], tailleSampling, maxTrial, maxTime, max_ratio_bv_gi, max_ratio_bv_pr)
                    push!(qualities,quality)
                    push!(nbcyclestotalVect,nbcyclestotal)
                    push!(nbcyclesMaxVect,nbcyclesMax)
                    push!(tpns,tpn)
                    push!(fpns,fpn)
                    push!(times,etime)
                    push!(feasibleReached,nbFeasible)
                    push!(maxTimeReached,nbMaxTime)
                    push!(maxTrialsReached,nbMaxTrials)
                    #=
                    For each column of the dataframe above, we have the evolution of the ratio of the number of 
                    floating point variables in the projection.
                    =#
            end
            output = DataFrame()
            #println(filenames)
            #println(qualities)
            #println(nbcyclestotalVect)
            #println(nbcyclesMaxVect)
            #println(tpns)
            #println(fpns)
            #println(times)
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
            CSV.write("./results/resultsBinVar/"*string(max_ratio_bv_pr)*"-"*string(max_ratio_bv_gi)*".csv",output;delim=";")
        end
    end 
end
#=
function test()
    filenames::Vector{String} = readdir(instanceDirectory)
    for instance in filenames 
        loadInstance2SPA(instance[4:end])
    end
end
=#
#=
#@time GM("sppaa02.txt", 6, 20, 20)
#@time GM("sppnw03.txt", 6, 20, 20) #pb glpk
#@time GM("sppnw10.txt", 6, 20, 20)
#@time GM("didactic5.txt", 5, 5, 10)
@time GM("sppnw29.txt", 6, 30, 20)
nothing

function GM( fname::String,
    tailleSampling::Int64,
    maxTrial::Int64,
    maxTime::Int64
  )

=#

# Output:
# instance's name | quality | number of cycle in total | maximum number of cycles for each cone | "true number" of non-dominated points | "estimated number of non-dominated points"
