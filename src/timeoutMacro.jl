# from : https://discourse.julialang.org/t/help-writing-a-timeout-macro/16591/7
macro timeout(seconds, expr, fail)
    quote
        tsk = @task $esc(expr)
        schedule(tsk)
        Timer($(esc(seconds))) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $(esc(fail))
        end
    end
end

function testTimeout(x,a,b)
    sleep(x)
    return a + b
end

function main()
    s = 10
    a = 1
    b = 2
    failValue = 0
    @timeout s testTimeout(a,b) failValue
end

