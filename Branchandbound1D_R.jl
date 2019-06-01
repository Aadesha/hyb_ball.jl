using Revise, TaylorModels, Plots, TaylorSeries, IntervalArithmetic, IntervalOptimisation


function branchandbound(p::Taylor1{T},dom::Interval{T},ϵ::Number)where{T}
    Rperv = evaluate(p,dom)
    D = bisect(dom)
    D = [ D[i][1] for i = 1:2]
    #println(D);
    R = [evaluate(p, D[i]) for i = 1:length(D)]
    Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                     maximum(R[i].hi for i=1:length(R)))
    println(Rperv)
    K = 1
    while  (Rperv.hi - Rnext.hi) <= ϵ*(Rnext.hi - Rnext.lo) &&
           (Rperv.lo - Rnext.lo) <= ϵ*(Rnext.hi - Rnext.lo) && (K <= 100)
        Rperv = Interval(minimum(R[i].lo for i=1:length(R)),
                         maximum(R[i].hi for i=1:length(R)))
        max = maximum(R[i].hi for i=1:length(R))
        for i = 1:length(R)
            if R[i].hi == max
                max = i
                break;
            end
        end
        min = minimum(R[i].lo for i=1:length(R))
        for i = 1:length(R)
            if R[i].lo == min
                min = i
                break;
            end
        end
        K = K + 1
        if min == max
            D1, D2 = bisect(D[max])
            D[max] = D1
            DD = append!(D[1:max], D2)
            D = append!(DD, D[(max+1):length(D)])
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
        else
            D1, D2 = bisect(D[max])
            D[max] = D1
            DD = append!(D[1:max], D2)
            D = append!(DD, D[(max+1):length(D)])
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            if max < min
                min = min + 1
            end
            D1, D2 = bisect(D[min])
            D[min] = D1
            DD = append!(D[1:min], D2)
            if (min + 1) == length(D)
                D = append!(DD, D[min+1])
            else
                D = append!(DD, D[(min+1):length(D)])
            end
            R[min] = evaluate(p,D[min])
            RR = append!(R[1:min], evaluate(p,D[min+1]))
            R = append!(RR, R[(min+1):length(R)])
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
        end
    end
    return Interval(minimum(R[i].lo for i=1:length(R)),
                    maximum(R[i].hi for i=1:length(R)))
end

ϵ=0.4

Dx=Interval(-4.5,-0.3)
dom=Dx
x=Taylor1(7)

sin=x - (x*x*x)/6.0 + (x*x*x*x*x)/120.0 - (x*x*x*x*x*x*x)/5040.0

println(branchandbound(sin,Dx,ϵ))

function _minmax(p, dom) #takes time
    global_min, _ = minimise(p, dom)
    minus_global_max, _ = minimise(-p, dom)
    global_max = -minus_global_max
    return global_min, global_max
end
global_min, global_max = _minmax(sin, Dx)

println(global_min)
println(global_max)
