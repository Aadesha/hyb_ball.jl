
using Revise, TaylorModels, Plots, TaylorSeries, IntervalArithmetic



function breanchandbound(P::Taylor1{T},dom::Interval{T},ϵ::Number)where{T}
    Rperv = evaluate(p,dom)
    H = (Float64(dom.hi) + Float64(dom.lo))/2
    D = [Interval(dom.lo,H),Interval(H,dom.hi)]
    println(D);println("1")
    R = [evaluate(p, D[i]) for i = 1:length(D)]

    println(R)

    Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                     maximum(R[i].hi for i=1:length(R)))
    K = 1
    while  (Rperv.hi - Rnext.hi) <= ϵ*(Rnext.hi - Rnext.lo) && (Rperv.lo - Rnext.lo) <= ϵ*(Rnext.hi - Rnext.lo) && (K <= 4)
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
        println(" \n")
        println(min)
        println(max)
        println(" \n")

        K = K + 1
        if min == max
            H = (D[max].hi + D[max].lo)/2
            D1 = Interval(D[max].lo,H)
            D2 = Interval(H,D[max].hi)
            D[max] = D1
            DD = append!(D[1:max], D2)
            D = append!(DD, D[(max+1):length(D)])
            println(" \n")
            println(D)
            println("2")
            println(" \n")
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            println(" r1r\n")
            println(R)
            println(" r2r\n")
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
            #R = R[2:length(R)]
        else
            H = (D[max].hi + D[max].lo)/2
            D1 = Interval(D[max].lo,H)
            D2 = Interval(H,D[max].hi)
            D[max] = D1
            DD = append!(D[1:max], D2)
            println("DD")
            println(DD)
            println("DD")
            D = append!(DD, D[(max+1):length(D)])
            println(" \n")
            println(D)
            println("3")
            println(" \n")
            #R = R[2:length(R)]
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            println(" r1r\n")
            println(R)
            println("r2r \n")
            if max < min
                min = min + 1
            end
            H = (D[min].hi + D[min].lo)/2
            D1 = Interval(D[min].lo,H)
            D2 = Interval(H,D[min].hi)
            D[min] = D1
            DD = append!(D[1:min], D2)
            println("DD")
            println(DD)
            println("DD")
            if (min + 1) == length(D)
                D = append!(DD, D[min+1])
            else
                D = append!(DD, D[(min+1):length(D)])
            end
            println(" \n")
            println(D)
            println("4")
            println(" \n")
            #R = R[2:length(R)]
            R[min] = evaluate(p,D[min])
            RR = append!(R[1:min], evaluate(p,D[min+1]))
            R = append!(RR, R[(min+1):length(R)])
            println(" r1r\n")
            println(R)
            println(" r2r\n")
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                              maximum(R[i].hi for i=1:length(R)))
        end
    end
    return Interval(minimum(R[i].lo for i=1:length(R)),
                            maximum(R[i].hi for i=1:length(R)))
end
p = Taylor1([0.0, 1.0, 1.0],3)
D = Interval(-3.0, 1.0)
ϵ=10
println(breanchandbound(p,D,ϵ))
