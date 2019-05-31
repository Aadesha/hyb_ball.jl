using Revise, TaylorModels, Plots, TaylorSeries, IntervalArithmetic
using LazySets:Hyperrectangle,split

function branchandbound(P::TaylorN{T},dom::IntervalBox{N,T},ϵ::Number) where {N,T}
    K = 1
    Rperv = evaluate(p,dom)
    D1,D2 = bisect(dom)
    D = [ D1 , D2 ]
    R = [evaluate(p, D[i]) for i = 1:length(D)]
    Rnext = Interval(minimum(R[i].lo for i = 1:length(R)),
                     maximum(R[i].hi for i = 1:length(R)))
    while  (Rperv.hi - Rnext.hi) <= ϵ*(Rnext.hi - Rnext.lo) &&
           (Rperv.lo - Rnext.lo) <= ϵ*(Rnext.hi - Rnext.lo) && (K <= 1000)
        Rperv = Interval(minimum(R[i].lo for i = 1:length(R)),
                         maximum(R[i].hi for i = 1:length(R)))
        R_x = [ R[i].hi for i = 1:length(R)]
        R_n = [ R[i].lo for i = 1:length(R)]
        max = maximum(R[i].hi for i = 1:length(R))
        max = findall(x->x == max, R_x)[1]
        min = minimum(R[i].lo for i = 1:length(R))
        min = findall(x->x == min, R_n)[1]
        l_D = length(D[1])
        println(D)
        println(R)
        K = K + 1
        if min == max
            BA = [ ((D[max][i]).hi - (D[max][i]).lo) for i = 1:l_D]
            Beta1 = maximum(BA)
            β = findall(x->x==Beta1,BA)
            Mat = ones(Int64,length(BA)); Mat[β] = 2
            H = convert(Hyperrectangle,D[max])
            D1,D2 = split(H,Mat)
            D[max] = convert(IntervalBox, D1)
            DD = append!(D[1:max],convert(IntervalBox, D2))
            D = append!(DD, D[(max+1):length(D)])
            println("\n")
            println(max)
            println(D)
            println("1")
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
            println(R)
        else
            BA = [ ((D[max][i]).hi - (D[max][i]).lo) for i = 1:l_D]
            Beta1 = maximum(BA)
            β = findall(x->x==Beta1,BA)[1]
            Mat = ones(Int64,length(BA)); Mat[β] = 2
            H = convert(Hyperrectangle,D[max])
            D1,D2 = split(H,Mat)
            D[max] = convert(IntervalBox, D1)
            DD = push!(D[1:max],convert(IntervalBox, D2)) #gives error
            D = append!(DD, D[(max+1):length(D)])
            println("\n")
            println(max)
            println(D)
            println("2")
            R[max] = evaluate(p,D[max])
            RR = append!(R[1:max], evaluate(p,D[max+1]))
            R = append!(RR, R[(max+1):length(R)])
            println(R)
            if max < min
                min = min + 1
            end
            BA = [ ((D[min][i]).hi - (D[min][i]).lo) for i = 1:l_D]
            Beta1 = maximum(BA)
            β = findall(x->x==Beta1,BA)[1]
            Mat = ones(Int64,length(BA)); Mat[β] = 2
            H = convert(Hyperrectangle,D[min])
            D1,D2 = split(H,Mat)
            D[min] = convert(IntervalBox, D1)
            DD = push!(D[1:min],convert(IntervalBox, D2))
            D = append!(DD, D[(min+1):length(D)])
            println("\n")
            println(min)
            println(D)
            println("3")
            R[min] = evaluate(p,D[min])
            RR = append!(R[1:min], evaluate(p,D[min+1]))
            R = append!(RR, R[(min+1):length(R)])
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
            println(R)
        end
    end
    println(K)
    return Interval(minimum(R[i].lo for i=1:length(R)),
                    maximum(R[i].hi for i=1:length(R)))
end

m = 4
x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=2*m)
p = 1 - x₂ + 2*x₁
Dx₁ = Interval(0.0, 3.0)
Dx₂ = Interval(-1.0, 1.0)
D = Dx₁ × Dx₂

println(branchandbound(p,D,0.1))
println(evaluate(p,D))
