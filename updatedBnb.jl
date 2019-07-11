
using Revise, TaylorModels


function branchandbound(p::TaylorN{T},dom::IntervalBox{N,T},ϵ::Number) where {N,T}
    K = 1
    Rperv = evaluate(p,dom)
    D1,D2 = bisect(dom)
    D = [ D1 , D2 ]
    R = [evaluate(p, D[i]) for i = 1:length(D)]
    Rnext = Interval(minimum(R[i].lo for i = 1:length(R)),
                     maximum(R[i].hi for i = 1:length(R)))
     #condition given in the algo and max no of steps K
    while  (Rperv.hi - Rnext.hi) <= ϵ*(Rnext.hi - Rnext.lo) &&
           (Rperv.lo - Rnext.lo) <= ϵ*(Rnext.hi - Rnext.lo) && (K <= 1000)
        println(D)
        #caclulate the max and min of range before further splitting
        Rperv = Interval(minimum(R[i].lo for i = 1:length(R)),
                         maximum(R[i].hi for i = 1:length(R)))
        R_x = [ R[i].hi for i = 1:length(R)]
        R_n = [ R[i].lo for i = 1:length(R)]
        max_range = maximum(R[i].hi for i = 1:length(R)) # max of range
        max_index = findall(x->x == max_range, R_x)[1]#index of max range giving domain
        min_range = minimum(R[i].lo for i = 1:length(R))# min of range
        min_index = findall(x->x == min_range, R_n)[1]#index of min range giving domain
        l_D = length(D[1])#number of dimension
        K = K + 1
        if min_index == max_index
            devide_dom(p,D,R,max_index)
#=          BA = [ ((D[max_index][i]).hi - (D[max_index][i]).lo) for i = 1:l_D] #width of the each dimention
            Beta1 = maximum(BA)#max element of BA
            β = findall(x->x==Beta1,BA)[1]
            D1,D2 = bisect(D[max_index],β)
            D[max_index] = D1
            DD = push!(D[1:max_index],D2)
            D = append!(DD, D[(max_index+1):length(D)])
            #update the range array
            R[max_index] = evaluate(p,D[max_index])
            RR = append!(R[1:max_index], evaluate(p,D[max_index+1]))
            R = append!(RR, R[(max_index+1):length(R)])  =#
            #max and min of range after split.
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
        else
            devide_dom(p,D,R,max_index)
#=          BA = [ ((D[max_index][i]).hi - (D[max_index][i]).lo) for i = 1:l_D] #width of the each dimention
            Beta1 = maximum(BA) #max element of BA
            β = findall(x->x==Beta1,BA)[1]#Index of split
            D1,D2 = bisect(D[max_index],β)
            D[max_index] = D1
            DD = push!(D[1:max_index],D2)
            D = append!(DD, D[(max_index+1):length(D)])
            #update the range array
            R[max_index] = evaluate(p,D[max_index])
            RR = append!(R[1:max_index], evaluate(p,D[max_index+1]))
            R = append!(RR, R[(max_index+1):length(R)]) =#
            #= if index 'max' occurs before the index 'min' then
            update 'min' as the indexing have changed after the update in the domain array.=#
            if max_index < min_index
                min_index = min_index + 1
            end
            devide_dom(p,D,R,min_index)
#=          BA = [ ((D[min_index][i]).hi - (D[min_index][i]).lo) for i = 1:l_D]#width of the each dimention
            Beta1 = maximum(BA)  #max element of BA
            β = findall(x->x == Beta1,BA)[1]#Index of split
            D1,D2 = bisect(D[min_index],β)
            D[min_index] = D1
            DD = push!(D[1:min_index],D2)
            D = append!(DD, D[(min_index + 1):length(D)])
            R[min_index] = evaluate(p,D[min_index])
            RR = append!(R[1:min_index], evaluate(p,D[min_index+1]))
            R = append!(RR, R[(min_index + 1):length(R)])  =#
            Rnext = Interval(minimum(R[i].lo for i=1:length(R)),
                             maximum(R[i].hi for i=1:length(R)))
        end
    end
    #println(K)
    return Interval(minimum(R[i].lo for i=1:length(R)),
                    maximum(R[i].hi for i=1:length(R)))
end

function devide_dom(p,D,R,index)
    BA = [ ((D[index][i]).hi - (D[index][i]).lo) for i = 1:length(D[1])] #width of the each dimention
    Beta1 = maximum(BA) #max element of BA
    β = findall(x->x==Beta1,BA)[1]#Index of split
    D1,D2 = bisect(D[index],β)
    D[index] = D1
    DD = push!(D[1:index],D2)
    D = append!(DD, D[(index+1):length(D)])
    #update the range array
    R[index] = evaluate(p,D[max_index])
    RR = append!(R[1:index], evaluate(p,D[index+1]))
    R = append!(RR, R[(index+1):length(R)])
end
