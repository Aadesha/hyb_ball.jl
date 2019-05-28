

function breanchandbound(P::Taylor1{T},dom::Interval{T},ϵ::Number)
    D = bisect(dom)
    D = [D[i][1] for i = 1:length(D)]
    R = [evaluate(p, D[i]) for i = 1:length(D)]
    while #conditon
        Rperv = Interval(minimum(R[i].lo for i=1:length(R)),
                         maximum(R[i].hi for i=1:length(R)))
        max = maximum(R)
        for i = 1:length(R)
            if R[i] == max
                max = i
                break;
            end
        end
        min = minimum(R)
        for i = 1:length(R)
            if R[i] == min
                min = i
                break;
            end
        end

        if min == max
            D1, D2 = bisect(D[max])
            D[max] = D1
            DD = append!(D[1:max], D2)
            D = append!(DD, D[(max+1):langth(D)])
            R = append!(R, evaluate(p, D[max]))
            R = append!(R, evaluate(p, D[max+1]))
        else
            D1, D2 = bisect(D[max])
            D[max] = D1
            DD = append!(D[1:max], D2)
            D = append!(DD, D[(max+1):langth(D)])
            R = append!(R, evaluate(p, D[max]))
            R = append!(R, evaluate(p, D[max+1]))
            D1, D2 = bisect(D[min])
            D[min] = D1
            DD = append!(D[1:min], D2)
            D = append!(DD, D[(min+1):langth(D)])
            R = append!(R, evaluate(p, D[min]))
            R = append!(R, evaluate(p, D[min+1]))
        end
    end
    return Interval(minimum(R[i].lo for i=1:length(R)),
                     maximum(R[i].hi for i=1:length(R)))
end


          
