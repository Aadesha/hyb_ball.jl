using TaylorModels

function lbd_1(f::Function, dom::IntervalBox{N}, order::Int) where {N}
    x0 = mid(dom)
    n = length(dom)
    set_variables(Float64, "x", order=2order, numvars=n)
    x = [TaylorModelN(i, order, IntervalBox(x0), dom) for i=1:n]
    return f(x...)
end

function lbd_2(TM::TaylorModelN)  (check again)
    domain = TM.dom
    if TM.pol.coeffs[1].coeffs[1] < 0
        TM.pol.coeffs[1].coeffs[1] =  -TM.pol.coeffs[1].coeffs[1]
    for i = 1:length(TM.pol.coeffs[2])
    # if i%2 == 0
        if TM.pol.coeffs[2].coeffs[i] < 0
            TM.pol.coeffs[2].coeffs[i] = -TM.pol.coeffs[2].coeffs[i] #for linear models
        end
    end
    return TM
end

function lbd_3(TM::TayolorModelN)
    I1 = linear_polynomial(TM.pol) + constant_term(TM.pol)
    Ih = TM.pol - I1
    bound = evaluate(Ih, TM.dom) + (evaluate(I1, TM.dom)).lo
    return bound
end
                                                                 #

function lbd_4(d::Interval,ϵ::Number,dom::IntervalBox,TM::TaylorModelN)
    if d <= ϵ
        return d
    else
        l = TM.pol.coeffs[2]
        for i = 1: length(TM.dom)
            if l[i] > 0 && dom[i] > d/l[i]
               dom = setindex(dom, dom[i] + d/l[i], i)
            end
        end
    end
    return dom
end

