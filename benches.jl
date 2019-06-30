# dependencies
using BenchmarkTools
using SumOfSquares, DynamicPolynomials, SemialgebraicSets, MosekTools, SDPA
using Plots
using TaylorModels
#using AffineArithmetic

# domains
a = Interval(-4.5, -0.3)
b = Interval(0.4, 0.9)
c = Interval(3.8, 7.8)
d = Interval(8.0, 10.0)
e = Interval(-10.0, 8.0)
f = Interval(1.0, 2.0)

# fix the order for the series expansions
REF = Dict()
RELPREC = Dict()
SUITE["Daisy"] = BenchmarkGroup()
RESULTS = SUITE["Daisy"]

REF["sin"] = Interval(-1.0, 0.977530117665097)
REF["bspline0"] = Interval(0.36616666666666675, 27.729166666666668)
REF["bspline1"] = Interval(-65.14583333333333, 0.5631666666666667)
REF["bspline2"] = Interval(0.07516666666666667, 53.604166666666664)
REF["bspline3"] = Interval(0.0045, 15.1875)
REF["himmilbeau"] = Interval(85.46830677734748, 221.7338939301446)
REF["kepler1"] = Interval(-5.255935494810441, 7.321362422825775)
REF["kepler2"] = Interval(-195.36974909125482, 78.3669520644375)
REF["kepler3"] = Interval(-309.8484155131222, 17.982082401462407)
REF["rigidbody1"] = Interval(-20.786552979420335, -0.540012836551535)
REF["rigidbody2"] = Interval(68.81138021006673, 359.98566570476504)

# This function measures the relative precision of the result in a more informative way than
# taking the scalar overestimation because it evaluates the precision of the lower and the
# upper range bounds separately, see Eq. (20) in [1].
function relative_precision(x, x_ref)
    x_low, x_high = inf(x), sup(x)
    x_ref_low, x_ref_high = inf(x_ref), sup(x_ref)
    rel_low = -(x_low - x_ref_low) / (x_ref_high - x_ref_low)
    rel_high = (x_high - x_ref_high) / (x_ref_high - x_ref_low)
    return 100 * Interval(rel_low, rel_high)
end

const available_algorithms = [ "IA", "TM", "TM_NORM", "ZOO"];

# This function returns an overapproximation of the range of `f`
# using different methods
function bounds(f::Function, dom::Union{Interval, IntervalBox}; algorithm="ZOO",
                                                order=2,name="sin")::Interval

    numvars = length(dom)
    if algorithm == "IA"
        bnd = _bounds_IA(f, dom)
        RESULTS[name] = BenchmarkGroup()
        RESULTS[name,order,algorithm] = @benchmarkable $_bounds_IA(f, dom)
    elseif algorithm == "TM"
        bnd, = _bounds_TM(f, dom, order)
        RESULTS[name] = BenchmarkGroup()
        RESULTS[name,order,algorithm] = @benchmarkable $_bounds_TM(f, dom)
    elseif algorithm == "TM_NORM"
        bnd = _bounds_TM_NORM(f, dom, order)
        RESULTS[name] = BenchmarkGroup()
        RESULTS[name,order,algorithm] = @benchmarkable $_bounds_TM_NORM(f, dom, order)
    elseif algorithm == "ZOO"
        algos = setdiff(available_algorithms, ["ZOO"])
        bnd_all = [bounds(f, dom, algorithm=a, order=order) for a in algos]
        bnd = reduce(∩, bnd_all)
    else
        error("algorithm $algorithm unknown")
    end
    return bnd
end

# interval arithmetic substitution
function _bounds_IA(f, dom)
    return f(dom...)
end

# affine arithmetic substitution
function _bounds_AA(f::Function, dom::Interval)
    x = AFF(dom, 1, 1)
    return interval(f(x))
end

# affine arithmetic in N variables
function _bounds_AA(f::Function, dom::IntervalBox{N}) where {N}
    x = [AFF(dom[i], N, i) for i in 1:N]
    return interval(f(x...))
end

# taylor model in one variable
function _bounds_TM(f::Function, dom::Interval, order::Int)
    x0 = Interval(mid(dom))
    x = TaylorModel1(order, x0, dom)
    return evaluate(f(x), dom - x0)
end

# normalized taylor model in one variable
function _bounds_TM_NORM(f::Function, dom::Interval, order::Int)
    x0 = Interval(mid(dom))
    x = TaylorModel1(order, x0, dom)
    xnorm = normalize_taylor(x.pol, dom - x0, true)
    xnormTM = TaylorModel1(xnorm, 0..0, 0..0, -1..1)
    return evaluate(f(xnormTM), -1..1)
end

# taylor model in N variables
function _bounds_TM(f::Function, dom::IntervalBox{N}, order) where {N}
    x0 = mid(dom)
    set_variables(Float64, "x", order=2order, numvars=N)
    x = [TaylorModelN(i, order, IntervalBox(x0), dom) for i=1:N]
    return evaluate(f(x...), dom - x0)
end

# normalized taylor model in N variables
function _bounds_TM_NORM(f::Function, dom::IntervalBox{N}, order::Int) where {N}
    x0 = mid(dom)
    set_variables(Float64, "x", order=2order, numvars=N)

    zeroBox = IntervalBox(0..0, N)
    symBox = IntervalBox(-1..1, N)

    x = [TaylorModelN(i, order, IntervalBox(x0), dom) for i=1:N]
    xnorm = [normalize_taylor(xi.pol, dom - x0, true) for xi in x]
    xnormTM = [TaylorModelN(xi_norm, 0..0, zeroBox, symBox) for xi_norm in xnorm]
    return evaluate(f(xnormTM...), symBox)
end


function SDSF(PP::Function)
    for ALGO in available_algorithms, m in [2,5,10]
        bnd = bounds(PP, dom, algorithm=ALGO,order = m)
        RELPREC[model_name,m,ALGO] = relative_precision(bnd, REF[model_name])
    end
end

model_name = "sin"
dom = a
SDSF(sin)

model_name = "bspline0"
dom = a
bspline0(x) = (1 - x) * (1 - x) * (1 - x) / 6.0
SDSF(bspline0)

model_name = "bspline1"
dom = a
bspline1(x) = (3*x*x*x - 6*x*x + 4) / 6.0
SDSF(bspline1)

model_name = "bspline2"
dom = a
bspline2(x) = (-3*x*x*x  + 3*x*x + 3*x + 1) / 6.0
SDSF(bspline2)

model_name = "bspline3"
dom = a
bspline3(x) = -x*x*x / 6.0
SDSF(bspline3)

model_name = "himmilbeau"
dom = a × b

himmilbeau(x1, x2) = (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) +
                                 (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)
SDSF(himmilbeau)

model_name = "kepler1"
dom = a × b × c × d × e × f
kepler1(x1, x2, x3, x4, x5, x6) = x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 +
                                         x1 * (-x1 + x2 + x3 - x4 + x5 + x6)
SDSF(kepler1)

model_name = "kepler2"
dom = a × b × c × d
kepler2(x1, x2, x3, x4) = x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4) +
                          x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4
SDSF(kepler2)

model_name = "kepler3"
dom = a × b × c × d × e × f
kepler3(x1, x2, x3, x4, x5, x6) = x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6) +
                                  x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) + x3 * x6 * (x1 + x2 - x3 + x4 + x5 - x6) -
                                  x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6
SDSF(kepler3)

model_name = "rigidbody1"
dom = a × b × c
rigidbody1(x1, x2, x3) = -x1*x2 - 2*x2*x3 - x1 - x3
SDSF(rigidbody1)

model_name = "rigidbody1"
dom = a × b × c
rigidbody2(x1, x2, x3) = 2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2
SDSF(rigidbody2)
