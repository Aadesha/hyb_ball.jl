# ==========
# Daisy
# ==========

# The benchmarks in this file are taken from [Project Daisy](https://github.com/malyzajko/daisy/blob/master/testcases/).

SUITE["Daisy"] = BenchmarkGroup()
DAISY = SUITE["Daisy"]

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

# Helper function to evaluate after normalization
normalize_and_evaluate(p, dom) = evaluate(normalize_taylor(p, dom, true), -1..1)

# Dictionary to hold the vector of relative precision intervals for each benchmark
RELPREC = Dict{String, Any}()

# The following intervals are used throughout the tests in Tables 3-5 in
# [1] Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018). Implementation of Taylor models in CORA 2018.
#     In Proc. of the 5th International Workshop on Applied Verification for Continuous and Hybrid Systems.
a = Interval(-4.5, -0.3)
b = Interval(0.4, 0.9)
c = Interval(3.8, 7.8)
d = Interval(8.0, 10.0)
e = Interval(-10.0, 8.0)
f = Interval(1.0, 2.0)

# ======
# sine
# ======

DAISY["sin"] = BenchmarkGroup()
RELPREC["sin"] = Dict()

# reference value for this test
ref = Interval(-0.9998434996853951, 2.724232859130124)

dom = a
for m in [2, 5, 10]
    DAISY["sin"]["order $m"] = BenchmarkGroup()
    RELPREC["sin"]["order $m"] = Dict()

    x = Taylor1(m)                                # Ref should be done again and as we have used the Ipopt for the Opimization? 
    p = sin(x)                                                                 # where is the equation? why function? 

    DAISY["sin"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["sin"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["sin"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["sin"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# bspline0
# ==========

DAISY["bspline0"] = BenchmarkGroup()
RELPREC["bspline0"] = Dict()

dom = a
for m in [2, 5, 10]
    DAISY["bspline0"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline0"]["order $m"] = Dict()
    ref = Interval(0.36616666073000836, 27.729167344785765)
    
    u = Taylor1(m)                                 # again m should be always greater than tha degree of euation? But here for 2 is not true?   
    p = (1 - u) * (1 - u) * (1 - u) / 6.0
   
    DAISY["bspline0"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline0"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline0"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline0"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end    

DAISY["bspline1"] = BenchmarkGroup()
RELPREC["bspline1"] = Dict()
ref = Interval(-65.14583510270246, 0.5631666715107633)
dom = a
for m in [2, 5, 10]
    DAISY["bspline1"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline1"]["order $m"] = Dict()
    p = (3 * u*u*u - 6 * u*u + 4) / 6.0
    
    DAISY["bspline1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end
DAISY["bspline2"] = BenchmarkGroup()
RELPREC["bspline2"] = Dict()
ref = Interval(-65.14583510270246, 0.5631666715107633)
dom = a
for m in [2, 5, 10]
    DAISY["bspline2"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline2"]["order $m"] = Dict()
    
    p = (-3 * u*u*u  + 3*u*u + 3*u + 1) / 6.0
    
    DAISY["bspline2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

DAISY["bspline3"] = BenchmarkGroup()
RELPREC["bspline3"] = Dict()
ref = Interval(0.004500002596107108, 15.187500453119108)
dom = a
for m in [2, 5, 10]
    DAISY["bspline3"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline3"]["order $m"] = Dict()
   
    p = -u*u*u / 6.0
    
    DAISY["bspline3"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline3"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline3"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline3"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end


# ==========
# himmilbeau
# ========== 
           

DAISY["himmilbeau"] = BenchmarkGroup()
RELPREC["himmilbeau"] = Dict()
ref = Interval(-2.58927, 110.461)
dom = a×b
for m in [2, 5, 10]
    DAISY["himmilbeau"]["order $m"] = BenchmarkGroup()
    RELPREC["himmilbeau"]["order $m"] = Dict()
    x1, x2 = set_variables(Float64,["x1","x2"],order = m)

    p = (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)
    
    DAISY["himmilbeau"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["himmilbeau"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["himmilbeau"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["himmilbeau"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end


# ==========
# kepler
# ==========

DAISY["kepler1"] = BenchmarkGroup()
RELPREC["kepler1"] = Dict()
ref = Interval(-68.62000208166202, 63.930001740164585)
dom = a×b×c×d×e×f
for m in [2, 5, 10]
    DAISY["kepler1"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler1"]["order $m"] = Dict()
    x1, x2, x3, x4, x5, x6 = set_variables(Float64,["x1","x2","x3",
                                                "x4","x5","x6"],order = m)
    p =  x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 +
      x1 * (-x1 + x2 + x3 - x4 + x5 + x6)

    DAISY["kepler1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

DAISY["kepler2"] = BenchmarkGroup()
RELPREC["kepler2"] = Dict()
dom = a×b×c×d
ref = Interval(-229.3700091323759, 89.34000554067282)
for m in [2, 5, 10]
    DAISY["kepler2"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler2"]["order $m"] = Dict()
    x1, x2, x3, x4 = set_variables(Float64,["x1","x2","x3","x4"],order = m)
    p = x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4) + 
        x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4

    DAISY["kepler2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

DAISY["kepler2"] = BenchmarkGroup()
RELPREC["kepler2"] = Dict()
dom = a×b×c×d×e×f
ref = Interval(-580.1000285855648, 89.34000554067282)
for m in [2, 5, 10]
    DAISY["kepler3"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler3"]["order $m"] = Dict()    
    
    x1, x2, x3, x4, x5, x6 = set_variables(Float64,["x1","x2","x3",
                                                "x4","x5","x6"],order = m)
    p =  x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6) + 
         x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6) -
         x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6
    DAISY["kepler3"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler3"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end    
    


# =========
# Rigidbody
# =========


DAISY["Rigidbody1"] = BenchmarkGroup()
RELPREC["Rigidbody1"] = Dict()
dom = a×b×c
ref = Interval(-21.27000038288167, -0.5399998451177108)
for m in [2, 5, 10]
    DAISY["Rigidbody1"]["order $m"] = BenchmarkGroup()
    RELPREC["Rigidbody1"]["order $m"] = Dict()

    x1, x2, x3 = set_variables(Float64,["x1","x2","x3"], order = m)

    p = -x1*x2 - 2*x2*x3 - x1 - x3

    DAISY["Rigidbody1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["Rigidbody1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["Rigidbody1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["Rigidbody1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end  


DAISY["Rigidbody2"] = BenchmarkGroup()
RELPREC["Rigidbody2"] = Dict()
dom = a×b×c
ref = Interval(68.81099823052828, 363.1424073663064)
for m in [2, 5, 10]
    DAISY["Rigidbody2"]["order $m"] = BenchmarkGroup()
    RELPREC["Rigidbody2"]["order $m"] = Dict()    
    
    x1, x2, x3 = set_variables(Float64,["x1","x2","x3"], order = m) 
        
    p = 2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2

    DAISY["Rigidbody2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["Rigidbody2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["Rigidbody2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["Rigidbody2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end   



