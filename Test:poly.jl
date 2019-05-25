
using TaylorModels, BenchmarkTools, TaylorSeries
using IntervalArithmetic, IntervalOptimisation

# helper function that returns interval enclosures for the global minimum
# and global maximum of a univariate or multivariate polynomial over a given domain
function minmax(p, dom)
    global_min, _ = minimise(p, dom)
    minus_global_max, _ = minimise(-p, dom)
    global_max = -minus_global_max
    return global_min, global_max
end

#sin___________________________________________________#

Dx=Interval(-1.57079632679,1.57079632679)
dom=Dx

s="sin=x - (x*x*x)/6.0 + (x*x*x*x*x)/120.0 - (x*x*x*x*x*x*x)/5040.0;"
sin=Meta.parse(s).args[1].args[2]
x=Taylor1(7)
sin=eval(sin)

evaluate(sin,dom)
@benchmark evaluate($sin, $dom)
global_min, global_max = minmax(sin, dom)

#bspline0___________________________________________________#

Du=Interval(0.0,1.0)
dom=Du

s="bspline0=(1 - u) * (1 - u) * (1 - u) / 6.0;"
bspline0=Meta.parse(s).args[1].args[2]
u=Taylor1(3)
bspline0=eval(bspline0)

evaluate(bspline0, dom)
@benchmark evaluate($bspline0, $dom)
global_min, global_max = minmax(bspline0, dom)


s="bspline1=(3 * u*u*u - 6 * u*u + 4) / 6.0;"
bspline1=Meta.parse(s).args[1].args[2]
bspline1=eval(bspline1)

evaluate(bspline1, dom)
@benchmark evaluate($bspline1, $dom)
global_min, global_max = minmax(bspline1, dom)


s="bspline2=(-3 * u*u*u  + 3*u*u + 3*u + 1) / 6.0;"
bspline2=Meta.parse(s).args[1].args[2]
bspline2=eval(bspline2)

evaluate(bspline2, dom)
@benchmark evaluate($bspline2, $dom)
global_min, global_max = minmax(bspline2, dom)


s="bspline3=-u*u*u / 6.0"
bspline3=Meta.parse(s).args[1].args[2]
bspline3=eval(bspline3)

evaluate(bspline3, dom)
@benchmark evaluate($bspline3, $dom)
global_min, global_max = minmax(bspline3, dom)

#Doppler___________________________________________________________________#

Du=Interval(-100.0,100.0)
Dv=Interval(20.0,20000.0)
Dt=Interval(-30.0,50.0)
dom=Du×Dv×Dt

s="Doppler=(- ((331.4 + 0.6 * T)) *v) / (((331.4 + 0.6 * T) + u)*((331.4 + 0.6 * T) + u));"
Doppler=Meta.parse(s).args[1].args[2]
Doppler=eval(Doppler)
v,u,T=set_variables(Float64,["v","u","T"],order=4)

evaluate(Doppler, dom)
@benchmark evaluate($Doppler, $dom)
global_min, global_max = minmax(Doppler, dom)

#himmilbeau___________________________________________________________________#

Dx1=Interval(-5.0,5.0)
Dx2=Dx1
dom=Dx1×Dx2

s="himmilbeau=(x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7);"
himmilbeau=Meta.parse(s).args[1].args[2]
himmilbeau=eval(himmilbeau)
x1,x2=set_variables(Float64,["x1","x2"],order=5)

evaluate(himmilbeau, dom)
@benchmark evaluate($himmilbeau , $dom)
global_min, global_max = minmax(himmilbeau , dom)

#kepler________________________________________________________________________#

Dx1=Interval(4.0,6.36)
Dx2=Dx3=Dx4=Dx5=Dx6=Dx1
dom=Dx1×Dx2×Dx3×Dx4×Dx5×Dx6


x1,x2,x3,x4,x5,x6=set_variables(Float64,["x1","x2","x3","x4","x5","x6"],order=3)
s="kepler0=  x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 + x1 * (-x1 + x2 + x3 - x4 + x5 + x6);"
kepler0=Meta.parse(s).args[1].args[2]
kepler0=eval(kepler0)


evaluate(kepler0, dom)
@benchmark evaluate($kepler0, $dom)
global_min, global_max = minmax(kepler0, dom)

dom=Dx1×Dx2×Dx3×Dx4
x1,x2,x3,x4=set_variables(Float64,["x1","x2","x3","x4"],order=3)
s="kepler1=  x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4) + x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4;"
kepler1=Meta.parse(s).args[1].args[2]
kepler1=eval(kepler1)

evaluate(kepler1, dom)
@benchmark evaluate($kepler1, $dom)
global_min, global_max = minmax(kepler1, dom)

dom=Dx1×Dx2×Dx3×Dx4×Dx5×Dx6
x1,x2,x3,x4,x5,x6=set_variables(Float64,["x1","x2","x3","x4","x5","x6"],order=3)
s="kepler2=  x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6) + x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6) - x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6;"
kepler2=Meta.parse(s).args[1].args[2]
kepler2=eval(kepler2)

evaluate(kepler2, dom)
@benchmark evaluate($kepler2, $dom)
global_min, global_max = minmax(kepler2, dom)

#Rigidbody________________________________________________________#
Dx1=Interval(-15.0,15.0)
Dx2=Dx3=Dx1
dom=Dx1×Dx2×Dx3
x1,x2,x3=set_variables(Float64,["x1","x2","x3"],order=3)

s="Rigidbody1=-x1*x2 - 2*x2*x3 - x1 - x3;"
Rigidbody1=Meta.parse(s).args[1].args[2]
Rigidbody1=eval(Rigidbody1)

evaluate(Rigidbody1, dom)
@benchmark evaluate($Rigidbody1, $dom)
global_min, global_max = minmax(Rigidbody1, dom)

x1,x2,x3=set_variables(Float64,["x1","x2","x3"],order=4)
s="Rigidbody2=2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2;"
Rigidbody2=Meta.parse(s).args[1].args[2]
Rigidbody2=eval(Rigidbody2)

evaluate(Rigidbody2, dom)
@benchmark evaluate($Rigidbody2, $dom)
global_min, global_max = minmax(Rigidbody2, dom)

#turbine_________________________________________________________#

Dv=Interval(-4.5,.3)
Dw=Interval(-0.3,0.9)
Dr=Interval(3.8,7.8)
dom=Dv×Dw×Dr
v,w,r=set_variables(Float64,["v","w","r"],order=6)

s="turbine1= 3+ 2/(r*r) - (0.125*(3-2*v)*(w*w*r*r))/(1-v) - 4.5 ;"
turbine1=Meta.parse(s).args[1].args[2]
turbine1=eval(turbine1)

evaluate(turbine1, dom)
@benchmark evaluate($turbine1, $dom)
global_min, global_max = minmax(turbine1, dom)


s="turbine2=6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5;"
turbine2=Meta.parse(s).args[1].args[2]
turbine2=eval(turbine2)

evaluate(turbine2, dom)
@benchmark evaluate($turbine2, $dom)
global_min, global_max = minmax(turbine2, dom)


s="turbine3= 3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5;"
turbine3=Meta.parse(s).args[1].args[2]
turbine3=eval(turbine3)
evaluate(turbine3, dom)
@benchmark evaluate($turbine3, $dom)
global_min, global_max = minmax(turbine3, dom)
