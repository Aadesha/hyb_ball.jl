using TaylorModels, Plots, TaylorSeries, IntervalArithmetic
using LazySets:Zonotope
function TM_TO_ZONO_1D(a::Array{TaylorModelN{S,T,T},N},b::IntervalBox{S,T}) where s<:Int ::Zonotope
   l_TM=length(a)
   Q= [normalize_taylor(TM2[i].pol,D,true) for i=1:l_TM]
   C_Z=[ (Q[i].coeffs[1].coeffs[1]) for i=1:l_TM ]
   G_Z= [ [(Q[i].coeffs[2].coeffs)] for i=1:l_TM ]
   Z=Zonotope(C_Z,G_Z)
end
