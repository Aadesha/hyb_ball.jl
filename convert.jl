
using Require,TaylorModels
using LazySets:Zonotope,Hyperrectangle,minkowski_sum

"""    convert(vTM::Vector{TaylorModelN{N, T, S}},
                    ::Zonotope)::Zonotope where{N, T, S}

Convert between TaylorModelN and in Zonotope.
### Input
- `vTM` -- target type
- `[P]` -- target set(Zonotope)
### Output
The Zonotope converted from TaylorModelN as the target set.
### example
julia> p=(x₁ + 3*x₁ + 3 -.1x₂ ,x₂ + 3*x₁ + 3 )
( 3.0 + 4.0 x₁ - 0.1 x₂ + 𝒪(‖x‖⁵),  3.0 + 3.0 x₁ + 1.0 x₂ + 𝒪(‖x‖⁵))
julia> TM=[TaylorModelN(p[1],I,x₀,D),TaylorModelN(p[2],I,x₀,D)]
2-element Array{TaylorModelN{2,Float64,Float64},1}:
  3.0 + 4.0 x₁ - 0.1 x₂ + [-0.5, 0.5]
  3.0 + 3.0 x₁ + 1.0 x₂ + [-0.5, 0.5]
julia> println(convert(TM,Zonotope))
  Zonotope{Float64}([8.85, 9.0],
    [1, 1]  =  10.0
    [2, 1]  =  7.5
    [1, 2]  =  -0.15
    [2, 2]  =  1.5
    [1, 3]  =  0.5
    [2, 4]  =  0.5)

"""
function convert(vTM::Vector{TaylorModelN{N, T, S}},
                    ::Type{Zonotope})::Zonotope where{N, T, S}

   TM_length=length(vTM)
   Q= [normalize_taylor(vTM[i].pol,vTM[i].dom,true) for i=1:TM_length]
   l_Vab=length(Q[1].coeffs[2].coeffs)
   Z=Zonotope([ (Q[i].coeffs[1].coeffs[1]) for i=1:TM_length],
               [ [(Q[i].coeffs[2].coeffs[j]) for i=1:TM_length] for j=1:l_Vab ])
  if TM_length==1
      I_Rem=TM[1].rem
  else
      for i=1:TM_length
          if i==1
              I_Rem=TM[1].rem
          else
              I_Rem=I_Rem×TM[i].rem
          end
       end
   end
   H=convert(Hyperrectangle,I_Rem)
   return minkowski_sum(Z, convert(Zonotope, H))
end
"""    convert(::Type{Zonotope}, vTM::Vector{TaylorModel1{T, S}}) where{T, S}

Convert between TaylorModel1 and in Zonotope.
### Input
- `vTM` -- Vector of TaylorModel1.
- `[P]` -- target set(Zonotope)
### Output
The Zonotope converted from TaylorModel1 as the target set.
### example
julia> p = Taylor1([2.0, 1.0], 2)
 2.0 + 1.0 t + 𝒪(t³)

julia> p1 = Taylor1([0.9, 3.0], 2)
 0.9 + 3.0 t + 𝒪(t³)

julia> TM = [TaylorModel1(p, I, x₀, D),TaylorModel1(p1, I, x₀, D)]
 2-element Array{TaylorModel1{Float64,Float64},1}:
   2.0 + 1.0 t + [-0.5, 0.5]
   0.9 + 3.0 t + [-0.5, 0.5]
julia> convert(Zonotope,TM)
   Zonotope{Float64}([1.0, -2.1],
     [1, 1]  =  2.0
     [2, 1]  =  6.0
     [1, 2]  =  0.5
     [2, 3]  =  0.5)

"""

function convert(::Type{Zonotope}, vTM::Vector{TaylorModel1{T, S}}) where{T, S}
   TM_length=length(vTM)
   Q= [normalize_taylor(vTM[i].pol,vTM[i].dom,true) for i=1:TM_length]
   #l_Vab=length(Q[1].coeffs[2].coeffs)
   Z=Zonotope([ (Q[i].coeffs[1]) for i=1:TM_length],
                [ vcat([Q[i].coeffs[2] for i=1:TM_length]) ])
   if TM_length==1
       I_Rem=TM[1].rem
   else
       for i=1:TM_length
          if i==1
              I_Rem=TM[1].rem
          else
              I_Rem=I_Rem×TM[i].rem
          end
       end
   end
   H=convert(Hyperrectangle,I_Rem)
   return minkowski_sum(Z, convert(Zonotope, H))
end




function overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope})where{T, S}
#(that the vector of TMs is linear, or explain in the docs until this assumption is lifted)
   # @assert TM.pol.order <=1
   l_vTM=length(vTM)
   Poly_L=[linear_polynomial(vTM[i].pol)+constant_term(vTM[i].pol) for i=1:l_vTM]
   Poly_NL=[(vTM[i].pol-Poly_L[i]) for i=1:l_vTM]
   Int_Poly_NL=[evaluate(Poly_NL[i],vTM[i].dom) for i=1:l_vTM]
   Rem=[(Int_Poly_NL[i]+vTM[i].rem) for i=1:l_vTM]
   TM=[ TaylorModel1(Poly_L[i],Rem[i],vTM[i].x0,vTM[i].dom) for i=1:l_vTM ]

    TM_length=length(TM)
    Q= [normalize_taylor(TM[i].pol,TM[i].dom,true) for i=1:TM_length]
    l_Vab=length(Q[1].coeffs[2].coeffs)
    Z=Zonotope([ (Q[i].coeffs[1]) for i=1:TM_length],
                [ vcat([Q[i].coeffs[2] for i=1:TM_length]) ])
    if TM_length==1
        I_Rem=TM[1].rem
    else
        for i=1:TM_length
            if i==1
                I_Rem=TM[1].rem
            else
                I_Rem=I_Rem×TM[i].rem
            end
        end
    end
    H=convert(Hyperrectangle,I_Rem)
    return minkowski_sum(Z, convert(Zonotope,H))
 end


 function overapproximate(::Type{Zonotope},
                       vTM::Vector{TaylorModelN{N, T, S}})where{N,T, S}

 #(that the vector of TMs is linear, or explain in the docs until this assumption is lifted)

    #@assert TM.pol.order <=1
    l_vTM=length(vTM)
    Poly_L=[linear_polynomial(vTM[i].pol)+constant_term(vTM[i].pol) for i=1:l_vTM]
    Poly_NL=[(vTM[i].pol-Poly_L[i]) for i=1:l_vTM]
    Int_Poly_NL=[evaluate(Poly_NL[i],vTM[i].dom) for i=1:l_vTM]
    Rem=[(Int_Poly_NL[i]+vTM[i].rem) for i=1:l_vTM]
    TM=[ TaylorModelN(Poly_L[i],Rem[i],vTM[i].x0,vTM[i].dom) for i=1:l_vTM ]

    TM_length=length(TM)
    Q= [normalize_taylor(TM[i].pol,TM[i].dom,true) for i=1:TM_length]
    l_Vab=length(Q[1].coeffs[2].coeffs)
    Z=Zonotope([ (Q[i].coeffs[1].coeffs[1]) for i=1:TM_length],
                [ [(Q[i].coeffs[2].coeffs[j]) for i=1:TM_length] for j=1:l_Vab ])
    if TM_length==1
       I_Rem=TM[1].rem
    else
       for i=1:TM_length
          if i==1
              I_Rem=TM[1].rem
          else
              I_Rem=I_Rem×TM[i].rem
          end
        end
    end
    H=convert(Hyperrectangle,I_Rem)
    return minkowski_sum(Z, convert(Zonotope, H))
 end
