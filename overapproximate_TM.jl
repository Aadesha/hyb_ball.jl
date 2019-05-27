"""
         overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope})where{T, S}

Convert between TaylorModel1 and in Zonotope.

### Input


- `vTM` -- Vector of TaylorModel1
- `[P]` -- target type

### Output

The Zonotope converted from vector of TaylorModel1 as the target set.

### example

julia> using TaylorModels,IntervalArithmetic,Revise

julia> using LazySets:Hyperrectangle,Zonotope,convert

julia> m = 2
2

julia> δ = 0.5; I = Interval(-δ, δ)
[-0.5, 0.5]

julia> x₀ = Interval(0.0)
[0, 0]

julia> D = Interval(-3.0, 1.0)
[-3, 1]

julia> p = Taylor1([2.0, 1.0], 2)
 2.0 + 1.0 t + 𝒪(t³)

julia> p1 = Taylor1([0.9, 3.0], 2)
 0.9 + 3.0 t + 𝒪(t³)

julia> TM = [TaylorModel1(p, I, x₀, D),TaylorModel1(p1, I, x₀, D)]
2-element Array{TaylorModel1{Float64,Float64},1}:
 2.0 + 1.0 t + [-0.5, 0.5]
 0.9 + 3.0 t + [-0.5, 0.5]

julia> overapproximate(Zonotope,TM)
Zonotope{Float64}([1.0, -2.1],
   [1, 1]  =  2.0
   [2, 1]  =  6.0
   [1, 2]  =  0.5
   [2, 3]  =  0.5)

"""
function overapproximate(vTM::Vector{TaylorModel1{T, S}},
                            ::Type{Zonotope})where{T, S}
#(that the vector of TMs is linear, or explain in the docs until this assumption is lifted)
   # @assert TM.pol.order <=1
   l_vTM = length(vTM)
   Poly_L = [ linear_polynomial(vTM[i].pol) + constant_term(vTM[i].pol) for i=1:l_vTM ]
   Poly_NL = [ (vTM[i].pol - Poly_L[i]) for i = 1:l_vTM ]
   Int_Poly_NL = [ evaluate(Poly_NL[i], vTM[i].dom) for i=1:l_vTM ]
   Rem = [ (Int_Poly_NL[i] + vTM[i].rem) for i = 1:l_vTM ]
   TM = [ TaylorModel1(Poly_L[i], Rem[i], vTM[i].x0, vTM[i].dom) for i=1:l_vTM ]
   TM_length = length(TM)
   Q = [ normalize_taylor(TM[i].pol, TM[i].dom, true) for i = 1:TM_length ]
   Z = Zonotope([ (Q[i].coeffs[1]) for i = 1:TM_length ],
                [ vcat([Q[i].coeffs[2] for i = 1:TM_length]) ])
   if TM_length == 1
        I_Rem = TM[1].rem
   else
       for i = 1:TM_length
          if i == 1
               I_Rem = TM[1].rem
          else
               I_Rem = I_Rem×TM[i].rem
          end
       end
    end
    H = convert(Hyperrectangle, I_Rem)
    return minkowski_sum(Z, convert(Zonotope, H))
 end

 """
          overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                             ::Type{Zonotope})where{N,T, S}


 Convert between TaylorModelN and in Zonotope.

 ### Input

 - `vTM` -- Vector of TaylorModelN
 - `[P]` -- target type

 ### Output

 The Zonotope converted from vector of TaylorModelN as the target set.

 ### example

julia> m = 4
4

julia> x₁, x₂ = set_variables(Float64, ["x₁", "x₂"], order=2*m)
2-element Array{TaylorN{Float64},1}:
  1.0 x₁ + 𝒪(‖x‖⁹)
  1.0 x₂ + 𝒪(‖x‖⁹)

julia> x₀ = Interval(0.0, 0.0) × Interval(0.0, 0.0)
[0, 0] × [0, 0]

julia> Dx₁ = Interval(0.0, 3.0)
[0, 3]

julia> Dx₂ = Interval(-1.0, 1.0)
[-1, 1]

julia> D = Dx₁ × Dx₂
[0, 3] × [-1, 1]

julia> δ = 0.5; I = Interval(-δ, δ)
[-0.5, 0.5]

julia> p = 1 + x₁^2 - x₂
 1.0 - 1.0 x₂ + 1.0 x₁² + 𝒪(‖x‖⁹)

julia> p1 = x₂^3 + 3x₁^4 + x₁ + 1
 1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + 𝒪(‖x‖⁹)

julia> TM = [ TaylorModelN(p, I, x₀, D), TaylorModelN(p1, I, x₀, D) ]
2-element Array{TaylorModelN{2,Float64,Float64},1}:
             1.0 - 1.0 x₂ + 1.0 x₁² + [-0.5, 0.5]
   1.0 + 1.0 x₁ + 1.0 x₂³ + 3.0 x₁⁴ + [-0.5, 0.5]

julia> overapproximate(TM, Zonotope)
Zonotope{Float64}([5.5, 124.0],
   [2, 1]  =  1.5
   [1, 2]  =  -1.0
   [1, 3]  =  5.0
   [2, 4]  =  123.0)

"""
 function overapproximate(vTM::Vector{TaylorModelN{N, T, S}},
                             ::Type{Zonotope})where{N,T, S}

 #(that the vector of TMs is linear, or explain in the docs until this assumption is lifted)
    #@assert TM.pol.order <=1
    l_vTM = length(vTM)
    Poly_L = [ linear_polynomial(vTM[i].pol) + constant_term(vTM[i].pol) for i = 1:l_vTM ]
    Poly_NL = [ (vTM[i].pol - Poly_L[i]) for i = 1:l_vTM ]
    Int_Poly_NL = [ evaluate(Poly_NL[i], vTM[i].dom) for i = 1:l_vTM ]
    Rem = [ (Int_Poly_NL[i] + vTM[i].rem) for i = 1:l_vTM ]
    TM = [ TaylorModelN(Poly_L[i], Rem[i], vTM[i].x0, vTM[i].dom) for i = 1:l_vTM ]
    TM_length = length(TM)
    Q = [ normalize_taylor(TM[i].pol, TM[i].dom, true) for i = 1:TM_length ]
    l_Vab = length(Q[1].coeffs[2].coeffs)
    Z = Zonotope([ (Q[i].coeffs[1].coeffs[1]) for i = 1:TM_length ],
                [ [ (Q[i].coeffs[2].coeffs[j]) for i = 1:TM_length ] for j = 1:l_Vab ])
    if TM_length == 1
       I_Rem = TM[1].rem
    else
       for i = 1:TM_length
          if i == 1
              I_Rem = TM[1].rem
          else
              I_Rem = I_Rem×TM[i].rem
          end
        end
    end
    H = convert(Hyperrectangle, I_Rem)
    return minkowski_sum(Z, convert(Zonotope, H))
 end
