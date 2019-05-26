function Linearizeonp(a::TaylorModelN)::TaylorModelN
   l_TM=length(a.pol)
   Poly_L=TaylorN(a.pol.coeffs[1:2])
   Poly_NL=TaylorN(a.pol.coeffs[3:l_TM])
   Int_Poly_NL=evaluate(Poly_NL,a.dom)
   Rem=Int_Poly_NL+a.rem
   TM=TaylorModelN(Poly_L,Rem,a.x0,a.dom)
end
