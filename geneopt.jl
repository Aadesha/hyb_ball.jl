using TaylorModels, IntervalArithmetic
using JuMP, Ipopt

p1 = Taylor1([0.0, 1.0, 1.0, 0.5], 4)
p2 = Taylor1([1.0, 3.0, -1.0], 4)
δ = 0.5; I = Interval(-δ, δ)
x₀ = Interval(0.0)
D = Interval(-3.0, 1.0)
TM = [TaylorModel1(p1, I, x₀, D), TaylorModel1(p2, I, x₀, D)]

poly1 = TM[1].pol

model = Model(with_optimizer(Ipopt.Optimizer))

@variable(model, x[1:length(TM)])

@variable(model, IN)

@variable(model, x1)

@variable(model, x2)

@NLparameter(model, l[1:length(TM)] == 10)

#@Nlexpression(model)

#my_expr = dot(l, x)

#register(my_expr)

#f1(x) = TM[1].pol(x)

#f2(x) = TM[2].pol(x)

@NLobjective(model, Max, l[1]*x[1]+l[2]*x[2])

@constraint(model, D.lo <= x1 <= D.hi)

@constraint(model, D.lo <= x2 <= D.hi)

@constraint(model, I.lo <= IN <= I.hi)

@constraint(model, x[1] == TM[1].pol(x1) + IN)

@constraint(model, x[2] == TM[2].pol(x2) + IN)

#@Nlexpression(model, Cons )

#@NLconstrain(model, )

Optimize!(model)
