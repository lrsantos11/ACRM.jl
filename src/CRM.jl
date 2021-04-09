"""
FindCircumcentermSet(X)

Finds the Circumcenter of linearly independent vectors ``x_0,x_1,…,x_m``, columns of matrix ``X``,
as described in [^Behling2018a] and [^Behling2018b].

[^Behling2018a]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018).
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018b]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018).
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)
"""

"""
    CRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function CRMiteration(xCRM::Vector, ReflectA::Function, ReflectB::Function)
    xCRM_RA = ReflectA(xCRM)
    xCRM_RBRA = ReflectB(xCRM_RA)
    if norm(xCRM_RA - xCRM)<ZERO_VAL
        xCRM = FindCircumcentermSet([xCRM, xCRM_RBRA])
    elseif norm(xCRM_RBRA - xCRM_RA)<ZERO_VAL
        xCRM =FindCircumcentermSet([xCRM,  xCRM_RA])
    else
        xCRM = FindCircumcentermSet([xCRM, xCRM_RA, xCRM_RBRA])
    end
    return xCRM  
end 

"""
    CRM()

Cirumcentered-Reflection method (CRM) as proposed in [^Behling2018a], [^Behling2018b] and [^Behling2021]

[^Behling2018a]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018).
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018b]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018).
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)
[^Behling2021]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the Circumcentered-Reflection Method for the Convex Feasibility Problem.  Numer. Algorithms. 86, 1475--1494 (2021).
[doi:10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6)
"""
function CRM(x₀::Vector,ProjectA::Function, ProjectB::Function; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolCRM = 1.
    xCRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    printoOnFile(filedir,hcat(k, tolCRM, xCRM'),deletefile=true)
    while tolCRM > EPSVAL && k < itmax
        xCRMOld = copy(xCRM)
        print_intermediate ?  printoOnFile(filedir,hcat(nothing,nothing,ProjectA(xCRM)')) : nothing
        xCRM  = CRMiteration(xCRM, ReflectA, ReflectB)
        tolCRM = gap_distance ? norm(ProjectA(xCRM)-xCRM) : Tolerance(xCRM,xCRMOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolCRM, xCRM'))

    end
    return Results(iter_total= k,final_tol=tolCRM,xApprox=xCRM,method=:CRM)
end
"""
    CRMprod(x₀, Projections)

    Cirumcentered-Reflection method on Pierra's product space reformulation  as proposed in [^Behling2021]

[^Behling2021]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the Circumcentered-Reflection Method for the Convex Feasibility Problem.  Numer. Algorithms. 86, 1475--1494 (2021).
[doi:10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6)
"""
function CRMprod(x₀::Vector{Float64},Projections::Vector; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolCRMprod = 1.
    num_sets = length(Projections)
    xCRMprod = Vector[]
    for i = 1:num_sets
        push!(xCRMprod,x₀)
    end
    ProjectAprod(x) = ProjectProdSpace(x,Projections)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    ReflectA(x) = Reflection(x,ProjectAprod)
    ReflectB(x) = Reflection(x,ProjectBprod)
    printoOnFile(filedir,hcat(k, tolCRMprod, xCRMprod[1]'),deletefile=true)
    while tolCRMprod > EPSVAL && k < itmax
        xCRMprodOld = copy(xCRMprod)
        print_intermediate ?  printoOnFile(filedir,hcat(nothing,nothing,(ProjectA(xCRMprod))[1]')) : nothing
        xCRMprod  = CRMiteration(xCRMprod, ReflectA, ReflectB)
        tolCRMprod = gap_distance ? norm(ProjectAprod(xCRMprod)-xCRMprod) : Tolerance(xCRMprod,xCRMprodOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolCRMprod, xCRMprod[1]'))
    end
    return Results(iter_total= k,
                  final_tol=tolCRMprod,xApprox=xCRMprod[1],method=:CRMprod)
end    

"""
    CRMprod(x₀, ProjectA, ProjectB)
"""
CRMprod(x₀::Vector{Float64},ProjectA::Function, ProjectB::Function;kwargs...) = CRMprod(x₀,[ProjectA,ProjectB],kwargs...) 