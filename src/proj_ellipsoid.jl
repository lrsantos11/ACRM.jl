

####################################
"""
        ProjEpiEllipsoid()
Project ``x₀ ∈ R^{n}`` onto the elllipsoid ``\mathcal{E} =: {x | q(x) = x^TAx + 2b^Tx ≤ α``.
using Algorithm 4 as described in [^Dai2006].


[^Dai2006]: Dai, Yu-Hong: Fast algorithms for projection onto an Ellipsoid. SIAM J. Optim. 16, 986–1006 (2006).
[doi:10.1137/040613305](https://doi.org/10.1137/040613305)
"""

 function ProjEpiEllipsoid(x₀,A,b,α,a; ϵ=1e-6, max_itr=120, 
                m₁ = 1, m₂ = 1, c₁ = 0.1, c₂ = 0.8 )
    xₖ=x₀
    gₐ=A*a + b
    Quadratic(x) = dot(x,A*x)+2*dot(b,x)
    qa_α = Quadratic(a)-α
    status = :MaxIter
    k=0
    while k<= max_itr 
        uₖ =A*xₖ+b
        vₖ =a-xₖ 
        uₖdotvₖ=dot(uₖ,vₖ)
        vₖdotvₖ=dot(vₖ,vₖ)
        normuₖ=norm(uₖ)
        normvₖ=sqrt(vₖdotvₖ)
        Mₖ=uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
        E=1-Mₖ 
        if E < ϵ
            status = :Optimal
            return status, xₖ,k,E         
        end
        zₖ =uₖ-((uₖdotvₖ)/(vₖdotvₖ))*vₖ
        normzₖ=norm(zₖ)
        zₖdotzₖ=dot(zₖ,zₖ)
        uₖdotzₖ=dot(uₖ,zₖ)
        Avₖ=A*vₖ 
        Azₖ=A*zₖ 
        vₖdotAvₖ=dot(vₖ,Avₖ)
        zₖdotAzₖ=dot(zₖ,Azₖ)
        vₖdotAzₖ=dot(vₖ,Azₖ)
        zₖdotAvₖ=dot(zₖ,Avₖ)
        
        #Alg1================================================
        γₖ¹ = compute_Alg1_γk(vₖdotAvₖ,vₖdotvₖ,vₖdotAzₖ,zₖdotAzₖ,zₖdotzₖ)   
        #Alg2================================================
        γₖ² = compute_Alg2_γk(normvₖ,normzₖ,uₖdotvₖ,uₖdotzₖ,vₖdotvₖ,vₖdotAvₖ,vₖdotAzₖ,zₖdotAvₖ,zₖdotAzₖ)
        #Safeguard for Alg2=================================
        if γₖ² < 0. || abs(γₖ²) <= ϵ
            γₖ² = γₖ¹
        end
        #Alg4================================================
        if mod(k,m₁+m₂) < m₁
            γₖ = γₖ²
        else
            γₖ = c₁* γₖ¹  + c₂* γₖ²
        end

        ωₖ=-(γₖ*uₖ)-vₖ
        ωₖdotAωₖ = dot(ωₖ,A*ωₖ)
        ηₖ = compute_ηk(gₐ,qa_α,ωₖ,ωₖdotAωₖ)
        xₖ=a.+(ηₖ*ωₖ)
        k+=1
    end
  
    return status, xₖ,k,E
end

function  selectRoot(β,h₁,h₂,λ₁,λ₂)   
    p = [0.,0]
    #solve Quadratic in lemma 1 to find p₁ then p₂ with following formula
    #P=((λ₁*(λ₁-λ₂)^2))*(p₁)^4+(2*λ₁*λ₂*h₁*(λ₁-λ₂))*(p₁)^3+((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))*(p₁^2)-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)*(p₁)-(β*λ₂^2*h₁^2)
    a4=((λ₁*(λ₁-λ₂)^2))
    a3=(2*λ₁*λ₂*h₁*(λ₁-λ₂))
    a2=((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))
    a1=-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)
    a0=-(β*λ₂^2*h₁^2)
    S=roots([a0, a1, a2, a3, a4])
    #choose just Real roots =============================
    indexreal=findall(x->abs.(x)<1e-12,imag.(S))
    Sol=real.(S[indexreal])
 
    #choose the right root ============================================
    !isempty(Sol) ? nothing :  @error "No solution for the quartic equation"
    for  p₁ in Sol
        sol1 = p₁
        μ₁=(h₁-p₁)/(λ₁*p₁)
        if μ₁<=0.
            continue
        else
            p₂=(λ₁*h₂)/((λ₁-λ₂)*p₁+λ₂*h₁)*p₁
            μ₂=(h₂-p₂)/(λ₂*p₂)
            if  μ₂<=0. && !(μ₂ ≈ μ₁) 
                continue
            end
        end
        sol2 = p₂
        p = [sol1,sol2]
    end
    return   p
end

function compute_Alg1_γk(vₖdotAvₖ,vₖdotvₖ,vₖdotAzₖ,zₖdotAzₖ,zₖdotzₖ)
    M1=(vₖdotAvₖ)/(vₖdotvₖ)
    M2=(zₖdotAzₖ)/(zₖdotzₖ)
    M3=((M1-M2)^2 )
    M4=4*(vₖdotAzₖ)^2/(vₖdotvₖ*zₖdotzₖ)
    M5=sqrt(M3+M4)
    ρₖ=0.5*(M1+M2+M5)
    γₖ=1/ρₖ
    return γₖ
end

function compute_ηk(gₐ,qa_α,wₖ,wₖdotAwₖ)
    N1=dot(gₐ,wₖ)/(wₖdotAwₖ)
    N2=(N1).^2
    N3=(qa_α)./(wₖdotAwₖ)
    ηₖ=-N1-sqrt(N2-N3)
    return ηₖ
end

function compute_Alg2_γk(normvₖ,normzₖ,uₖdotvₖ,uₖdotzₖ,vₖdotvₖ,vₖdotAvₖ,vₖdotAzₖ,zₖdotAvₖ,zₖdotAzₖ)
       aₗ=[normvₖ,0]
       Aₖ=[(vₖdotAvₖ)/(normvₖ^2)  (vₖdotAzₖ)/(normvₖ*normzₖ) ; 
             (zₖdotAvₖ)/(normvₖ*normzₖ)  (zₖdotAzₖ)/(normzₖ^2)]
       bₖ=[(uₖdotvₖ)/(normvₖ);(uₖdotzₖ)/(normzₖ)]
       D, Q = eigen(Aₖ)
       Qbₖ=Q*bₖ
       D⁻¹Qbₖ=D.\Qbₖ
       β=dot(Qbₖ,D⁻¹Qbₖ)
       aₕ=Q*aₗ+D⁻¹Qbₖ
       #Solve Quadratic=========================================
       aₕ⁽ᵖ⁾  = selectRoot(β,aₕ[1],aₕ[2],D[1],D[2])
       aₗ⁽ᵖ⁾=Q'*(aₕ⁽ᵖ⁾- D⁻¹Qbₖ)
       ηₖ² = 1 - (aₗ⁽ᵖ⁾[1]/normvₖ- (uₖdotvₖ/vₖdotvₖ)*aₗ⁽ᵖ⁾[2]/normzₖ)
       γₖ²=(- aₗ⁽ᵖ⁾[2]/normzₖ)/(ηₖ²)
    return γₖ²
end
