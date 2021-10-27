maybe_reverse(x, cond) = cond ≥ zero(cond) ? x : 2π-x

strip_ephemeris(x) = SVector(
    ustrip(yr, x.epoch),
    ustrip(AU, x.a),
    x.e,
    ustrip(rad, x.i),
    ustrip(rad, x.Ω),
    ustrip(rad, x.ω),
    ustrip(rad, x.M),
)

# propagate(t::Quantity, x::DataFrameRow; kwargs...) = propagate(t, strip_to_rad(@view(x[2:end-1])); kwargs...)
function propagate(t::Quantity, x; mu=μ)
    x = SVector(
        ustrip(yr, x.epoch),
        ustrip(AU, x.a),
        x.e,
        ustrip(rad, x.i),
        ustrip(rad, x.Ω),
        ustrip(rad, x.ω),
        ustrip(rad, x.M),
    )
    return propagate(ustrip(yr, t), x; mu=ustrip(AU^3/yr^2, mu))
end		

function propagate(t, eph; mu=μ)
    x = [
        eph.t0,
        eph.a,
        eph.e,
        eph.i,
        eph.Ω,
        eph.ω,
        eph.M0,
    ]
    return propagate(t, x; mu=ustrip(AU^3/yr^2, mu))

end

function propagate(t,x::AbstractVector; mu=ustrip(μ))
    t0,a,e,i,W,w,M0 = x
    n = sqrt((mu)/a^3)
    M = n*(t-t0)+M0
 
    #Using Newton's method to solve Keplers Equation for E Eccentric Anomaly  
    f(x) = x-e*sin(x)-M
    fp(x) = 1-e*cos(x)
    E = 0.1 #first guess
    epsilon = 1f-6
    F = 1 # just to get started
    ctr = 0    
    while abs(F) > epsilon
        ctr = ctr+1
        if ctr >= 1000            
            break            
        end
        F = f(E)
        E = E - f(E)/fp(E)
    end    
 
    #True Anomaly
    f = 2*atan(tan(E/2)/sqrt((1-e)/(1+e)))
   
    #Flight Path Angle
    gamma = atan(e*sin(f)/(1+e*cos(f)))
 
    #Norm of the radius vector
    r = a*(1-e^2)/(1+e*cos(f))
 
    #Norm of the velocity vector
    v = sqrt(2*mu/r-mu/a)
 
    x = r*(cos(f+w)*cos(W)-sin(f+w)*cos(i)*sin(W))
    y = r*(cos(f+w)*sin(W)+sin(f+w)*cos(i)*cos(W))
    z = r*(sin(f+w)*sin(i))
 
    vx = v*(-sin(f+w-gamma)*cos(W)-cos(f+w-gamma)*cos(i)*sin(W))
    vy = v*(-sin(f+w-gamma)*sin(W)+cos(f+w-gamma)*cos(i)*cos(W))
    vz = v*(cos(f+w-gamma)*sin(i))
 
    # rf = SVector(x,y,z) #.|> ustrip
    # vf = SVector(vx,vy,vz) #.|> ustrip
   
    return SVector(x,y,z,vx,vy,vz)
end

# THIS DOESN'T WORK
function kepler_to_cartesian(df::DataFrame; out_fun=ustrip, tol=1e-12)
    named_tuples = NamedTuple.(eachrow(df[:,2:end-1]))
    state_vectors = [out_fun.(kepler_to_cartesian(; nt...)) for nt in named_tuples]
    return collect(reduce(hcat, state_vectors)')
end
function kepler_to_cartesian(; a, e, i, Ω, ω, M, epoch=0yr, t0=0yr, μ=μ, tol=1e-12)
	Δt = epoch-t0
	M = M + Δt*sqrt(μ/a^3)
	
	Eⱼ = Inf
	Eⱼ₊₁ = M
	while abs(Eⱼ - Eⱼ₊₁) > tol
		Eⱼ = Eⱼ₊₁
		sE, cE = sincos(Eⱼ)
		Eⱼ₊₁ = Eⱼ - (Eⱼ - e*sE - M)/(1 - e*cE)
	end
	E = Eⱼ₊₁
	sE, cE = sincos(E)
	sE2, cE2 = sincos(E/2)
	
	# ν = 2atan(sqrt(1+e)*sE2, sqrt(1-e)cE2)
	# ν = 2*atan(tan(E/2), sqrt((1-e)/(1+e)))
    ν = 2*atan(tan(E/2)/sqrt((1-e)/(1+e)))
	sν, cν = sincos(ν)
	
	rc = a*(1-e*cE)
	
	o⃗ = rc*SVector(cν, sν, 0)
	o⃗̇ = sqrt(μ*a)/rc*SVector(-sE, sqrt(1-e^2)*cE, 0)
	
	R = RotZXZ(-Ω, -i, -ω)
	r⃗ = R * o⃗
	r⃗̇ = R * o⃗̇
	
	return [r⃗; r⃗̇]
end

function cartesian_to_kepler(r⃗r⃗̇, μ=ustrip(μ))
    r⃗, r⃗̇ = @views r⃗r⃗̇[1:3], r⃗r⃗̇[4:6]
	r⃗, r⃗̇ = SVector{3}(r⃗), SVector{3}(r⃗̇)
	r = norm(r⃗)
	ṙ = norm(r⃗̇)
	
	h⃗ = r⃗×r⃗̇
	e⃗ = (r⃗̇×h⃗)/μ - r⃗/r
	n⃗ = SVector(0,0,1) × h⃗
	h = norm(h⃗)
	e = norm(e⃗)
	n = norm(n⃗)
	
	ν = maybe_reverse(acos(e⃗⋅r⃗/(e*r)), r⃗⋅r⃗̇)
	i = acos(h⃗[3]/h)
	E = 2atan(tan(ν/2), sqrt((1+e)/(1-e)))
	Ω = maybe_reverse(acos(n⃗[1]/n), n⃗[2])
	ω = maybe_reverse(acos(n⃗⋅e⃗/(n*e)), e⃗[3])
	M = E - e*sin(E)
	a = 1 / (2/r - ṙ^2/μ)
	
	return (; a, e, i, Ω, ω, M)
end

