# Need this because GTOC 11 uses a different definition of an AU than the standard
@unit AU "AU" AstronomicalUnit 149_597_870_691*m false

vec3(v) = ComponentArray(x=v[1], y=v[2], z=v[3])

state_vec(x) = ComponentArray(; r=vec3(x[1:3]), ṙ=vec3(x[4:6]))

propagate(t, x::DataFrameRow; kwargs...) = propagate(t, NamedTuple(@view(x[2:end-1])); kwargs...)
function propagate(t,x; mu=μ)
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
 
    xf = SVector(x,y,z,vx,vy,vz)
   
    return xf
end

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

function cartesian_to_kepler(r⃗r⃗̇, μ=μ)
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

maybe_reverse(x, cond) = cond ≥ zero(cond) ? x : 2π-x

function read_asteroids_file(fname)
    data = readdlm(fname; header=true)
    values = data[1]
    colnames = ["ID", "epoch", "a", "e", "i", "Ω", "ω", "M", "mass"]
    df = DataFrame(values, colnames)
    df[!,1] = df[:,1] .|> Int
    df[!,2] = df[:,2] .* d .|> yr
    df[!,3] = df[:,3] .* AU
    df[!,5:8] = df[:,5:8] .* °
    df[!,9] = df[:,9] .* kg
    return df
end

