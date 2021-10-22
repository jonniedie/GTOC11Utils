# Need this because GTOC 11 uses a different definition of an AU than the standard
@unit AU "AU" AstronomicalUnit 149_597_870_691*m false

vec3(v) = ComponentArray(x=v[1], y=v[2], z=v[3])

state_vec(x) = ComponentArray(; r=vec3(x[1:3]), ṙ=vec3(x[4:6]))

quat_from_vec_angle(û, θ) = UnitQuaternion(cos(θ/2), sin(θ/2).*(û)...)

function delta_v_with_radius(v⃗ₘ, v⃗₁, v⃗₂, r)
    v⃗ₘ, v⃗₁, v⃗₂ = SVector{3}.((v⃗ₘ, v⃗₁, v⃗₂))
	l⃗₁ = v⃗ₘ - v⃗₁
	l⃗₂ = v⃗₂ - v⃗₁
	l₁ = norm(l⃗₁)
	l₂ = norm(l⃗₂)
	
	cosθ = (l⃗₁⋅l⃗₂)/(l₁*l₂)
	α = acos(l₂ / sqrt((l₂*cosθ-l₁)^2 + l₂^2))
	n̂ = normalize(l⃗₁×l⃗₂)
	
	q = quat_from_vec_angle(n̂, α)
	
	r⃗ = r*normalize(q*l⃗₁)
	
	Δv⃗₁ₘ = l⃗₁ - r⃗
	Δv⃗₁₂ = l⃗₂ - r⃗
	
	return Δv⃗₁ₘ, Δv⃗₁₂
end

function read_asteroids_file(fname; time_unit=yr, distance_unit=AU, angle_unit=°, mass_unit=kg, keep_units=true)
    data = readdlm(fname; header=true)

    # Make dataframe
    values = data[1]
    colnames = ["ID", "epoch", "a", "e", "i", "Ω", "ω", "M", "mass"]
    df = DataFrame(values, colnames)

    # Convert units
    df[!,1] =     df[:,1]     .|> Int
    df[!,2] =     df[:,2].*d  .|> time_unit
    df[!,3] =     df[:,3].*AU .|> distance_unit
    df[!,5:8] = df[:,5:8].*°  .|> angle_unit
    df[!,9] =     df[:,9].*kg .|> mass_unit

    return keep_units ? df : ustrip.(df)
end

