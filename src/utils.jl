# Need this because GTOC 11 uses a different definition of an AU than the standard
@unit AU "AU" AstronomicalUnit 149_597_870_691*m false

vec3(v) = ComponentArray(x=v[1], y=v[2], z=v[3])

state_vec(x) = ComponentArray(; r=vec3(x[1:3]), ṙ=vec3(x[4:6]))