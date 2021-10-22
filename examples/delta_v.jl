using LinearAlgebra: norm, normalize, ⋅, ×
using GeometryBasics: Point, Sphere
using GLMakie: mesh, quiver!, annotations!
using GTOC11Utils: delta_v_with_radius


##
vm = [8, 2, 0]
v1 = [5, 0, 5]
v2 = [-1, 3, 3]
r = 1


##
Δvm, Δv2 = delta_v_with_radius(vm, v1, v2, r)


##
vec_wrap(x) = [x]

plt = mesh(Sphere(Point3(v1), r); color=:violet)
for v in [vm, v1, v2]
    quiver!([0], [0], [0], vec_wrap.(v)...)
end
quiver!(vec_wrap.(v1)..., vec_wrap.(vm-v1)...; color=:blue)
quiver!(vec_wrap.(v1)..., vec_wrap.(v2-v1)...; color=:blue)
quiver!(vec_wrap.(vm-Δvm)..., vec_wrap.(Δvm)...; color=:green)
quiver!(vec_wrap.(v2-Δv2)..., vec_wrap.(Δv2)...; color=:green)
annotations!(["vm", "v1", "v2"], Point3.([vm.+0.2, v1, v2.+0.2]))
display(plt)