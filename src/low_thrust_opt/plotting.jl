chase_plot(plt, sol; kwargs...) = chase_plot!(plt, sol; kwargs...)
chase_plot(sol; kwargs...) = chase_plot!(plot(), sol; kwargs...)

chase_plot!(sol; kwargs...) = chase_plot!(plot(), sol; kwargs...)
function chase_plot!(plt, sol; show_accel=true, kwargs...)
    plot_args = (;
        linewidth=3,
        aspect_ratio=:equal,
        autosize=false,
        xlims=(-3,3),
        ylims=(-3,3),
        zlims=(-3,3),
        kwargs...,
    )
    plot!(plt, sol, vars=(1,2,3); label="chaser", plot_args...)
    plot!(plt, sol, vars=(7,8,9); label="target", plot_args...)
    if show_accel
        t = range(extrema(sol.t)...; length=20)
        a = get_accelerations(sol, t) ./ 2
        r = reduce(hcat, [u.x.chaser.r for u in sol(t).u])'
        for (r,a) in zip(eachrow(r), eachrow(a))
            L = r + a
            plot!([r[1], L[1]], [r[2], L[2]], [r[3], L[3]]; color=1, label=false, plot_args...)
        end
    else
    end
    scatter!(plt, [0], [0], [0]; label="sun", color=:goldenrod, aspect_ratio=1)
end


acceleration_plot(plt, sol; kwargs...) = acceleration_plot!(plt, sol; kwargs...)
acceleration_plot(sol; kwargs...) = acceleration_plot!(plot(), sol; kwargs...)

acceleration_plot!(sol; kwargs...) = acceleration_plot!(plot(), sol; kwargs...)
function acceleration_plot!(plt, sol; kwargs...)
    t = range(extrema(sol.t)...; length=1000)
    a = get_accelerations(sol, t)
    plot!(plt, t, a; labels=["a_x" "a_y" "a_z"])
end