@kwdef @concrete struct OneVehicleSimState
    x
    λ
end

const sim_proto_one_vehicle = ComponentArray(
    x = state_vec(zeros(6)),
    λ = state_vec(zeros(6)),
)

ComponentArrays.ComponentArray(x::OneVehicleSimState) = ComponentArray([x.x; x.λ], getaxes(sim_proto_one_vehicle))


@kwdef @concrete struct TwoVehicleSimState
    chaser
    target
    λ
end

const sim_proto_two_vehicles = ComponentArray(
    x = (
        chaser = state_vec(zeros(6)),
        target = state_vec(zeros(6)),
    ),
    λ = state_vec(zeros(6)),
)

ComponentArrays.ComponentArray(x::TwoVehicleSimState) = ComponentArray([x.chaser; x.target; x.λ], getaxes(sim_proto_two_vehicles))



@kwdef @concrete struct OptInput
    λ
    t
end

const opt_proto = ComponentArray(
    λ = state_vec(zeros(6)),
    t = 0,
)

ComponentArrays.ComponentArray(x::OptInput) = ComponentArray([x.λ; x.t], getaxes(opt_proto))
