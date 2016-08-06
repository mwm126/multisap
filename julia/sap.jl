module sap

type Endpoint
    #  owner of the endpoint
    #  if Min endpoint, id = -box_id
    #  if Max endpoint, id = box_id
    box_id::Int

    #  actual value of the endpoint; single precision to save space
    value::Float32
end

type AABB
    maxendpoint::Array{Float32}
    minendpoint::Array{Float32}
end

type Box{T}
    maxendpoint_id::Array{Int}
    minendpoint_id::Array{Int}
    particle_id::T
end

type Sap
    # ! id to identify within a multisap
    id::Int

    # ! list of endpoints of boxes
    x_endpoints::Array{Endpoint}
    y_endpoints::Array{Endpoint}
    z_endpoints::Array{Endpoint}

    ht::Dict

    # ! list of boxes
    boxes::Array{Box}
    boxes_len::Integer

    # ! list of indices in boxes that have value 0 (deleted boxes)
    deleted_boxes::Array{Integer}

    # ! list of pairs of particle ids that have intersecting boxes
    # type(hashtable_t) :: hashtable
end

function add_box(sap, aabb, id)

    push!(sap.x_endpoints, Endpoint( -id, aabb.minendpoint(1)))
    push!(sap.y_endpoints, Endpoint( -id, aabb.minendpoint(2)))
    push!(sap.z_endpoints, Endpoint( -id, aabb.minendpoint(3)))
    MIN = length(sap.x_endpoints)

    push!(sap.x_endpoints, Endpoint( id, aabb.maxendpoint(1)))
    push!(sap.y_endpoints, Endpoint( id, aabb.maxendpoint(2)))
    push!(sap.z_endpoints, Endpoint( id, aabb.maxendpoint(3)))
    MAX = length(sap.x_endpoints)

    push!(sap.boxes, Box( (MAX, MAX, MAX), (MIN, MIN, MIN), id))
end

function update_box(sap, id, aabb)

    sap.x_endpoints(sap.boxes[id].maxendpoint_id[1]).value = aabb.maxendpoint[1]
    sap.y_endpoints(sap.boxes[id].maxendpoint_id[2]).value = aabb.maxendpoint[2]
    sap.z_endpoints(sap.boxes[id].maxendpoint_id[3]).value = aabb.maxendpoint[3]
    sap.x_endpoints(sap.boxes[id].minendpoint_id[1]).value = aabb.minendpoint[1]
    sap.y_endpoints(sap.boxes[id].minendpoint_id[2]).value = aabb.minendpoint[2]
    sap.z_endpoints(sap.boxes[id].minendpoint_id[3]).value = aabb.minendpoint[3]
end

function get_pairs(bodies)
    # first step, update all these bodies
end

function update_body(sap::Sap, body::Body)
end

# Usage:
# call add_box at beginning for each body
#
# Each timestep:
#
