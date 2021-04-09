__precompile__()
# module ACRM
# export CRM, CRMprod, MAP, MAPprod
include("methods_utils.jl")
include("proj_ellipsoid.jl")
global const ZERO_VAL = 1e-15
include("CRM.jl")
include("MAP.jl")
# end
