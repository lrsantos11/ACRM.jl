__precompile__()
# module ACRM
export CRM, CRMprod, MAP, MAPprod
include("crm_utils.jl")
global const ZERO_VAL = 1e-15
include("CRM.jl")
include("MAP.jl")
# end
