
# include("find_value.jl")
# include("LeftMove.jl")
# include("RightMove.jl")
# include("UpMove.jl")
# include("DownMove.jl")
include("renormalisation.jl")

# module movesN
#   using ITensors
#   using ..ShastryGS
#   include("invert_diag_sqrt.jl")
#   include("find_value.jl")
#   include("LeftMove.jl")
#   include("RightMove.jl")
#   include("UpMove.jl")
#   include("DownMove.jl")
#   export LeftMove!, RightMove!, UpMove!, DownMove!
# end

# module movesN2
#   using ITensors
#   using ..ShastryGS
#   include("invert_diag_sqrt.jl")
#   include("find_value.jl")
#   include("movesN2/LeftMove.jl")
#   include("movesN2/RightMove.jl")
#   include("movesN2/UpMove.jl")
#   include("movesN2/DownMove.jl")
#   export LeftMove!, RightMove!, UpMove!, DownMove!
# end


# function load_moves(N)
#   if N == 2
#     module_name = :movesN2  # Example: decide dynamically which module to load
#     eval(Meta.parse("using .$(module_name)"))  # Dynamically import the module
#   else
#     module_name = :movesN  # Example: decide dynamically which module to load
#     eval(Meta.parse("using .$(module_name)"))  # Dynamically import the module
#   end
#   return nothing
# end

function update!(C,T,tens,gt::Matrix{Symbol},chi::Int64,N::Int64)


  # load_moves(N)

  renormalisation!(C,T,gt,N)  

  LeftMove!(C,T,tens,gt,chi,4,N)
  renormalisation!(C,T,gt,N)  

  RightMove!(C,T,tens,gt,chi,2,N)
  renormalisation!(C,T,gt,N)    

  UpMove!(C,T,tens,gt,chi,1,N)
  renormalisation!(C,T,gt,N)  
  
  DownMove!(C,T,tens,gt,chi,3,N)
  renormalisation!(C,T,gt,N)  
 
  nothing

end