
include("LeftMove.jl")
include("RightMove.jl")
include("UpMove.jl")
include("DownMove.jl")
include("renormalisation.jl")


function update!(C,T,tens,gt::Matrix{Symbol},chi::Int64)

  renormalisation!(C,T,gt)  

  LeftMove!(C,T,tens,gt,chi,4)
  renormalisation!(C,T,gt)  

  RightMove!(C,T,tens,gt,chi,2)
  renormalisation!(C,T,gt)    

  UpMove!(C,T,tens,gt,chi,1)
  renormalisation!(C,T,gt)  
  
  DownMove!(C,T,tens,gt,chi,3)
  renormalisation!(C,T,gt)  
 
  nothing

end