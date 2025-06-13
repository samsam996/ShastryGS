


include("renormalisation.jl")

function update!(C,T,tens,gt::Matrix{Symbol},chi::Int64,N::Int64)

  """
    Perform one single iteration for the CTMRG; left-right-up and down moves.

    ARGS:   
        C, and T : the environement
        tens: the local tensors
        gt: the grid of the tensors
        chi: the bond dimension of the environement
        N: the size of the unit cell.
  """

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