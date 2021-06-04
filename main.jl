#Set the directory
cd("C:\\Users\\engme\\Desktop\\Upload to GitHub")

using Distributed

#Relaxation tightening level
@everywhere RT=2 

#Number of pieces in piecewise relaxation
@everywhere K=1

#Number of processors
@everywhere np = 1

#Optimality gap tolerance
@everywhere  gap_tol = 0.0007

if nprocs() < np
   addprocs(np-nprocs());
end

@everywhere begin
  using CSV, JuMP, CPLEX, NLsolve, LightGraphs, Ipopt, Statistics, SharedArrays, DataFrames
end

@everywhere include("load_data.jl")
@everywhere include("LBmodel.jl");
include("UBmodel.jl");
include("BT.jl");

function Global(bounds, XL)
    UB = 1e20
    LB = 0
    #Number of iterations
    niter = 10
    for iter = 1:niter
        FBBT!(bounds)
	t0 = time();
        LB_trial,~, XL, Qi = get_LB(bounds,RT,K)
	println("LB time  ", time()-t0)
        LB = max(LB, LB_trial)
	println("LB:    ",LB,"   ",XL)
        println("Qi","   ",Qi)
    	UB = min(UB, get_UB_AC_OPF_model(bounds, XL, Qi));
        gap=(UB-LB)/UB;
	println("UB  ", UB)
	println("gap ",gap)
        if gap < gap_tol
           break;
	end
	if iter < niter
	   t0 = time();
    	   OBBT!(bounds, UB)
	   println("OBBT time  ", time()-t0)
	end
    end
end
Global(bounds,XL)
