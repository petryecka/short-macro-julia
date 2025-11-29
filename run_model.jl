using Dynare;
using ForwardDiff;
using LinearAlgebra;
using Optim;
using BenchmarkTools;
using PATHSolver;


context = @dynare "HS_nk_loglin.mod"

#acces objects within context
#res = context.results.model_results[1]
#propertynames(res)
#decision rules: oo_.ghx in
#res.linearrationalexpectations.g1_1
#oo_.ghu in
#res.linearrationalexpectations.g1_2


# estimation + benchmark
#function run_mod_once()
#    context = @dynare "HS_nk_loglin.mod"   
#    return context
#end

# Warm-up to compile
#run_mod_once();

# Benchmark full execution of the .mod file (after warm-up)
#@btime run_mod_once();