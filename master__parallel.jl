### Description:
	# this file is the master - executing it calls all relevant files, so no single script
	# must be called;
	# note that the section loading packages to workers must be executed twice (first run yields an error)
	# it is normal for the "Packages and setup"-part to take a while (1-5 min); first run might take >10 min
	# while parameter setup must always be executed, Integration-nodes and the PF (all for
	# given parametrization) are stored in the resp "internal_data" folders and can just be loaded;
	# There are naturally dependencies between different sections (i.e. loading data), but as long 
	# as one executes the code sequentially, everything should be fine.
	# Contents:
		# 1) Setting up Julia environment
		# 2) Loading model parameters & elementary functions 
		# 3) Doing PFI for the IRRCE-household-problem
		# 4) Computing the IRRCE 
		# 5) Solving for the stationary Aiyagari equilibrium 
		# 6) Computing IRRCE when started at Aiyagari equilibrium 
		# 7) Computing IRFs to a one-off MIT shock to TFP 
		# 8) Generating graphs and tables for comparison of stat eqa 
		# 9) Generating graphs for IRFs 
	# Dependencies (x ← y ⇔ "x requires running y"; x ← y* ⇔ "x requires only stored results from y")
		# 2 ← 1; 3 ← 1,2; 4 ← 1,2,3*; 5 ← 1,2; 6 ← 1,2,3*,5*; 7 ← 1,2,3*,4*,5*; 8 ← 1,2,4*,5*,6*; 9 ← 1,2,4*,5*,7*;


### 1) Packages and setup
	using Pkg
	Pkg.activate("C:\\PATH")
	Pkg.instantiate()
	using QuantEcon, LinearAlgebra, Statistics, Plots, Interpolations, NLsolve, Optim, Random, BenchmarkTools, Parameters, Roots, LaTeXStrings, Distributions, JLD2, FileIO
	cd("C:/PATH")

	# Setup for parallelization
	using Distributed
	# Sys.cpu_info()		# give info about the number of virtual CPU cores that can be added
	addprocs(5; exeflags="--project")
	@everywhere begin 			# for some reason, this section needs to be executed twice to make sure all workers have the Packages
		using Pkg
		Pkg.activate("C:\\PATH")
		Pkg.instantiate()
		using QuantEcon, LinearAlgebra, Statistics, Plots, Interpolations, NLsolve, Optim, Random, BenchmarkTools, Parameters, Roots, LaTeXStrings, Distributions, JLD2, FileIO
	end


### 2) Setup: loading parameters (must always be done) and overwriting shocks and integration nodes & weights; 

	# This must always be executed
	@everywhere include("setup/Params.jl")
	include("setup/TransMat.jl")

	# out-comment to prevent existing structures from being overwritten
	include("setup/Quad_nodes_weights.jl")
	Quad_nodes_weights(calib_1)
	Quad_nodes_weights(calib_2)
	Quad_nodes_weights(calib_3)
	Quad_nodes_weights(calib_4)
	Quad_nodes_weights(calib_5)
	Quad_nodes_weights(calib_6)


### 3) Solving household problem
	# Note: Depending on size of grid, this can take several hours
	# (~ 1hr 37min per calib for: gridsize (a,r,z,d): 100×100×7×100 = 7mn on Dell XPS15 w/ 16 GB RAM, Intel i7-8750H 
	# CPU @ 2.20GHz, 2208 Mhz, 6 Core(s), 12 Logical Processors)
	
	include("IRRCE_solution/PFI__parallel.jl")
	compute_PF__parallel(calib_1, load_init = true)
	compute_PF__parallel(calib_2, load_init = true)
	compute_PF__parallel(calib_3, load_init = true)
	compute_PF__parallel(calib_4, load_init = true)
	compute_PF__parallel(calib_5, load_init = true)
	compute_PF__parallel(calib_6, load_init = true)


### 4) Finding equilibrium 
	# Purpose: for showing convergence to Aiyagari eqa from 'arbitrary' initial rate 
	# (i.e. to show learning)

	include("IRRCE_solution/equilibrium.jl")
	IRRCE_sim(calib_1)
	IRRCE_sim(calib_2)
	IRRCE_sim(calib_3)
	IRRCE_sim(calib_4)
	IRRCE_sim(calib_5)
	IRRCE_sim(calib_6)


### 5) Solving the Aiyagari model

	include("aiyagari_solution/Aiy_PFI__parallel.jl")
		# Note: Policy function on params.r_grid (for comparison) - this includes r's above λ, but PF_r 
		# are independent from one another for different r, so no problem
		# Note: this step is not necessary for solving the equilibrium since the bisection
		# calls PFI on its own -- this computation is purely for illustrative purposes, i.e. 
		# to analyze the policy-function of households in an RE-economy (for Aiyagari, with 
		# perceived fix rate r)
	compute_PF_aiy__parallel(calib_1)
	compute_PF_aiy__parallel(calib_2)
	compute_PF_aiy__parallel(calib_3)
	# compute_PF_aiy__parallel(calib_4)	# here, only parameters from IRRCE changed
	# compute_PF_aiy__parallel(calib_5)
	compute_PF_aiy__parallel(calib_6)
	#
	include("aiyagari_solution/Aiy_equilibrium.jl")
	find_eqm_aiy(calib_1)
	find_eqm_aiy(calib_2)
	find_eqm_aiy(calib_3)
	# find_eqm_aiy(calib_4)	# here, only parameters from IRRCE changed
	# find_eqm_aiy(calib_5)
	find_eqm_aiy(calib_6)


### 6) Computing equilibrium: initiated at Aiyagari SREE
	# Purpose: for welfare comparisons: what if agents switched from Aiyagari equilibrium 
	# to IRRCE?
	# Note: be careful not to change the identifier 'ident',
	# so that loadings inside IRRCE_sim execute correctly!

	# defining new calibration objects
	const calib_1_fromAiy = Setup(ident = "calib_1", r_init = load("aiyagari_solution/internal_data/[calib_1] aiy_equilibrium.jld2")["aiy_r_star"]);
	const calib_2_fromAiy = Setup(ident = "calib_2", r_init = load("aiyagari_solution/internal_data/[calib_2] aiy_equilibrium.jld2")["aiy_r_star"],
									σz_sq = (0.2)^2, ρz = 0.6);
	const calib_3_fromAiy = Setup(ident = "calib_3", r_init = load("aiyagari_solution/internal_data/[calib_3] aiy_equilibrium.jld2")["aiy_r_star"], 
									σ = 2.0);
	const calib_4_fromAiy = Setup(ident = "calib_4", r_init = load("aiyagari_solution/internal_data/[calib_3] aiy_equilibrium.jld2")["aiy_r_star"], 
									σθ_sq = (1/3) * (60.0 / 10_000)^2);
	const calib_5_fromAiy = Setup(ident = "calib_5", r_init = load("aiyagari_solution/internal_data/[calib_3] aiy_equilibrium.jld2")["aiy_r_star"], 
									ϕ = 0.983, σr_sq = (58.17 / 10_000)^2);
	const calib_6_fromAiy = Setup(ident = "calib_6", r_init = load("aiyagari_solution/internal_data/[calib_3] aiy_equilibrium.jld2")["aiy_r_star"], 
									b = 3.0);
	# computing 
	include("IRRCE_solution/equilibrium.jl")
	IRRCE_sim(calib_1_fromAiy, storage_stamp = "calib_1_fromAiy")
	IRRCE_sim(calib_2_fromAiy, storage_stamp = "calib_2_fromAiy")
	IRRCE_sim(calib_3_fromAiy, storage_stamp = "calib_3_fromAiy")
	IRRCE_sim(calib_4_fromAiy, storage_stamp = "calib_4_fromAiy")
	IRRCE_sim(calib_5_fromAiy, storage_stamp = "calib_5_fromAiy")
	IRRCE_sim(calib_6_fromAiy, storage_stamp = "calib_6_fromAiy")


### 7) Solving for reactions to one-off MIT-shock

	include("IRRCE_solution/PFIforIRF__parallel.jl")
		# Note: Depending on size of grid, this can take several hours
		# (~ 11 hrs per calib for: gridsize (a,r,z,A,d): 100×100×7×10×29 = 20mn on Dell XPS15 w/ 16 GB RAM, Intel i7-8750H 
		# CPU @ 2.20GHz, 2208 Mhz, 6 Core(s), 12 Logical Processors)
	PFIforIRF__parallel(calib_1, load_init = true)
	PFIforIRF__parallel(calib_2, load_init = true)
	PFIforIRF__parallel(calib_3, load_init = true)
	PFIforIRF__parallel(calib_4, load_init = true)
	PFIforIRF__parallel(calib_5, load_init = true)
	PFIforIRF__parallel(calib_6, load_init = true)
	include("IRRCE_solution/IRF.jl")
	IRF(calib_1)
	IRF(calib_2)
	IRF(calib_3)
	IRF(calib_4)
	IRF(calib_5)
	IRF(calib_6)
	include("aiyagari_solution/Aiy_IRF.jl")
	IRF_Aiy(calib_1)
	IRF_Aiy(calib_2)
	IRF_Aiy(calib_3)
	# IRF_Aiy(calib_4)	# here, only parameters from IRRCE changed
	# IRF_Aiy(calib_5)
	IRF_Aiy(calib_6)


### 8) Analysis: comparing stationary equilibria

	using PGFPlotsX, LaTeXTabulars, StatsBase, Formatting
	include("analysis/stat_eqa_comp.jl")
	equilibrium_graph(calib_1, path = "C:/PATH2LATEX")
	equilibrium_graph(calib_2, path = "C:/PATH2LATEX")
	equilibrium_graph(calib_3, path = "C:/PATH2LATEX")
	equilibrium_graph(calib_4, path = "C:/PATH2LATEX")
	equilibrium_graph(calib_5, path = "C:/PATH2LATEX")
	equilibrium_graph(calib_6, path = "C:/PATH2LATEX")
	# 
	key_eqm_stats(calib_1, path = "C:/PATH2LATEX/tables")
	key_eqm_stats(calib_2, path = "C:/PATH2LATEX/tables")
	key_eqm_stats(calib_3, path = "C:/PATH2LATEX/tables")
	key_eqm_stats(calib_4, path = "C:/PATH2LATEX/tables")
	key_eqm_stats(calib_5, path = "C:/PATH2LATEX/tables")
	key_eqm_stats(calib_6, path = "C:/PATH2LATEX/tables")

### 9) Analysis: comparing post-MIT adjustment paths

	using PGFPlotsX, LaTeXTabulars, StatsBase, Formatting
	include("analysis/IRF_comp.jl")
	IRF_graph(calib_1, path = "C:/PATH2LATEX")
	IRF_graph(calib_2, path = "C:/PATH2LATEX")
	IRF_graph(calib_3, path = "C:/PATH2LATEX")
	IRF_graph(calib_4, path = "C:/PATH2LATEX")
	IRF_graph(calib_5, path = "C:/PATH2LATEX")
	IRF_graph(calib_6, path = "C:/PATH2LATEX")