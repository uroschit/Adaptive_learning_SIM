### Description:
	# in this script, you can change all relevant parameters
	# be sure to keep all parameters as Float64 or Int64 
	# (depending on special case), since storing as const 
	# prohibits changing the type later 


### code
	Setup = @with_kw (  # identifier
						ident = " ",
						# Free model parameters
						r_init = 0.015,
						β = 0.96,
						σ = 3.0,
						b_adhoc = 0.0,
						α = 0.36,
						A = 1.0,
						δ = 0.08,
						ϕ = 0.90,
						σr_sq = (60.0 / 10_000)^2,	# specified as std dev in basis points;
													# loosely based on Backus et al (1998), 
													# table 1 (p. 32), 12 mo. maturity:
													# (1- .983^2)^(1/2)*316.8 ≈ 58.17 bpt
						σz_sq = (0.4)^2,
						ρz = 0.9,
						A_shock = 1.05, 	# 5% shock on whatever A is 
						ρA = (0.979)^4,
						# Free numerical parameters
						r_low = - δ + 1e-2,
						r_high = 0.1,
						r_size = 100,
						z_size = 7,
						A_size = 10,
						s_low = -2.5,	# 2.5 std dev.s cover 98.8% of p-mass 
						s_high = 2.5,
						s_h_approx = 0.25, # approx interval size for simpson integration
						tol = 1e-4,		# tol ≤ 1e-3 bcs o/w Aiyagari r* is up to 5 bpt away from actual level! 
										# values for tol and maxiter are inspired by QuantEcon:
										# http://quantecon.github.io/QuantEcon.jl/stable/api/QuantEcon.html#QuantEcon.solve-Union{Tuple{Algo},%20Tuple{QuantEcon.DiscreteDP{T,NQ,NR,Tbeta,Tind,TQ}%20where%20TQ%3C:AbstractArray{T,NQ}%20where%20Tind%20where%20Tbeta%3C:Real%20where%20NR%20where%20NQ,Type{Algo}},%20Tuple{QuantEcon.DiscreteDP{T,NQ,NR,Tbeta,Tind,TQ}%20where%20TQ%3C:AbstractArray{T,NQ}%20where%20Tind%20where%20Tbeta%3C:Real%20where%20NR%20where%20NQ},%20Tuple{T}}%20where%20T%20where%20Algo%3C:QuantEcon.DDPAlgorithm
						maxiter = 300,
						N_hh = 20_000, 
						T = 1_000, # should be at least 501, better ≥ 1_000
						# dependent model parameters:
						λ = β^(-1) - 1,
						a_init = ((r_init + δ)/(α*A))^(1/(α - 1)),
						d_init = log(1 + r_init),
						σθ_sq = (1/10) * σr_sq,
						ξ = (4*(1+ϕ^2))^(-1) * (-(1+ϕ)^2 *σθ_sq  +  ((1+ϕ)^4 * σθ_sq^2  - 8*(1+ϕ^2)*(σθ_sq^2 * ϕ^2 - σθ_sq*σr_sq))^(1/2)),
						g = ((1-ϕ)*ξ + σθ_sq) / ((1-ϕ)^2 * ξ + σθ_sq + σr_sq),
						w = r -> (1 - α) * A * ((α*A)/(r + δ))^(α / (1-α)),	
						# dependent numerical parameters:
						d_low = log(1 + r_low),
						d_high = log(1 + r_high),
						r_grid = range(r_low, r_high; length = r_size),
						d_grid = range(d_low, d_high; length = r_size),
						z_chain = tauchen(z_size, ρz, sqrt(σz_sq * (1-(ρz)^2)), 0.0),
						Π = permutedims(z_chain.p, [2, 1]),	# → watch out: in params.z_chain.p rows sum to one!
						z_statdist = stationary_distributions(z_chain)[1],
		  				z_grid = exp.(z_chain.state_values)./dot(exp.(z_chain.state_values),z_statdist), # renormalization so that Expectation is 1.0; see Aiyagari WP, footnote 33;
		  				A_shock_grid = range(1.0, A_shock, length = A_size),
		  				n_nodes_simps = ceil(Int64,(s_high - s_low)/s_h_approx) + 
										iseven(ceil(Int64,(s_high - s_low)/s_h_approx))*1,
						# asset parameters
						b = min(b_adhoc, minimum(
								[r > 0.0 ? minimum(z_grid)*w(r)/r : +Inf   for r in r_grid])),
						a_low = -b + 1e-5,
						a_high = 70.0,
						a_size = 100,		# Note: a_lownodes nodes are between 0 and 10
						a_distnodes_size = 4*a_size,
						a_grid = a_low .+ (a_high-a_low)*(range(0.0,step=0.5,length=a_size).^4.0/maximum(range(0.0,step=0.5,length=a_size).^4.0)),
						a_distnodes = a_low .+ (a_high-a_low)*(range(0.0,step=0.5,length=a_distnodes_size).^4.0/maximum(range(0.0,step=0.5,length=a_distnodes_size).^4.0)),
							# grid spacing taken from A. McKay: https://github.com/amckay/McKayLectureCodes/blob/master/Julia/Aiyagari/Aiyagari.jl - search for 'grid_fun'
		  				)

	const calib_1 = Setup(ident = "calib_1");
	const calib_2 = Setup(ident = "calib_2", σz_sq = (0.2)^2, ρz = 0.6);
	const calib_3 = Setup(ident = "calib_3", σ = 2.0);
	const calib_4 = Setup(ident = "calib_4", σθ_sq = (1/3) * (60.0 / 10_000)^2);
	const calib_5 = Setup(ident = "calib_5", ϕ = 0.983, σr_sq = (58.17 / 10_000)^2);
	const calib_6 = Setup(ident = "calib_6", ρA = 0.0);
		# declaring constant helps type inference in functions of these globals, tremendously speeding up computations



	println("parameters loaded. ")
