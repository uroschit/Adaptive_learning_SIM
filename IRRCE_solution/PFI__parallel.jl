### Description:
	# this file executes the time iteration 
	# as the suffix "__parallel" indicates, this code is parallelized on multiple cores of 
	# the CPU; it is essential to define all relevant functions with the "@everywhere"-macro 


### code 
	function compute_PF__parallel(params; load_init::Bool = false)
		exe_log = vcat([["log-file for execution of 'compute_PF__parallel' [master]"]], 
					[["log-file for execution of 'compute_PF__parallel' [d_ind = $(d_ind)]"] for d_ind in 1:params.r_size])
		push!(exe_log[1], "[$(params.ident)] computing policy function:")
		println("[$(params.ident)] computing policy function:")

		# Loading 
			push!(exe_log[1], "[$(params.ident)] loading;")
			println("[$(params.ident)] loading; ")
			n_simps::Array{Float64,1}	= load(string("IRRCE_solution/internal_data/[", params.ident, "] Quad_nodes.jld2"))["n_simps"]
			w_simps::Array{Float64,1}	= load(string("IRRCE_solution/internal_data/[", params.ident, "] Quad_nodes.jld2"))["w_simps"]
		
		# Setting up
			push!(exe_log[1], "[$(params.ident)] preallocating policy array and interm container; ")
			println("[$(params.ident)] preallocating policy array and interm container; ")
			g_star = zeros(Float64, params.a_size, params.r_size, params.z_size, params.r_size)
			container = [zeros(Float64, params.a_size, params.r_size, params.z_size) for d in params.d_grid]
			if load_init	# if available, load output from past computations to speed up convergence 
				g_init::Array{Float64,4} = load(string("IRRCE_solution/internal_data/[", params.ident, "] policy.jld2"))["g_star"]
			else 
				g_init = [max((1+r)*a + params.w(r)*z, 1e-4)
							for a in params.a_grid, r in params.r_grid, z in params.z_grid, d in params.d_grid]
			end 

		# Defining Coleman operator & associated functions
			push!(exe_log[1], "[$(params.ident)] defining Coleman operator; ")
			println("[$(params.ident)] defining Coleman operator; ")
			# aux functions 
				@inline function r′(R̃′, r, d)
					return exp((2*(params.ξ*(1- params.ϕ)^2 + params.σr_sq + params.σθ_sq))^(1/2) * R̃′ + params.ϕ*log(1+r) + 
											(1- params.ϕ)*d - (params.σr_sq + params.σθ_sq)/2) - 1.0
				end 
			# Define integral
				function expec(a′, r, z, d, policy_d)
					z_ind = minimum(findall(params.z_grid .== z))
		            expec_vals = zeros(Float64, params.n_nodes_simps, params.z_size)	# pre-allocating 
		            # Defining lower bound of censored integration domain in terms of transformed variable
		            critval = (log(1-params.δ) - params.ϕ*log(1+r) - (1- params.ϕ)*d + (params.σr_sq + params.σθ_sq)/2)/
		                        ((2*(params.ξ*(1- params.ϕ)^2 + params.σr_sq + params.σθ_sq))^(1/2))
		            # Finding lower bound on grid of Simpson integration nodes 
		            n_simps_crit = minimum(findall(n_simps.>critval)) == 1 ? 1 : minimum(findall(n_simps.>critval)) - 1
		            # Computing values at nodes 
		            for (i, R̃′) in enumerate(n_simps[n_simps_crit:end])
		                for (j, z′) in enumerate(params.z_grid)
		                    expec_vals[i,j] = params.β * policy_d(a′, r′(R̃′, r, d), z′)^(-params.σ) *
		                    					(1.0+r′(R̃′, r, d)) * exp(- (R̃′)^2) * π^(-1/2) *
		                                        w_simps[i] * params.Π[j,z_ind] 		# this line contains the weights
		                end
		            end 
				    return sum(expec_vals)
		        end
			# Defining Coleman operator (d fixed); using EGM
				function Coleman_fix_d(g_d, d) 
					Kg_d = similar(g_d) 	# preallocating
					# interpolate input policy
					ḡ_d = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid), g_d, 
							extrapolation_bc = Line()); 
					for ((r_ind,r), (z_ind,z)) in Base.Iterators.product(enumerate(params.r_grid), enumerate(params.z_grid))
						# Aux function 
						expec_a(a′) = expec(a′, r, z, d, ḡ_d)
						# Find endogenous grid for slackness region
						c_endo = (expec_a.(params.a_grid)).^(- 1.0/(params.σ))
						a_endo = ((1.0+r)^(-1.0)).*(c_endo + params.a_grid .- params.w(r)*z)
						c_interp = LinearInterpolation(a_endo, c_endo, extrapolation_bc = Line())
						# compute policy
						for (a_ind,a) in enumerate(params.a_grid)
							if a < minimum(a_endo) 		# constraint is binding
								Kg_d[a_ind, r_ind, z_ind] = (1.0+r)*a + params.w(r)*z - minimum(params.a_grid)		# note: minimum(params.a_grid) = -b
							else						# constraint is slack 
								Kg_d[a_ind, r_ind, z_ind] = c_interp(a)
							end  
						end 
					end
					return Kg_d
				end

		# Defining time iteration (outer loop)
			function time_iteration_fix_d(d, g_init)
				d_ind = minimum(findall(x -> x == d, params.d_grid))	
					# note: minimum() needed so that output is Int64 (needed for creation of g_d)
				dist = params.tol + 1.0
				g_d = copy(g_init[:, :, :, d_ind])
				for i in 1:params.maxiter
					if dist < params.tol 
						push!(exe_log[d_ind+1], "[$(params.ident); d_ind = $d_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")
						println("[$(params.ident); d_ind = $d_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")
						break 
					end 
					g_d_new = Coleman_fix_d(g_d, d)
					dist = maximum(abs.(view(g_d_new, :, :, :) - view(g_d, :, :, :)))
					g_d = g_d_new
					if i == params.maxiter
						push!(exe_log[d_ind+1], "[$(params.ident); d_ind = $d_ind] no convergence at i = $i; dist = $dist; ")
						println("[$(params.ident); d_ind = $d_ind] no convergence at i = $i; dist = $dist; ")
					elseif mod(i, 50) ≈ 0 || i < 4	# update about status only every X iterations to declutter feed
						push!(exe_log[d_ind+1], "[$(params.ident); d_ind = $d_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")
						println("[$(params.ident); d_ind = $d_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")
					end 
				end
				return g_d
			end
			# fixing initial g 
			time_iteration_fix_d_given_g(d) = time_iteration_fix_d(d, g_init)

		# Computing
			push!(exe_log[1], "[$(params.ident)] starting time iteration; ")
			println("[$(params.ident)] starting time iteration; ")
			@time container = pmap(time_iteration_fix_d_given_g, params.d_grid)
			for d_ind in 1:params.r_size
				g_star[ :, :, :, d_ind] = container[d_ind]
			end	

		# Storing
			push!(exe_log[1], "[$(params.ident)] saving results; ")
			println("[$(params.ident)] saving results; ")
			save(string("IRRCE_solution/internal_data/[", params.ident, "] policy.jld2"), "g_star", g_star, "g_init", g_init, "exe_log", exe_log)
		
		println("[$(params.ident)] results saved; operation completed.")
	end