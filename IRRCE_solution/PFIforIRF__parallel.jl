### Description:
	# this file computes the policy function depending also 
	# on TFP (parameter A);


### code 
	
	# PFI for IRRCE with A_shock 
		function PFIforIRF__parallel(params; load_init::Bool = false)
			exe_log = vcat([["log-file for execution of 'PFIforIRF__parallel' [master]"]], 
						[["log-file for execution of 'PFIforIRF__parallel' [d_ind = $(d_ind)]"] for d_ind in 1:params.r_size])
			push!(exe_log[1], "[$(params.ident)] computing policy function with A_shock:")
			println("[$(params.ident)] computing policy function with A_shock:")

			# Loading & Setup
				push!(exe_log[1], "[$(params.ident); PFI4IRF] loading data;")
				println("[$(params.ident); PFI4IRF] loading data;")
				n_simps::Array{Float64,1}	= load(string("IRRCE_solution/internal_data/[", params.ident, "] Quad_nodes.jld2"))["n_simps"]
				w_simps::Array{Float64,1}	= load(string("IRRCE_solution/internal_data/[", params.ident, "] Quad_nodes.jld2"))["w_simps"]
				g_star::Array{Float64,4} 	= load(string("IRRCE_solution/internal_data/[", params.ident, "] policy.jld2"))["g_star"]
				# 
				wage(r, A)		            = A^(1/(1-params.α))*params.w(r)
				d_grid_small				= params.d_grid[minimum(findall(params.d_grid .> 0.0)):maximum(findall(params.d_grid .< log(1.05)))]
					# so that we don't unneccessarily compute PF for d that aren't needed for IRF 
					# if you change d_grid_small, also change it in IRF computation!

			# pre-allocating
				g_star_IRF = zeros(Float64, params.a_size, params.r_size, params.z_size, params.A_size, length(d_grid_small))
				container  = [zeros(Float64, params.a_size, params.r_size, params.z_size, params.A_size) for d in d_grid_small]
				if load_init
					g_init_IRF::Array{Float64,5} = load(string("IRRCE_solution/internal_data/[", params.ident, "] policy_IRF.jld2"))["g_star_IRF"]
				else 
					g_init_IRF = [g_star[a_ind, r_ind, z_ind, d_ind]
								for (a_ind, a) in enumerate(params.a_grid), (r_ind, r) in enumerate(params.r_grid), 
								(z_ind, z) in enumerate(params.z_grid), (A_ind, A) in enumerate(params.A_shock_grid), 
								(d_ind, d) in enumerate(d_grid_small)]
				end 
			
			# aux functions 
				@inline function r′(R̃′, r, d)
					return exp((2*(params.ξ*(1- params.ϕ)^2 + params.σr_sq + params.σθ_sq))^(1/2) * R̃′ + params.ϕ*log(1+r) + 
											(1- params.ϕ)*d - (params.σr_sq + params.σθ_sq)/2) - 1.0
				end 
		        function expec(a′, r, z, A, d, policy_d)
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
		                    expec_vals[i,j] = params.β * policy_d(a′, r′(R̃′, r, d), z′, 
		                    					exp(params.ρA * log(A)) )^(-params.σ) *
		                    					(1+r′(R̃′, r, d)) * exp(- (R̃′)^2) * π^(-1/2) *
		                                        w_simps[i] * params.Π[j,z_ind]  		# this line contains the weights
		                end
		            end 
				    return sum(expec_vals)
		        end
			
			# Coleman-Operator 
				function Coleman_fix_d(g_d, d) # inputs: policy-function, current belief parameter
					Kg_d = similar(g_d) 	# preallocating
					# interpolate input policy
					ḡ_d = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid, params.A_shock_grid), g_d, 
							extrapolation_bc = Line()); 
					for ((r_ind,r), (z_ind,z), (A_ind,A)) in Base.Iterators.product(enumerate(params.r_grid),
						enumerate(params.z_grid), enumerate(params.A_shock_grid))
						# Aux function 
						expec_a(a′) = expec(a′, r, z, A, d, ḡ_d)
						# Find endogenous grid for slackness region
						c_endo = (expec_a.(params.a_grid)).^(- 1/(params.σ))
						a_endo = ((1+r)^(-1)).*(c_endo + params.a_grid .- wage(r, A)*z)
						c_interp = LinearInterpolation(a_endo, c_endo, extrapolation_bc = Line())
						# compute policy
						for (a_ind,a) in enumerate(params.a_grid)
							if a < minimum(a_endo) 		# constraint is binding
								Kg_d[a_ind, r_ind, z_ind, A_ind] = (1+r)*a + wage(r, A)*z - 
																	minimum(params.a_grid)	# note: minimum(params.a_grid) = -b
							else						# constraint is slack 
								Kg_d[a_ind, r_ind, z_ind, A_ind] = c_interp(a)
							end  
						end 
					end
					return Kg_d
				end
			
			# time iteration
				function time_iteration_fix_d(d, g_init)
					d_ind = minimum(findall(x -> x == d, d_grid_small))	
						# note: minimum() needed so that output is Int64 (needed for creation of g_d)
					dist = params.tol + 1.0
					g_d = copy(g_init[:, :, :, :, d_ind])
					for i in 1:params.maxiter
						if dist < params.tol 
							push!(exe_log[d_ind+1], "[$(params.ident); PFI4IRF; d_ind = $d_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")
							println("[$(params.ident); PFI4IRF; d_ind = $d_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")
							break 
						end 
						g_d_new = Coleman_fix_d(g_d, d)
						dist = maximum(abs.(view(g_d_new, :, :, :, :) - view(g_d, :, :, :, :)))
						g_d = g_d_new
						if i == params.maxiter
							push!(exe_log[d_ind+1], "[$(params.ident); PFI4IRF; d_ind = $d_ind] no convergence at i = $i; dist = $dist; ")
							println("[$(params.ident); PFI4IRF; d_ind = $d_ind] no convergence at i = $i; dist = $dist; ")
						elseif mod(i, 50) ≈ 0 || i < 4	# update about status only every X iterations to declutter feed
							push!(exe_log[d_ind+1], "[$(params.ident); PFI4IRF; d_ind = $d_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")
							println("[$(params.ident); PFI4IRF; d_ind = $d_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")
						end 
					end
					return g_d
				end
				# fixing initial g 
				time_iteration_fix_d_given_g(d) = time_iteration_fix_d(d, g_init_IRF)
			
			# Computing 
				push!(exe_log[1], "[$(params.ident); PFI4IRF] starting time iteration;")
				println("[$(params.ident); PFI4IRF] starting time iteration;")
				@time container = pmap(time_iteration_fix_d_given_g, d_grid_small)
				for d_ind in 1:length(d_grid_small)
					g_star_IRF[ :, :, :, :, d_ind] = container[d_ind]
				end	
			
			# Storing
				push!(exe_log[1], "[$(params.ident); PFI4IRF] saving PFI results; ")
				println("[$(params.ident); PFI4IRF] saving PFI results; ")
				save(string("IRRCE_solution/internal_data/[", params.ident, "] policy_IRF.jld2"), "g_star_IRF", 
					g_star_IRF, "g_init_IRF", g_init_IRF, "exe_log", exe_log)

			println("[$(params.ident); PFI4IRF] PFI for IRF done.")
		end;