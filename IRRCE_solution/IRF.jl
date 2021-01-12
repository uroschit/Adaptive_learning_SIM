### Description:
	# this file computes the transition path after an unannounced, one-off 
	# shock to TFP (parameter A), for the IRRCE, and stores the results


### code 
	
	# Transition path IRRCE
		function IRF(params)
			exe_log = ["log-file for execution of 'IRF' [master]"]
			push!(exe_log, "[$(params.ident)] computing IRF:")
			println("[$(params.ident)] computing IRF:")

			# Loading & Setup
				push!(exe_log, "[$(params.ident)] loading data;")
				println("[$(params.ident)] loading data;")
				d_grid_small = params.d_grid[minimum(findall(params.d_grid .> 0.0)):maximum(findall(params.d_grid .< log(1.05)))]
					# this grid ↑ should be the same as the one in the PFI computation with A_shock 
				g_star_IRF::Array{Float64,5} = load(string("IRRCE_solution/internal_data/[", params.ident, "] policy_IRF.jld2"))["g_star_IRF"]
				aggr_state::Array{Float64,2} = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["aggr_state"]
				dist::Array{Float64,2}       = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["dist"]
				wage(r, A)          		 = A^(1/(1-params.α))*params.w(r)
				a_star_IRF					 = [(1+r)*a + wage(r, A)*z - g_star_IRF[a_ind, r_ind, z_ind, A_ind, d_ind]   
												for (a_ind, a) in enumerate(params.a_grid), (r_ind, r) in enumerate(params.r_grid), 
												(z_ind, z) in enumerate(params.z_grid), (A_ind, A) in enumerate(params.A_shock_grid), 
												(d_ind, d) in enumerate(d_grid_small)]
				ā_star_IRF					 = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid, params.A_shock_grid, 
												d_grid_small), a_star_IRF, extrapolation_bc = Line())

			# pre-allocating output arrays 
				push!(exe_log, "[$(params.ident)] computing; ")
				println("[$(params.ident)] computing; ")
				aggr_state_IRF       = zeros(Float64, 5, params.T) 	# aggreg state is 5-dim: K_t, r_t, w_t, d_{t-1}, A_shock_t
				dist_IRF 		     = zeros(Float64, params.a_distnodes_size * params.z_size, params.T)		# sequence of household distributions
				r_shock				 = params.A_shock * params.α * params.A * aggr_state[1,end]^(params.α - 1.0) - params.δ
				aggr_state_IRF[:, 1] = [aggr_state[1,end], r_shock,
										params.A_shock * params.A * (1-params.α) * aggr_state[1,end]^(params.α), 
										aggr_state[4,end-1] - (params.σθ_sq/2) + params.g * (log(1+r_shock) - params.ϕ*
										log(1+aggr_state[2,end-1]) - (1-params.ϕ)*aggr_state[4,end-1] + (params.σr_sq + params.σθ_sq)/2), 
										params.A_shock]
				dist_IRF[:,1]        = dist[:,end]

			# computing 
				for t in 2:params.T
					prev_ā_policy(a, z) = ā_star_IRF(a, aggr_state_IRF[2,t-1], z, aggr_state_IRF[5,t-1], aggr_state_IRF[4,t-1])
					transmat 			= trans_mat(params.Π, params.a_distnodes, params.z_grid, prev_ā_policy)
					dist_IRF[:,t] 	   	= transmat*dist_IRF[:,t-1]
					aggr_state_IRF[5,t] = exp(params.ρA * log(aggr_state_IRF[5,t-1]))
					aggr_state_IRF[1,t] = dot(dist_IRF[:,t], repeat(params.a_distnodes, params.z_size))
					aggr_state_IRF[2,t] = aggr_state_IRF[5,t] * params.α * params.A * aggr_state_IRF[1,t]^(params.α - 1.0) - params.δ
					aggr_state_IRF[3,t] = wage(aggr_state_IRF[2,t], aggr_state_IRF[5,t])
					aggr_state_IRF[4,t] = aggr_state_IRF[4,t-1] - (params.σθ_sq/2) + params.g * (log(1+aggr_state_IRF[2,t]) - params.ϕ*
											log(1+aggr_state_IRF[2,t-1]) - (1-params.ϕ)*aggr_state_IRF[4,t-1] + (params.σr_sq + params.σθ_sq)/2)
					if mod(t, 50) ≈ 0 
						push!(exe_log, "[$(params.ident); IRF computation] done: $t of $(params.T)")
						println("[$(params.ident); IRF computation] done: $t of $(params.T)")
					end 
				end

			# Storing
				push!(exe_log, "[$(params.ident); IRF computation] storing; ")
				println("[$(params.ident); IRF computation] storing; ")
				save(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium_IRF.jld2"), "aggr_state_IRF", 
						aggr_state_IRF, "dist_IRF", dist_IRF, "exe_log", exe_log)
			
			println("[$(params.ident)] IRF done.")			
		end