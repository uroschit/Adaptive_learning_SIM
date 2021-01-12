### Description:
	# this script implements the computation of the IRRCE;


### code
	function IRRCE_sim(params; storage_stamp::String = params.ident)
		exe_log = ["log-file for execution of 'IRRCE_sim' [master]"]
		push!(exe_log, "[$(storage_stamp)] equilibrium computation: ")
		println("[$(storage_stamp)] equilibrium computation: ")

		# Loading & setup 
			push!(exe_log, "[$(storage_stamp)] loading; ")
			println("[$(storage_stamp)] loading; ")
			g_star::Array{Float64,4} = load(string("IRRCE_solution/internal_data/[", params.ident, "] policy.jld2"))["g_star"]
			a_star					 = [(1+r)*a + params.w(r)*z - g_star[a_ind, r_ind, z_ind, d_ind]   
										for (a_ind, a) in enumerate(params.a_grid), (r_ind, r) in enumerate(params.r_grid), 
											(z_ind, z) in enumerate(params.z_grid), (d_ind, d) in enumerate(params.d_grid)]
			ā_star 					 = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid, params.d_grid), a_star, 
										extrapolation_bc = Line())

		# pre-allocations
			push!(exe_log, "[$(storage_stamp)] computing; ")
			println("[$(storage_stamp)] computing; ")
			aggr_state       = zeros(Float64, 4, params.T) 				# aggreg state is 4-dim: K_t, r_t, w_t, d_{t-1}
			dist 		     = zeros(Float64, params.a_distnodes_size * params.z_size, params.T)		# sequence of household distributions
			aggr_state[:, 1] = [params.a_init, params.r_init, params.w(params.r_init), params.d_init]
			for z_ind in 1:params.z_size		# implementing mass point at a_init in first distribution 
				dist[(z_ind-1)*params.a_distnodes_size + maximum(findall(params.a_distnodes .≤ params.a_init)), 1] = params.z_statdist[z_ind]
			end

		# updating forward
			push!(exe_log, "[$(storage_stamp)] updating forward; ")
			println("[$(storage_stamp)] updating forward; ")
			for t in 2:params.T
				prev_ā_policy(a, z)  = ā_star(a, aggr_state[2,t-1], z, aggr_state[4,t-1])
				transmat 			 = trans_mat(params.Π, params.a_distnodes, params.z_grid, prev_ā_policy)
				dist[:,t] 			 = transmat*dist[:,t-1]
				aggr_state[1,t] 	 = dot(dist[:,t], repeat(params.a_distnodes, params.z_size))
				aggr_state[2,t] 	 = params.α * params.A * aggr_state[1,t]^(params.α - 1.0) - params.δ
				aggr_state[3,t] 	 = params.w(aggr_state[2,t])
				aggr_state[4,t] 	 = aggr_state[4,t-1] - (params.σθ_sq/2) + params.g * (log(1+aggr_state[2,t]) - params.ϕ*
										log(1+aggr_state[2,t-1]) - (1-params.ϕ)*aggr_state[4,t-1] + (params.σr_sq + params.σθ_sq)/2)
				if mod(t, 50) ≈ 0 
					push!(exe_log, "[$(storage_stamp); eqm computation] done: $t of $(params.T);")
					println("[$(storage_stamp); eqm computation] done: $t of $(params.T);")
				end 
			end
		
		# saving results
			push!(exe_log, "[$(storage_stamp); eqm computation] storing; ")
			println("[$(storage_stamp); eqm computation] storing; ")
			save(string("IRRCE_solution/internal_data/[", storage_stamp, "] equilibrium.jld2"), "aggr_state", aggr_state, 
					"dist", dist, "exe_log", exe_log)

		println("[$(storage_stamp)] equilibrium done, results stored.")
	end
