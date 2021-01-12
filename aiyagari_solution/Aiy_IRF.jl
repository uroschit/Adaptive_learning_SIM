### Description:
	# this file computes the transition path after an unannounced, one-off 
	# shock to TFP (parameter A), for the Aiyagari economy, and stores the 
	# results; the approach used on computing the IRF is that of adaptive 
	# dampening, inspired by the Auberbach-Kotlikoff Fixed-Point Solver with 
	# Adaptive Dampening; source: https://notes.quantecon.org/submission/5b3faf1fb9eab00015b89f9a

### code 

	function EulerBack(params, policy′, r, r′, w)
		# defining integral
			policy_interp′ = LinearInterpolation((params.a_grid, params.z_grid), policy′, 
								extrapolation_bc = Line());
			function expec_aiy(a′, z)
				z_ind = minimum(findall(params.z_grid .== z))
				expec_vals = [params.β * policy_interp′(a′, z′)^(-params.σ) * (1+r′) * 
											params.Π[i,z_ind] 		# this line contains the weights
								for (i, z′) in enumerate(params.z_grid)]
				return sum(expec_vals)
			end

		# computing 
			g = zeros(Float64, params.a_size, params.z_size)
			for (z_ind,z) in enumerate(params.z_grid)
				# Aux function 
				expec_a(a′) = expec_aiy(a′, z)
				# Find endogenous grid for slackness region
				c_endo = (expec_a.(params.a_grid)).^(- 1/(params.σ))
				a_endo = ((1+r)^(-1)).*(c_endo + params.a_grid .- w*z)
				c_interp = LinearInterpolation(a_endo, c_endo, extrapolation_bc = Line())
				# compute policy
				for (a_ind,a) in enumerate(params.a_grid)
					if a < minimum(a_endo) 		# constraint is binding
						g[a_ind, z_ind] = (1+r)*a + w*z - minimum(params.a_grid)		# note: minimum(params.a_grid) = -b
					else						# constraint is slack 
						g[a_ind, z_ind] = c_interp(a)
					end  
				end 
			end

		return g
	end;

	function update_backward(K_path_guess, params, wage, interest, A_path, c_ss)
		# generate state variable paths 
			w_path = [wage(K_path_guess[t], A_path[t])   for t in 1:params.T]
			r_path = [interest(K_path_guess[t], A_path[t])   for t in 1:params.T]

		# pre-allocate policy arrays 
			c_path = zeros(Float64, params.a_size, params.z_size, params.T)
			c_path[:,:,end] = c_ss
			a_path = zeros(Float64, params.a_size, params.z_size, params.T)
			a_path[:,:,end] = [(1+r_path[end])*a + w_path[end]*z - c_ss[a_ind, z_ind]
								for (a_ind, a) in enumerate(params.a_grid), (z_ind, z) in enumerate(params.z_grid)]

		# iterate backwards
			for t in params.T:(-1):2
				c_path[:,:,t-1] = EulerBack(params, c_path[:,:,t], r_path[t-1], r_path[t], w_path[t-1])
				a_path[:,:,t-1] = [(1+r_path[t-1])*a + w_path[t-1]*z - c_path[a_ind, z_ind ,t-1]
									for (a_ind, a) in enumerate(params.a_grid), (z_ind, z) in enumerate(params.z_grid)]
			end

		return a_path 
	end;

	function update_forward(a_path, dist_ss, params)
		# pre-allocate distribution and K arrays
			K_path_new = zeros(Float64, params.T)
			K_path_new[1] = dot(dist_ss, repeat(params.a_distnodes, params.z_size))
			dist_path_new = zeros(Float64, length(dist_ss), params.T)
			dist_path_new[:,1] = dist_ss

		# update forward 
			for t in 2:params.T 
				ā_policy = LinearInterpolation((params.a_grid, params.z_grid), a_path[:,:,t-1], extrapolation_bc = Line())
				transmat = trans_mat(params.Π, params.a_distnodes, params.z_grid, ā_policy)
				dist_path_new[:,t] = transmat*dist_path_new[:,t-1]
				K_path_new[t] = dot(dist_path_new[:,t], repeat(params.a_distnodes, params.z_size))
			end

		return K_path_new, dist_path_new
	end;

	function IRF_Aiy(params; maxiter = 60, tol = 1e-4, ω = 0.2, expand = 1.1, contract = 0.5)
		exe_log = ["log-file for execution of 'IRF_Aiy' [master]"]
		push!(exe_log, "[$(params.ident); Aiyagari] computing IRF:")
		println("[$(params.ident); Aiyagari] computing IRF:")
		
		# Loading & Setup
			push!(exe_log, "[$(params.ident); Aiyagari] loading;")
			println("[$(params.ident); Aiyagari] loading;")
			aiy_r_star::Float64 		= load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_r_star"]
			aiy_K::Float64 				= load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_K"]
			aiy_dist::Array{Float64,1}  = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_dist"]
			aiy_PF::Array{Float64,2} 	= load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_PF"]
			wage(K, A)          		= A * params.A * (1-params.α) * K^(params.α)
			interest(K, A)            	= A * params.α * params.A * K^(params.α - 1.0) - params.δ
		
		# pre-allocating output arrays 
			aiy_aggr_state_IRF = zeros(Float64, 4, params.T)	# aggreg state is 4-dim: K_t, r_t, w_t, A_shock_t
			aiy_aggr_state_IRF[4,1] = params.A_shock 
			for t in 2:params.T 
				aiy_aggr_state_IRF[4,t] = exp(params.ρA * log(aiy_aggr_state_IRF[4,t-1]))
			end
			a_policy_path = zeros(Float64, params.a_size, params.z_size, params.T)
			aiy_dist_IRF = zeros(Float64, length(aiy_dist), params.T)
		
		# computing 
			push!(exe_log, "[$(params.ident); Aiyagari] computing transition path;")
			println("[$(params.ident); Aiyagari] computing transition path;")
			dist = tol + 1.0; dist_old = tol + 1.0; iter = 0; damp = ω; K_path = repeat([aiy_K], params.T)
			for i in 1:(maxiter+1)
				iter += 1
				a_path = update_backward(K_path, params, wage, interest, aiy_aggr_state_IRF[4,:], aiy_PF)
				K_path_new = update_forward(a_path, aiy_dist, params)[1]
				dist = maximum(abs.(view(K_path, :) - view(K_path_new, :)))
				if dist > dist_old	# decrease dampening factor to avoid divergence
					damp = max(min(damp * contract, 1.0-eps()), eps())
				else 	# increase dampening factor to speed up convergence 
					damp = max(min(damp * expand, 1.0-eps()), eps())
				end 
				if dist < tol 
					K_path = K_path_new
					aiy_aggr_state_IRF[1,:] = K_path 
					aiy_aggr_state_IRF[2,:] = [interest(aiy_aggr_state_IRF[1,t], aiy_aggr_state_IRF[4,t])   for t in 1:params.T]
					aiy_aggr_state_IRF[3,:] = [wage(aiy_aggr_state_IRF[1,t], aiy_aggr_state_IRF[4,t])   for t in 1:params.T]
					a_policy_path = update_backward(K_path, params, wage, interest, aiy_aggr_state_IRF[4,:], aiy_PF)
					aiy_dist_IRF = update_forward(a_policy_path, aiy_dist, params)[2] 
					push!(exe_log, "[$(params.ident); Aiyagari] K-path done: dist = $(dist), iter = $(iter);")
					println("[$(params.ident); Aiyagari] K-path done: dist = $(dist), iter = $(iter);")
					break
				elseif i == (maxiter+1)
					K_path = K_path_new
					aiy_aggr_state_IRF[1,:] = K_path 
					aiy_aggr_state_IRF[2,:] = [interest(aiy_aggr_state_IRF[1,t], aiy_aggr_state_IRF[4,t])   for t in 1:params.T]
					aiy_aggr_state_IRF[3,:] = [wage(aiy_aggr_state_IRF[1,t], aiy_aggr_state_IRF[4,t])   for t in 1:params.T]
					a_policy_path = update_backward(K_path, params, wage, interest, aiy_aggr_state_IRF[4,:], aiy_PF)
					aiy_dist_IRF = update_forward(a_policy_path, aiy_dist, params)[2] 
					push!(exe_log, "[$(params.ident); Aiyagari] K-path no convergence at maxiter = $(maxiter): dist = $(dist);")
					println("[$(params.ident); Aiyagari] K-path no convergence at maxiter = $(maxiter): dist = $(dist);")
				else 
					K_path = damp*K_path_new + (1.0-damp)*K_path
					dist_old = dist 
					push!(exe_log, "[$(params.ident); Aiyagari] K-path iter $(iter): dist = $(dist), damp = $(damp);")
					println("[$(params.ident); Aiyagari] K-path iter $(iter): dist = $(dist), damp = $(damp);")
				end 
			end

		# Storing 
			save(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium_IRF.jld2"),
				"aiy_aggr_state_IRF", aiy_aggr_state_IRF, "aiy_dist_IRF", aiy_dist_IRF, "a_policy_path", a_policy_path, "exe_log", exe_log)

		println("[$(params.ident); Aiyagari] IRF done.")
	end;