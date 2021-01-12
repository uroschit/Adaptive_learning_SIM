### Description:
	# This script executes the bisection method to find the stationary 
	# equilibrium of the Aiyagari model 


### code 
	function find_eqm_aiy(params)	
		exe_log = ["log-file for execution of 'find_eqm_aiy' [master]"]
		push!(exe_log, "[Aiyagari; $(params.ident)] equilibrium computation:")
		println("[Aiyagari; $(params.ident)] equilibrium computation:")	

		# Defining Coleman operator & associated functions
			push!(exe_log, "[Aiyagari; $(params.ident)] defining Coleman operator; ")
			println("[Aiyagari; $(params.ident)] defining Coleman operator; ")
			# Defining integral
				function expec_aiy(a′, r, z, policy_r)
					z_ind = minimum(findall(params.z_grid .== z))
					expec_vals = [params.β * policy_r(a′, z′)^(-params.σ) * (1+r) * 
												params.Π[i,z_ind] 		# this line contains the weights
									for (i, z′) in enumerate(params.z_grid)]
					return sum(expec_vals)
				end
	        # Defining Coleman operator (r fixed); using EGM
	        	function Coleman_fix_r(g_r, r) # input: policy-function, interest rate 
					Kg_r = similar(g_r) 	# preallocating
					# interpolate input policy
					ḡ_r = LinearInterpolation((params.a_grid, params.z_grid), g_r, 
							extrapolation_bc = Line()); 
					for (z_ind,z) in enumerate(params.z_grid)
						# Aux function 
						expec_a(a′) = expec_aiy(a′, r, z, ḡ_r)
						# Find endogenous grid for slackness region
						c_endo = (expec_a.(params.a_grid)).^(- 1/(params.σ))
						a_endo = ((1+r)^(-1)).*(c_endo + params.a_grid .- params.w(r)*z)
						c_interp = LinearInterpolation(a_endo, c_endo, extrapolation_bc = Line())
						# compute policy
						for (a_ind,a) in enumerate(params.a_grid)
							if a < minimum(a_endo) 		# constraint is binding
								Kg_r[a_ind, z_ind] = (1+r)*a + params.w(r)*z - minimum(params.a_grid)		# note: minimum(params.a_grid) = -b
							else						# constraint is slack 
								Kg_r[a_ind, z_ind] = c_interp(a)
							end  
						end 
					end
					return Kg_r
				end

		# Defining time time iteration	
			function time_iteration(r)
				g_r = [(1+r)*a + params.w(r)*z 
						for a in params.a_grid, z in params.z_grid]
				dist = params.tol + 1.0
				for i in 1:params.maxiter
					if dist < params.tol 
						break 
					end 
					g_r_new = Coleman_fix_r(g_r, r)
					dist = maximum(abs.(view(g_r_new, :, :) - view(g_r, :, :)))
					g_r = g_r_new
				end
				return g_r
			end

		# Defining long run asset supply 
			function longrun_assets(r; maxiter_transmat = 1_000, tol_transmat = 1e-10)
				# find policy 
				c_policy = time_iteration(r)
				a_policy = [(1+r)*a + params.w(r)*z - c_policy[a_ind, z_ind]   
							for (a_ind, a) in enumerate(params.a_grid), (z_ind, z) in enumerate(params.z_grid)];
				ā_policy = LinearInterpolation((params.a_grid, params.z_grid), a_policy, extrapolation_bc = Line());
				# find transition matrix 
				Γ = trans_mat(params.Π, params.a_distnodes, params.z_grid, ā_policy)
				# arbitrary initial distribution 
				λ = ones(Float64, params.a_distnodes_size * params.z_size)/(params.a_distnodes_size * params.z_size)
				# iterating on transition matrix until convergence
				iter = 0; Δ = tol_transmat + 1.0
				for i in 1:maxiter_transmat
					iter += 1
					if Δ < tol_transmat
						break
					end 
					λ_new = Γ*λ
					Δ = maximum(abs.(λ - λ_new))
						# Note: be careful with selecting tol_transmat! Setting it to 1e-5 is already too generous: it 
						# makes the above line converge too early so that LR assets are beyond O(10^-1) of 
						# the mean-K that results from iterating over Γλ 1_000 times;
					λ = λ_new
				end
				# return result				
				return dot(repeat(params.a_distnodes, params.z_size), λ), λ, Γ
			end

		# Defining bisection
			function bisection(supply, demand; tol = 1e-4, maxiter = 25)
				low = minimum(params.r_grid); high = params.λ - 1e-4; guess_val = 10.0; iter = 0; r = 0.0
				diff(r) = supply(r) - demand(r)
				for i in 1:maxiter
					iter += 1
					r = 0.5*(high + low)
					guess_val = diff(r)
					if abs(guess_val) < tol 
						push!(exe_log, "[Aiyagari; $(params.ident)] eqm done: r_star = $r, K-dist = $(abs(guess_val)), iter = $iter. ")
						println("[$(params.ident); Aiyagari, maxit = $maxiter] eqm done: r_star = $r, K-dist = $(abs(guess_val)), iter = $iter.")
						break 
					elseif signbit(guess_val) == true 	# signbit(-1) = true
						low = r 
						push!(exe_log, "[Aiyagari; $(params.ident)] iter $iter: r = $r, K-dist = $(abs(guess_val)), corr. bound = lower; ")
						println("$(params.ident); iter $iter: r = $r, K-dist = $(abs(guess_val)), corr. bound = lower;")
					else 
						high = r 
						push!(exe_log, "[Aiyagari; $(params.ident)] iter $iter: r = $r, K-dist = $(abs(guess_val)), corr. bound = upper; ")
						println("$(params.ident); iter $iter: r = $r, K-dist = $(abs(guess_val)), corr. bound = upper;")
					end
				end 
				return r
			end

		# Computing
			longrun_assets_onlyK(r) = longrun_assets(r)[1]
			push!(exe_log, "[Aiyagari; $(params.ident)] computing; ")
			println("[Aiyagari; $(params.ident)] computing; ")	
			@time aiy_r_star = bisection(longrun_assets_onlyK, 
					(r -> ((r+params.δ)/(params.α*params.A))^(-1.0/(1.0-params.α))))
			aiy_K = ((aiy_r_star+params.δ)/(params.α*params.A))^(-1.0/(1.0-params.α))
			aiy_dist = longrun_assets(aiy_r_star)[2]
			aiy_transmat = longrun_assets(aiy_r_star)[3]
			aiy_PF = time_iteration(aiy_r_star)
			push!(exe_log, "[Aiyagari; $(params.ident)] storing; ")
			println("[Aiyagari; $(params.ident)] storing; ")	
			save(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"), 
				"aiy_r_star", aiy_r_star, "aiy_K", aiy_K, "aiy_dist", aiy_dist, "aiy_transmat", 
				aiy_transmat, "aiy_PF", aiy_PF, "exe_log", exe_log)

		println("[Aiyagari; $(params.ident)] equilibrium computed and stored.")	
	end;