### Description:
	# this file executes the time iteration for Aiyagari, for
	# different values of r 


### code 
	function compute_PF_aiy__parallel(params)
		exe_log = vcat([["log-file for execution of 'compute_PF_aiy__parallel' [master]"]], 
						[["log-file for execution of 'compute_PF_aiy__parallel' [r_ind = $(r_ind)]"] for r_ind in 1:params.r_size])
		push!(exe_log[1], "[Aiyagari; $(params.ident)] computing policy function:")
		println("[Aiyagari; $(params.ident)] computing policy function:")

		# Setting up
			push!(exe_log[1], "[Aiyagari; $(params.ident)] preallocating policy array and interm container; ")
			println("[Aiyagari; $(params.ident)] preallocating policy array and interm container; ")
			aiy_g_star	= zeros(Float64, params.a_size, params.r_size, params.z_size)		# Pre-allocating 
			aiy_container = [zeros(Float64, params.a_size, params.z_size) for r in params.r_grid]
			push!(exe_log[1], "[Aiyagari; $(params.ident)] setting initial guess for g; ")			
			println("[Aiyagari; $(params.ident)] setting initial guess for g; ")
			aiy_g_init = [(1+r)*a + params.w(r)*z 
							for a in params.a_grid, r in params.r_grid, z in params.z_grid]
			
		# Defining Coleman operator & associated functions
			push!(exe_log[1], "[Aiyagari; $(params.ident)] defining Coleman operator; ")			
			println("[Aiyagari; $(params.ident)] defining Coleman operator; ")
			# Defining integral
				function expec_aiy(a′, r, z, policy_r)
					z_ind = minimum(findall(params.z_grid .== z))
					expec_vals = [params.β * policy_r(a′, z′)^(-params.σ) * (1.0+r) * 
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
						c_endo = (expec_a.(params.a_grid)).^(- 1.0/(params.σ))
						a_endo = ((1.0+r)^(-1.0)).*(c_endo + params.a_grid .- params.w(r)*z)
						c_interp = LinearInterpolation(a_endo, c_endo, extrapolation_bc = Line())
						# compute policy
						for (a_ind,a) in enumerate(params.a_grid)
							if a < minimum(a_endo) 		# constraint is binding
								Kg_r[a_ind, z_ind] = (1.0+r)*a + params.w(r)*z - minimum(params.a_grid)		# note: minimum(params.a_grid) = -b
							else						# constraint is slack 
								Kg_r[a_ind, z_ind] = c_interp(a)
							end  
						end 
					end
					return Kg_r
				end

		# Defining time iteration (outer loop)
			function time_iteration_fix_r(r, g_init)
				r_ind = minimum(findall(x -> x == r, params.r_grid))	# note: minimum() needed so that output is Int64 (needed for creation of g_d)
				dist = params.tol + 1.0
				g_r = copy(g_init[:, r_ind, :])
				for i in 1:params.maxiter
					if dist < params.tol 
						push!(exe_log[r_ind+1], "[$(params.ident); r_ind = $r_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")			
						println("[$(params.ident); r_ind = $r_ind] done with $(i-1) of $(params.maxiter), dist = $dist; ")
						break 
					end 
					g_r_new = Coleman_fix_r(g_r, r)
					dist = maximum(abs.(view(g_r_new, :, :) - view(g_r, :, :)))
					g_r = g_r_new
					if i == params.maxiter
						push!(exe_log[r_ind+1], "[$(params.ident); r_ind = $r_ind] no convergence at i = $i; dist = $dist; ")			
						println("[$(params.ident); r_ind = $r_ind] no convergence at i = $i; dist = $dist; ")
					elseif mod(i, 50) ≈ 0 || i < 4	# update about status only every X iterations to declutter feed
						push!(exe_log[r_ind+1], "[$(params.ident); r_ind = $r_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")			
						println("[$(params.ident); r_ind = $r_ind] iter $i of $(params.maxiter), curr_dist = $dist; ")
					end 
				end
				return g_r
			end
			# fixing initial g 
			time_iteration_fix_r_given_g(r) = time_iteration_fix_r(r, aiy_g_init);

		# Computing
			push!(exe_log[1], "[Aiyagari; $(params.ident)] starting time iteration; ")
			println("[Aiyagari; $(params.ident)] starting time iteration; ")
			@time aiy_container = pmap(time_iteration_fix_r_given_g, params.r_grid)		# pmap implements parallel computation
			for r_ind in 1:params.r_size
				aiy_g_star[:, r_ind, :] = aiy_container[r_ind];
			end 

		# Storing
			push!(exe_log[1], "[Aiyagari; $(params.ident)] storing; ")
			println("[Aiyagari; $(params.ident)] storing; ")
			save(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_policy.jld2"), 
				"aiy_g_star", aiy_g_star, "aiy_g_init", aiy_g_init, "exe_log", exe_log)
		
		println("[Aiyagari; $(params.ident)] PFI done & stored.")
	end