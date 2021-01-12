### Description:
	# this file produces the results from various comparisons 
	# between IRRCE and the Aiyagari equilibrium, stores key 
	# data, and outputs selected results as graphics and tables


### code 

	# visual equilibrium results
		function equilibrium_graph(params; path::String = "C:/Users/ulric/Desktop/Sonstige Uniarbeiten/MSc/Master thesis/Code/RoschitschUlrich_MasterThesis/fifth version/analysis/output_data")
			# Loading & Setup
				aiy_r_star::Float64 		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_r_star"]
				aiy_K::Float64 				 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_K"]
				aiy_dist::Array{Float64,1} 	 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_dist"]
				aggr_state::Array{Float64,2} = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["aggr_state"]
				dist::Array{Float64,2} 		 = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["dist"]
				Output 						 = params.A .* aggr_state[1,:].^(params.α)
				aiy_Output					 = params.A .* aiy_K.^(params.α)

			# Plots
				plt_upleft    = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="north east",
						title = "Short rate"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:params.T, aggr_state[2,:])
					),
					LegendEntry(raw"IRRCE $r_{t}$"),
					Plot(
						{color = "red"},
						Coordinates(1:params.T, repeat([aiy_r_star], params.T))
						),
					LegendEntry(raw"Aiyagari $r^{*}$"),
					Plot(
						{
							no_marks,
							color = "blue",
							dashed
						},
						Coordinates(1:params.T, exp.(aggr_state[4,:]) .-1.0)
					),
					LegendEntry(raw"Mean beliefs $\exp(d_{t-1})-1$"),
				)

				plt_upright   = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="south east",
						title = "Aggr. capital"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:params.T, aggr_state[1,:])
					),
					LegendEntry(raw"IRRCE $K_t$"),
					Plot(
						{color = "red"},
						Coordinates(1:params.T, repeat([aiy_K], params.T))
						),
					LegendEntry(raw"Aiyagari $K^{*}$"),
				)

				plt_downleft  = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="south east",
						title = "Output"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:params.T, Output)
					),
					LegendEntry(raw"IRRCE"),
					Plot(
						{color = "red"},
						Coordinates(1:params.T, repeat([aiy_Output], params.T))
						),
					LegendEntry(raw"Aiyagari"),
				)

				# arranging distribution arrays 
				aiy_dens		 = sum(reshape(aiy_dist, params.a_distnodes_size, params.z_size), dims = 2)[:,1]
				IRRCE_2_dens	 = sum(reshape(dist[:,2], params.a_distnodes_size, params.z_size), dims = 2)[:,1]
				IRRCE_6_dens	 = sum(reshape(dist[:,6], params.a_distnodes_size, params.z_size), dims = 2)[:,1]
				IRRCE_51_dens	 = sum(reshape(dist[:,51], params.a_distnodes_size, params.z_size), dims = 2)[:,1]
				IRRCE_final_dens = sum(reshape(dist[:,end], params.a_distnodes_size, params.z_size), dims = 2)[:,1]
				upper_lim        = maximum(vcat(IRRCE_final_dens,aiy_dens)) + 1e-2
				plt_downright = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xmin = -0.5,
						xmax = 30,
						ymin = 0.0,
						ymax = upper_lim,
						ytick = 0:0.01:upper_lim,
						scaled_y_ticks="base 10:2",
						xtick = 0:5:30,
						xlabel = raw"assets $a_{t}^{i}$",
						legend_pos="north east",
						title = "Wealth Distributions"
					},
					PlotInc(
						{
							ybar,
							bar_width="0.01pt",
							color = "blue!50!white",
							solid
						},
						Coordinates([params.a_init], [1.0])
					),
					LegendEntry(raw"IRRCE, $t=0$"),
					Plot(
						{
							no_marks,
							color = "blue!50!white",
							dashed
						},
						Coordinates(params.a_distnodes, IRRCE_2_dens)
					),
					LegendEntry(raw"IRRCE, $t=1$"),
					Plot(
						{
							no_marks,
							color = "blue!50!white",
							loosely_dotted
						},
						Coordinates(params.a_distnodes, IRRCE_6_dens)
					),
					LegendEntry(raw"IRRCE, $t=5$"),
					Plot(
						{
							no_marks,
							color = "blue!50!white",
							dashdotted
						},
						Coordinates(params.a_distnodes, IRRCE_51_dens)
					),
					LegendEntry(raw"IRRCE, $t=50$"),
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(params.a_distnodes, IRRCE_final_dens)
					),
					LegendEntry(raw"IRRCE, $t=T$"),
					Plot(
						{
							no_marks,
							color = "red"
						},
						Coordinates(params.a_distnodes, aiy_dens)
					),
					LegendEntry(raw"Aiyagari"),
				)

				eqm_graph = @pgf GroupPlot(
					{ group_style = { group_size="2 by 2", vertical_sep="60pt"},
					no_markers,
					legend_style={font=raw"\tiny"},
					title_style={font=raw"\footnotesize"},
					},
					plt_upleft, plt_upright, plt_downleft, plt_downright)

			# Storing 
				pgfsave(string(path, "/", params.ident, "_eqmgraph.tikz"), eqm_graph)
		end

	# key statistics
		function key_eqm_stats(params; path::String = "C:/Users/ulric/Desktop/Sonstige Uniarbeiten/MSc/Master thesis/Code/RoschitschUlrich_MasterThesis/fifth version/analysis/output_data")
			# Loading & Setup
				aiy_r_star::Float64 		 		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_r_star"]
				aiy_K::Float64 				 		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_K"]
				aiy_dist::Array{Float64,1} 	 		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_dist"]
				aiy_g::Array{Float64,2} 	 		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_PF"]
				aiy_transmat::Array{Float64,2}		 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_transmat"]
				aggr_state::Array{Float64,2} 		 = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["aggr_state"]
				dist::Array{Float64,2} 		 		 = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["dist"]
				g_star::Array{Float64,4} 	 		 = load(string("IRRCE_solution/internal_data/[", params.ident, "] policy.jld2"))["g_star"]
				aggr_state_fromAiy::Array{Float64,2} = load(string("IRRCE_solution/internal_data/[", params.ident, "_fromAiy] equilibrium.jld2"))["aggr_state"]
				println("[Key eqm stats; $(params.ident)] loading done; ")
				aiy_ḡ						 		 = LinearInterpolation((params.a_grid, params.z_grid), aiy_g, extrapolation_bc = Line())
				ḡ_star						 		 = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid, params.d_grid), g_star, extrapolation_bc = Line())
				aiy_a								 = [(1+aiy_r_star)*a + params.w(aiy_r_star)*z - aiy_g[a_ind, z_ind]   
														for (a_ind, a) in enumerate(params.a_grid), (z_ind, z) in enumerate(params.z_grid)]
				aiy_ā								 = LinearInterpolation((params.a_grid, params.z_grid), aiy_a, extrapolation_bc = Line())
				a_star					 			 = [(1+r)*a + params.w(r)*z - g_star[a_ind, r_ind, z_ind, d_ind]   
														for (a_ind, a) in enumerate(params.a_grid), (r_ind, r) in enumerate(params.r_grid), 
															(z_ind, z) in enumerate(params.z_grid), (d_ind, d) in enumerate(params.d_grid)]
				ā_star 					 			 = LinearInterpolation((params.a_grid, params.r_grid, params.z_grid, params.d_grid), a_star, 
														extrapolation_bc = Line())
				distmat_az							 = reshape(dist[:,end], params.a_distnodes_size, params.z_size)
				aiy_distmat_az				 		 = reshape(aiy_dist[:,end], params.a_distnodes_size, params.z_size)

			# Quantities needed: r, K, Y, SR, c, a
				# IRRCE
				r = aggr_state[2,end]
				K = aggr_state[1,end]
				Y = params.A .* K.^(params.α);
				SR = mean([(aggr_state[1,(500+t)] - (1-params.δ)*aggr_state[1,(500+t-1)]) /
										(params.A * aggr_state[1,(500+t-1)]^(params.α))   for t in 1:(params.T - 500)])
				Eassets = dot(params.a_distnodes, sum(distmat_az, dims=2)[:,1])
				STDassets = dot((params.a_distnodes .- Eassets).^2, sum(distmat_az, dims=2)[:,1])^(1/2)
				Econs = sum([ḡ_star(a, aggr_state[2,end], z, aggr_state[4,end])*distmat_az[a_ind, z_ind]
							for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)])
				STDcons = sum([(ḡ_star(a, aggr_state[2,end], z, aggr_state[4,end]) - Econs)^2 *distmat_az[a_ind, z_ind]
							for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)])^(1/2)
				# Aiyagari
				aiy_r = aiy_r_star
				aiy_Y = params.A * aiy_K^(params.α)
				aiy_SR = (params.α*params.δ) / (params.δ + aiy_r)
				Eaiy_assets = dot(params.a_distnodes, sum(aiy_distmat_az, dims=2)[:,1])
				STDaiy_assets = dot((params.a_distnodes .- Eaiy_assets).^2, sum(aiy_distmat_az, dims=2)[:,1])^(1/2)
				Eaiy_cons = sum([aiy_ḡ(a, z)*aiy_distmat_az[a_ind, z_ind]
								for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)])
				STDaiy_cons = sum([(aiy_ḡ(a, z) - Eaiy_cons)^2 *aiy_distmat_az[a_ind, z_ind]
								for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)])^(1/2)

			# Producing table "key equilibrium statistics"
				latex_tabular(string(path, "/", params.ident, "_equilibrium_table.tex"),
					Tabular("l l l l l l l l l"),
					[[MultiColumn(9, :c, "Key Equilibrium Statistics")],
					Rule(:top),
					vcat("Aggregates", L"r_{t}", L"K_{t}", L"Y_{t}", L"SR_{t}", L"c_{t}^{i}", " ", L"a_{t}^{i}", " "),
					vcat(L"\vphantom{\HUGE A}", "(sl)", "(sl)", "(sl)", "(sl)", "(m)", raw"(sd)$\times 10^{2}$", "(m)", raw"(sd)$\times 10^{2}$"),
					Rule(:mid),
					vcat("IRRCE", string(format(100*r, precision = 2), L" \%"), format(K, precision = 2),  format(Y, precision = 2), 
						string(format(SR*100, precision = 2), L" \%"), format(Econs, precision = 2), format(STDcons*100, precision = 2), 
						format(Eassets, precision = 2), format(STDassets*100, precision = 0)),
					vcat("Aiyagari", string(format(aiy_r*100, precision = 2), L" \%"), format(aiy_K, precision = 2), format(aiy_Y, precision = 2), string(format(aiy_SR*100, precision = 2), L" \%"),
					format(Eaiy_cons, precision = 2), format(STDaiy_cons*100, precision = 2), format(Eaiy_assets, precision = 2), format(STDaiy_assets*100, precision = 0)),
					Rule(:bottom)]
				)
				println("[Key eqm stats; $(params.ident)] table 1/3 done; ")
				
			# distributional stats 
				Qs = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]	# displayed quantiles - maximally 8 (or else adjust table formatting)
				# IRRCE
				a_dist = DiscreteNonParametric(params.a_distnodes, sum(distmat_az, dims=2)[:,1])
				E = mean(a_dist)
				S = std(a_dist)
				skw = skewness(a_dist)	# Reference: https://en.wikipedia.org/wiki/Skewness
				krt = kurtosis(a_dist)	# Reference: https://en.wikipedia.org/wiki/Kurtosis
				gini = 1.0/(2*E) * sum([probs(a_dist)[a_i_ind] *probs(a_dist)[a_j_ind] *abs(a_i - a_j)
						for (a_i_ind, a_i) in enumerate(params.a_distnodes), (a_j_ind, a_j) in enumerate(params.a_distnodes)])
						# Reference: https://en.wikipedia.org/wiki/Gini_coefficient
				quants = quantile.(a_dist, Qs)
				wquants = [dot(params.a_distnodes[params.a_distnodes.≤quants[i]], probs(a_dist)[params.a_distnodes.≤quants[i]]) / 
							(E*sum(probs(a_dist)[params.a_distnodes.≤quants[i]]))   for i in 1:length(quants)]
				# Aiyagari
				aiy_a_dist = DiscreteNonParametric(params.a_distnodes, sum(aiy_distmat_az, dims=2)[:,1])
				aiy_E = mean(aiy_a_dist)
				aiy_S = std(aiy_a_dist)
				aiy_skw = skewness(aiy_a_dist)	# Reference: https://en.wikipedia.org/wiki/Skewness
				aiy_krt = kurtosis(aiy_a_dist)	# Reference: https://en.wikipedia.org/wiki/Kurtosis
				aiy_gini = 1.0/(2*aiy_E) * sum([probs(aiy_a_dist)[a_i_ind] *probs(aiy_a_dist)[a_j_ind] *abs(a_i - a_j)
							for (a_i_ind, a_i) in enumerate(params.a_distnodes), (a_j_ind, a_j) in enumerate(params.a_distnodes)])
							# Reference: https://en.wikipedia.org/wiki/Gini_coefficient
				aiy_quants = quantile.(aiy_a_dist, Qs)
				aiy_wquants = [dot(params.a_distnodes[params.a_distnodes.≤aiy_quants[i]], probs(aiy_a_dist)[params.a_distnodes.≤aiy_quants[i]]) / 
								(E*sum(probs(aiy_a_dist)[params.a_distnodes.≤aiy_quants[i]]))   for i in 1:length(aiy_quants)]

			# Producing table "Statistics On The Wealth Distribution"
				latex_tabular(string(path, "/", params.ident, "_distribution_table.tex"),
					Tabular("l l l l l l l l l"),
					[[MultiColumn(9, :c, "Statistics On The Wealth Distribution")],
					Rule(:top),
					vcat("Moments", "mean", "std. dev.", "skew.", "kurt.", "Gini", " ", " ", " "),
					Rule(:mid),
					vcat("IRRCE", format(E, precision = 2), format(S, precision = 2), format(skw, precision = 2), format(krt, precision = 2), format(gini, precision = 2), " ", " ", " "),
					vcat("Aiyagari", format(aiy_E, precision = 2), format(aiy_S, precision = 2), format(aiy_skw, precision = 2), format(aiy_krt, precision = 2), format(aiy_gini, precision = 2), " ", " ", " "),
					Rule(:mid),
					vcat("Quant.", string("Q", floor(Int64, Qs[1]*100)), string("Q", floor(Int64, Qs[2]*100)), string("Q", floor(Int64, Qs[3]*100)), 
						string("Q", floor(Int64, Qs[4]*100)), string("Q", floor(Int64, Qs[5]*100)), string("Q", floor(Int64, Qs[6]*100)), 
						string("Q", floor(Int64, Qs[7]*100)), string("Q", floor(Int64, Qs[8]*100))),
					Rule(:mid),
					vcat("IRRCE", format.(quants, precision = 2)),
					vcat("Aiyagari", format.(aiy_quants, precision = 2)),
					Rule(:mid),
					vcat(string("Wealth held ", L"(\%)"), string("wQ", floor(Int64, Qs[1]*100)), string("wQ", floor(Int64, Qs[2]*100)), string("wQ", floor(Int64, Qs[3]*100)), 
						string("wQ", floor(Int64, Qs[4]*100)), string("wQ", floor(Int64, Qs[5]*100)), string("wQ", floor(Int64, Qs[6]*100)), 
						string("wQ", floor(Int64, Qs[7]*100)), string("wQ", floor(Int64, Qs[8]*100))),
					Rule(:mid),
					vcat("IRRCE", format.(wquants, precision = 2)),
					vcat("Aiyagari", format.(aiy_wquants, precision = 2)),
					Rule(:bottom)]
				)
				println("[Key eqm stats; $(params.ident)] table 2/3 done; ")

			# Welfare function 
				function welfares(a_init)
					# IRRCE
					V = 0.0
					dist_container =  zeros(Float64, params.a_distnodes_size * params.z_size)
					for z_ind in 1:params.z_size		# implementing mass point at a_init in first distribution 
						dist_container[(z_ind-1)*params.a_distnodes_size + maximum(findall(params.a_distnodes .≤ a_init))] = params.z_statdist[z_ind]
					end
					for t in 1:params.T
						distmat = reshape(dist_container, params.a_distnodes_size, params.z_size)
						util_mat = [ḡ_star(a, aggr_state_fromAiy[2,t], z, aggr_state_fromAiy[4,t])^(1.0-params.σ) / 
									(1.0-params.σ) *distmat[a_ind, z_ind]
									for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)]
						V += params.β^(t-1) * sum(util_mat)
						ā_policy(a, z)  = ā_star(a, aggr_state_fromAiy[2,t], z, aggr_state_fromAiy[4,t])
						transmat = trans_mat(params.Π, params.a_distnodes, params.z_grid, ā_policy)
						dist_container = transmat * dist_container
					end
					# Aiyagari
					aiy_V = 0.0
					aiy_dist_container =  zeros(Float64, params.a_distnodes_size * params.z_size)
					for z_ind in 1:params.z_size		# implementing mass point at a_init in first distribution 
						aiy_dist_container[(z_ind-1)*params.a_distnodes_size + maximum(findall(params.a_distnodes .≤ a_init))] = params.z_statdist[z_ind]
					end
					for t in 1:params.T
						aiy_distmat = reshape(aiy_dist_container, params.a_distnodes_size, params.z_size)
						aiy_util_mat = [aiy_ḡ(a, z)^(1.0-params.σ) / (1.0-params.σ) *aiy_distmat[a_ind, z_ind]
										for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)]
						aiy_V += params.β^(t-1) * sum(aiy_util_mat)
						aiy_dist_container = aiy_transmat * aiy_dist_container
					end

					return ((V / aiy_V)^(1.0/(1.0-calib_1.σ)) -1.0)*100
				end
				# w/o transition
				function welfares_wo_trans(a_init)
					# IRRCE
					V = 0.0
					dist_container =  zeros(Float64, params.a_distnodes_size * params.z_size)
					for z_ind in 1:params.z_size		# implementing mass point at a_init in first distribution 
						dist_container[(z_ind-1)*params.a_distnodes_size + maximum(findall(params.a_distnodes .≤ a_init))] = params.z_statdist[z_ind]
					end
					ā_policy(a, z)  = ā_star(a, aggr_state[2,end], z, aggr_state[4,end])
					transmat = trans_mat(params.Π, params.a_distnodes, params.z_grid, ā_policy)
					for t in 1:params.T
						distmat = reshape(dist_container, params.a_distnodes_size, params.z_size)
						util_mat = [ḡ_star(a, aggr_state[2,end], z, aggr_state[4,end])^(1.0-params.σ) / 
									(1.0-params.σ) *distmat[a_ind, z_ind]
									for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)]
						V += params.β^(t-1) * sum(util_mat)
						dist_container = transmat * dist_container
					end
					# Aiyagari
					aiy_V = 0.0
					aiy_dist_container =  zeros(Float64, params.a_distnodes_size * params.z_size)
					for z_ind in 1:params.z_size		# implementing mass point at a_init in first distribution 
						aiy_dist_container[(z_ind-1)*params.a_distnodes_size + maximum(findall(params.a_distnodes .≤ a_init))] = params.z_statdist[z_ind]
					end
					for t in 1:params.T
						aiy_distmat = reshape(aiy_dist_container, params.a_distnodes_size, params.z_size)
						aiy_util_mat = [aiy_ḡ(a, z)^(1.0-params.σ) / (1.0-params.σ) *aiy_distmat[a_ind, z_ind]
										for (a_ind, a) in enumerate(params.a_distnodes), (z_ind, z) in enumerate(params.z_grid)]
						aiy_V += params.β^(t-1) * sum(aiy_util_mat)
						aiy_dist_container = aiy_transmat * aiy_dist_container
					end

					return ((V / aiy_V)^(1.0/(1.0-calib_1.σ)) -1.0)*100
				end

			# Welfare data 
				println("[Key eqm stats; $(params.ident)] computing welfare changes - this may take 2-5 min, get a coffee :) ; ")
				aiy_quants_welf 		 = [welfares(q)   for q in aiy_quants]
				aiy_quants_welf_wo_trans = [welfares_wo_trans(q)   for q in aiy_quants]

			# Producing table "Consumption-Equivalent Welfare Changes"
				latex_tabular(string(path, "/", params.ident, "_welfare_changes.tex"),
					Tabular("l l l l l l l l l"),
					[[MultiColumn(9, :c, "Consumption-Equivalent Welfare Changes")],
					Rule(:top),
					vcat("Quantile", string("Q", floor(Int64, Qs[1]*100)), string("Q", floor(Int64, Qs[2]*100)), string("Q", floor(Int64, Qs[3]*100)), 
						string("Q", floor(Int64, Qs[4]*100)), string("Q", floor(Int64, Qs[5]*100)), string("Q", floor(Int64, Qs[6]*100)), 
						string("Q", floor(Int64, Qs[7]*100)), string("Q", floor(Int64, Qs[8]*100))),
					Rule(:mid),
					vcat(string(L"\Delta" ,"Welf. ", L"(\%)"), format.(aiy_quants_welf, precision = 2)),
					vcat(string(L"\Delta" ,"Welf. w/o trans. ", L"(\%)"), format.(aiy_quants_welf_wo_trans, precision = 2)),
					Rule(:bottom)]
				)
				println("[Key eqm stats; $(params.ident)] table 3/3 done; ")
		end;