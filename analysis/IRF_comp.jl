### Description:
	# this file produces the results from various comparisons 
	# between IRRCE and the Aiyagari equilibrium, stores key 
	# data, and outputs selected results as graphics and tables

### code 

	# visual equilibrium results
		function IRF_graph(params; path::String = "C:/Users/ulric/Desktop/Sonstige Uniarbeiten/MSc/Master thesis/Code/RoschitschUlrich_MasterThesis/fifth version/analysis/output_data", 
						   horizon = 300)
			# Loading & Setup
				aiy_aggr_state_IRF::Array{Float64,2} = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium_IRF.jld2"))["aiy_aggr_state_IRF"]
				aiy_r_star::Float64 		 	 	 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_r_star"]
				aiy_K::Float64 				 	 	 = load(string("aiyagari_solution/internal_data/[", params.ident, "] aiy_equilibrium.jld2"))["aiy_K"]
				aggr_state_IRF::Array{Float64,2} 	 = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium_IRF.jld2"))["aggr_state_IRF"]
				aggr_state::Array{Float64,2} 	 	 = load(string("IRRCE_solution/internal_data/[", params.ident, "] equilibrium.jld2"))["aggr_state"]
				Output 						 	 	 = params.A * aggr_state[1,end]^(params.α)
				aiy_Output					 	 	 = params.A * aiy_K^(params.α)
				Output_IRF					 	 	 = aggr_state_IRF[5,:] .* aggr_state_IRF[1,:].^(params.α)
				aiy_Output_IRF				 	 	 = aiy_aggr_state_IRF[4,:] .* aiy_aggr_state_IRF[1,:].^(params.α)

			# Plots
				plt_upleft    = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="north east",
						title = raw"TFP shock ($\%$ dev. from sS)"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:horizon, aggr_state_IRF[5,1:horizon] .- 1.0)
					),
					LegendEntry(raw"TFP $A_{t}$")
				)
				
				plt_upright    = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="north east",
						title = raw"Short rate ($\%$ dev. from sS)"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:horizon, aggr_state_IRF[2,1:horizon]./aggr_state[2,end] .- 1.0)
					),
					LegendEntry(raw"IRRCE $r_{t}$"),
					Plot(
						{color = "red"},
						Coordinates(1:horizon, aiy_aggr_state_IRF[2,1:horizon]./aiy_r_star .- 1.0)
						),
					LegendEntry(raw"Aiyagari $r^{*}_{t}$"),
					Plot(
						{
							no_marks,
							color = "blue",
							dashed
						},
						Coordinates(1:horizon, (exp.(aggr_state_IRF[4,1:horizon]) .-1.0)./(exp(aggr_state[4,end]) -1.0) .- 1.0)
					),
					LegendEntry(raw"Mean beliefs $\exp(d_{t-1})-1$"),
				)

				plt_downleft   = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="north east",
						title = raw"Aggr. capital ($\%$ dev. from sS)"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:horizon, aggr_state_IRF[1,1:horizon]./aggr_state[1,end] .- 1.0)
					),
					LegendEntry(raw"IRRCE $K_t$"),
					Plot(
						{color = "red"},
						Coordinates(1:horizon, aiy_aggr_state_IRF[1,1:horizon]./aiy_K .- 1.0)
						),
					LegendEntry(raw"Aiyagari $K^{*}_{t}$"),
				)

				plt_downright  = @pgf Axis(
					{
						xmajorgrids,
						ymajorgrids,
						major_grid_style = "{color=black!5!white}",
						xlabel = L"t",
						legend_pos="north east",
						title = raw"Output ($\%$ dev. from sS)"
					},
					Plot(
						{
							no_marks,
							color = "blue"
						},
						Coordinates(1:horizon, Output_IRF[1:horizon]./Output .-1.0)
					),
					LegendEntry(raw"IRRCE"),
					Plot(
						{color = "red"},
						Coordinates(1:horizon, aiy_Output_IRF[1:horizon]./aiy_Output .-1.0)
						),
					LegendEntry(raw"Aiyagari"),
				)

				IRFgraph = @pgf GroupPlot(
					{ group_style = { group_size="2 by 2", vertical_sep="60pt"},
					no_markers,
					legend_style={font=raw"\tiny"},
					title_style={font=raw"\footnotesize"},
					},
					plt_upleft, plt_upright, plt_downleft, plt_downright)
			
			# Storing 
				pgfsave(string(path, "/", params.ident, "_IRFgraph.tikz"), IRFgraph)
		end