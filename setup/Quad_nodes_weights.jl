### Description:
	# this file generates standard Gauss-Hermite and Simpson
	# nodes and weights and stores them 
	

### code
	function Quad_nodes_weights(params)
		# Nodes and weights for expectation integration: r 
			n_simps = range(params.s_low, params.s_high, length = params.n_nodes_simps);
			w_simps = [if i==1 || i==params.n_nodes_simps
								 step(n_simps)/3 
							 elseif isodd(i)
								 step(n_simps) * 2/3
							 else 
								 step(n_simps) * 4/3
							 end 
							for i in 1:params.n_nodes_simps];
		save(string("IRRCE_solution/internal_data/[", params.ident, "] Quad_nodes.jld2"), 
			"n_simps", n_simps, "w_simps", w_simps)
		println("[$(params.ident)] Quadrature nodes and weights stored.")
	end

	