### Description:
	# this file contains a function that writes a (a,z)-transition
	# matrix from a z-transition matrix and some a′(a,z)-policy function 
	

### code
	function trans_mat(Π, a_distnodes, z_distnodes, a_policy)
		# note: all columns of Π should sum to 1
		a_len = length(a_distnodes); z_len = length(z_distnodes)
		Γ = zeros(Float64, a_len*z_len, a_len*z_len)
		for z_ind in 1:z_len
			for a_ind in 1:a_len
				a′ = a_policy(a_distnodes[a_ind], z_distnodes[z_ind])
				a_indlow = searchsortedlast(a_distnodes, a′)
					# 'searchsortedlast' returns the index of the last value in a_distnodes less than or equal to a′;
					#  Returns 0 if a′ is less than all values in a_distnodes; link: http://www.jlhub.com/julia/manual/en/function/searchsortedlast 
				if a_indlow == 0 || a_indlow == a_len
					Γ[((1:z_len).*a_len .- a_len .+ min(max(a_indlow,1), a_len)), ((z_ind - 1)*a_len + a_ind)] = Π[:,z_ind]
				else 
					ω_low = (a_distnodes[a_indlow+1] - a′)/(a_distnodes[a_indlow+1] - a_distnodes[a_indlow])
					Γ[((1:z_len).*a_len .- a_len .+ a_indlow), ((z_ind - 1)*a_len + a_ind)] 		= ω_low.*Π[:,z_ind]
					Γ[((1:z_len).*a_len .- a_len .+ a_indlow .+ 1), ((z_ind - 1)*a_len + a_ind)] 	= (1.0 - ω_low).*Π[:,z_ind]
				end 
			end 
		end 
		return Γ
	end