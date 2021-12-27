module XICOR
	using Base: length, sum, abs, zeros, size, diff
	using StatsBase: sortperm, tiedrank
	export xicor

	function xicor(x::Vector{T},y::Vector{T}) where {T<:Real}
		n = length(x);
		if n != length(y)
			error("Vectors must be the same length")
		end
		xord = sortperm(x);
		r_i = tiedrank(y[xord]);
		l_i = tiedrank(y[xord],rev=true);
		RS =  n*sum(abs,diff(r_i));
		LS =  2*l_i'*(n .- l_i);
		# return 1 - n*sum(abs,r_i[2:end] .- r_i[1:end-1])/(2*sum(l_i.*(n .- l_i))),r_i,l_i;
		return 1 - RS/LS;
	end
	

	xicor(X::Matrix{T},y::Vector{T}) where {T<:Real} = [xicor(X[:,col],y) for col in 1:size(X,2)];
	xicor(x::Vector{T},Y::Matrix{T}) where {T<:Real} = [xicor(x,Y[:,col]) for col in 1:size(Y,2)];
	xicor(X::Matrix{T},Y::Matrix{T}) where {T<:Real} = [xicor(X[:,colX],Y[:,colY]) for colX in 1:size(X,2),colY in 1:size(Y,2)];

	function xicor(X::Matrix{T}) where {T<:Real}
		nRows,nColumns = size(X);
		ximat = zeros(nColumns,nColumns);

		for jColumn in 1:nColumns
			for iColumn in 1:jColumn
				xloop = X[:,iColumn];
				yloop = X[:,jColumn];
				ximat[iColumn,jColumn]= xicor(xloop,yloop);
				ximat[jColumn,iColumn]= xicor(yloop,xloop);
			end
		end

		return ximat
	end
end # module
