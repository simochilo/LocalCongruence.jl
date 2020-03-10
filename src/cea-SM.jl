"""
	vertCongruence(V::Lar.Array{Float64,2}; err=1e-6)

Evaluates the Vertex Congruence for 3D-points ``V ∈ ℳ(3,n)``.

The function determines the points of ``V`` closer than ``err`` and
builds a new Vertex Set made of the representative of each point cluster.

The method returns:
 - the new Vertex Set
 - a map that, for every new vertex, identifies the old vertices it is made of
"""
#function vertCongruence(
#		V::Lar.Array{Float64,2};
#		err=1e-6
#	)::Tuple{Lar.Points, Array{Array{Int,1},1}}
function vertCongruence(V; err=1e-6)

	Vcls    = []
	let nidx = 1,
		visited = [],
		kdtree  = NearestNeighbors.KDTree(V);
		for vidx = 1 : size(V, 2)  if !(vidx in visited)
			nearvs = NearestNeighbors.inrange(kdtree, V[:, vidx], err)
			push!(Vcls, nearvs)
			push!(visited, nearvs...)
			nidx += 1
		end  end
	end

	V = hcat([sum(V[:, cl], dims=2)/length(cl) for cl in Vcls]...)
	return V, Vcls
end


"""
	cellCongruence(
		cop::Lar.ChainOp,
		lo_cls::Array{Array{Int,1},1},
		lo_sign::Array{Array{Int8,1},1};
		imp = false,
		d = 0
	)

Evaluates the Cell Congruence for a Cochain ``cop`` with classes ``lo_cls``.

The function determines the new Cochain Operator built from ``cop`` where
the lower order cells are merged according to ``lo_cls`` map.
``lo_sign`` specifies whether a cell must be considered in reverse order.

If optional paramenter ``imp`` is set to ``true`` then FP imprecisions
are taken into account in the sense that lower order cells may have collided.
The parameter ``d`` represent then the order of the cell (that also is the
least number of lower order cells a current order cell is made of).

The method returns:
 - the new Cochain Operator
 - a map that, for every new cell, identifies the old cells it is made of
 - a map that, for every new cell, specify if old cells have changed ordering.
"""
function cellCongruence(cop, lo_cls, lo_sign; imp = false, d = 0)
#function cellCongruence(
#		cop::Lar.ChainOp,
#		lo_cls::Array{Array{Int,1},1},
#		lo_sign::Array{Array{Int8,1},1};
#		imp = false,
#		d = 0
#	)::Tuple{
#		Lar.ChainOp,
#		Array{Array{Int,1},1},
#		Array{Array{Int8,1},1}
#	}

	# Collide columns
	copCols = []
	for i = 1 : length(lo_cls)
		col = sum([
			cop[:, lo_cls[i][j]] .* lo_sign[i][j]
			for j = 1 : length(lo_cls[i])
		])
#		if imp
#			#TODO remove zeros
#			if length(col.nzind) > d  push!(copCols, col)  end
#		end
		push!(copCols, col)
	end

	# Retrieve Matrix Form and extract rows
#	[x for i = 1 : length(lo_cls) if (x = sum(cop[:, lo_cls[i]], dims = 2))>d]
#	rows = [[copCols[col][row] for row = 1 : cop.m] for col = 1:length(lo_cls)]
	cop = hcat(copCols...)
	rows = [cop[row, :] for row = 1 : cop.m]

	# Adjustinfg signs such that the first value always is `-1`
	sign = ones(Int8, cop.m)
	for row = 1 : cop.m  if rows[row].nzval[1] > 0
		rows[row] = -rows[row]
		sign[row] = -1
	end  end

	# Sort rows in order to collide them quickly
	rows_ord = sortperm(rows)
	nrows = unique(rows[rows_ord])
	ho_cls = [Array{Int,1}() for i = 1 : length(nrows)]
	nidx = 1

	# Collide cells with same structure
	for cidx = 1 : cop.m
		if rows[rows_ord[cidx]] != nrows[nidx]  nidx = nidx + 1;  end
		push!(ho_cls[nidx], rows_ord[cidx])
	end

	# Reshaping Sign and Cochain
	ho_sign = [[sign[el] for el in cell] for cell in ho_cls]
#	cop = convert(Lar.ChainOp, hcat(nrows...)')

	return hcat(nrows...)', ho_cls, ho_sign
end

"""
	chainCongruence(G::Lar.Points, T::Array{Lar.ChainOp}; imp=false, err=1e-6)

Perform the Geometry ``G`` congruence and coherently reshape the Topology ``T``
"""
function chainCongruence(G, T; imp=false, err=1e-6)
#function chainCongruence(
#		G::Lar.Points,
#		T::Array{Lar.ChainOp};
#		imp=false,
#		err=1e-6
#	)::Tuple{Lar.Points, Array{Lar.ChainOp}}

	# Perform the Geometry Congruence
	G, cls = vertCongruence(G; err=err)
	# Build default sign
	sign = [ones(Int8, length(cl)) for cl in cls]
	# Update the Topology coherently
	Tn = Array{Lar.ChainOp, 1}(undef, length(T))
	for d = 1 : length(T)
		Tn[d], cls, sign = cellCongruence(T[d], cls, sign; imp=imp, d=d)
	end
	# Euler Characteristic
	# @show size(G,2) - size(T[1],1) + size(T[2],1)
	return G, Tn
end
