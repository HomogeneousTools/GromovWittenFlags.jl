"""
    _independent_indices(matrix)

Select the first linearly independent rows of a finite exact rational matrix,
returning their row indices in original order. Maintain a normalized row-echelon
basis incrementally; each input row is reduced once. Pivot rows are kept in column
order even when a later input row introduces an earlier pivot column.

Each row is converted directly into an owned `Rational{BigInt}` scratch vector;
the input and its views are never mutated. For a symmetric intersection matrix,
the selected rows also select a nonsingular principal submatrix: their
coordinate classes form a basis modulo the pairing's radical.
"""
function _independent_indices(matrix::AbstractMatrix{<:Rational})
  selected = Int[]
  pivots = Int[]
  basis = Vector{_BigRational}[]
  n_columns = size(matrix, 2)
  for i in axes(matrix, 1)
    row = collect(_BigRational, @view matrix[i, :])
    for (pivot, previous) in zip(pivots, basis)
      coefficient = row[pivot]
      iszero(coefficient) && continue
      row[pivot] = zero(_BigRational)
      for column in (pivot + 1):n_columns
        row[column] -= coefficient * previous[column]
      end
    end
    pivot = findfirst(!iszero, row)
    pivot === nothing && continue
    coefficient = row[pivot]
    row[pivot] = one(_BigRational)
    for column in (pivot + 1):n_columns
      row[column] /= coefficient
    end
    position = searchsortedfirst(pivots, pivot)
    insert!(pivots, position, pivot)
    insert!(basis, position, row)
    push!(selected, i)
    length(selected) == n_columns && break
  end
  selected
end

"""Compute exact rational rank using the same incremental elimination as row selection."""
_exact_rank(matrix::AbstractMatrix{<:Rational}) = length(_independent_indices(matrix))
