"""
    _ExactSum(T)

Streaming balanced sum of exact scalars of type `T`. Binary digit `i` of `count`
indicates whether `partials[i]` holds a block of `2^(i-1)` nonzero summands.
Adding a term merges equal-size blocks like a binary carry. Only logarithmically
many partial sums are retained.

For rational localization, sequential addition repeatedly takes a gcd against
the large denominator of the entire prefix. Balanced blocks keep most gcds
between much smaller operands without changing the exact result. This is an
arithmetic optimization, with no tolerance or approximate cancellation.
"""
mutable struct _ExactSum{T}
  partials::Vector{T}
  count::UInt
end

_ExactSum(::Type{T}) where {T} = _ExactSum{T}(T[], zero(UInt))

"""Add an exact scalar to an owned balanced accumulator, skipping zero contributions."""
function _add_exact!(accumulator::_ExactSum{T}, value::T) where {T}
  iszero(value) && return accumulator
  next_count = Base.Checked.checked_add(accumulator.count, one(UInt))
  count = accumulator.count
  position = 1
  while isodd(count)
    value += accumulator.partials[position]
    count >>= 1
    position += 1
  end
  if position > length(accumulator.partials)
    push!(accumulator.partials, value)
  else
    accumulator.partials[position] = value
  end
  accumulator.count = next_count
  accumulator
end

"""Read the exact sum without consuming its blocks; further additions remain valid."""
function _sum_value(accumulator::_ExactSum{T}) where {T}
  total = zero(T)
  count = accumulator.count
  position = 1
  while !iszero(count)
    isodd(count) && (total += accumulator.partials[position])
    count >>= 1
    position += 1
  end
  total
end
