# Constant-time free-tree generation of Wright--Richmond--Odlyzko--McKay.
# This compact iterator is adapted from the MIT-licensed GKMtools.jl
# implementation by Csaba Schneider and Giosuè Muratore.
"""
    _FreeTreeIterator(n)

Enumerate one preorder depth sequence for every unrooted tree on `n` vertices
using the Wright–Richmond–Odlyzko–McKay constant-delay algorithm. Each yielded
vector is owned by the caller; the mutable successor state never aliases it.
Depth one is the root, and edge parents are the last preceding vertices one
level closer to the root. See the retained GKMtools notice for provenance.
"""
struct _FreeTreeIterator
  n::Int
  function _FreeTreeIterator(n::Integer)
    n > 0 || throw(ArgumentError("a free tree must have at least one vertex"))
    new(Int(n))
  end
end

Base.eltype(::Type{_FreeTreeIterator}) = Vector{Int}
Base.IteratorSize(::Type{_FreeTreeIterator}) = Base.SizeUnknown()

"""
    _FreeTreeState

Mutable integer bookkeeping for the WROM successor algorithm. The published
indices `p`, `q`, `h1`, `h2`, and `r` are retained to keep its case analysis
auditable against the original algorithm. `widths` is its auxiliary work array,
not a parent array. `center == typemax(Int)` represents the absent boundary
previously encoded by floating-point infinity. A zero `q` terminates iteration.
"""
mutable struct _FreeTreeState
  q::Int
  levels::Vector{Int}
  widths::Vector{Int}
  p::Int
  h1::Int
  h2::Int
  center::Int
  r::Int
end

function Base.iterate(iterator::_FreeTreeIterator)
  n = iterator.n
  if n < 4
    levels = Int[min(2, i) for i in 1:n]
    return copy(levels), _FreeTreeState(0, levels, Int[], 0, 0, 0, 0, 0)
  end
  r = n ÷ 2 + 1
  levels = Int[1:r; 2:(n - r + 1)]
  widths = Int[0:(r - 1); 1:(n - 1)]
  center = isodd(n) ? typemax(Int) : n + 1
  copy(levels), _FreeTreeState(n - 1, levels, widths, n, r, n, center, r)
end

function Base.iterate(iterator::_FreeTreeIterator, state::_FreeTreeState)
  iszero(state.q) && return nothing
  n = iterator.n
  q, levels, widths = state.q, state.levels, state.widths
  p, h1, h2, center, r = state.p, state.h1, state.h2, state.center, state.r
  fixit = false
  need_r = false
  need_center = false
  need_h2 = false

  if center == n + 1 ||
    (
    p == h2 &&
    (
      (levels[h1] == levels[h2] + 1 && n - h2 > r - h1) ||
      (levels[h1] == levels[h2] && n - h2 + 1 < r - h1)
    )
  )
    if levels[r] > 3
      p = r
      q = widths[r]
      h1 == r && (h1 -= 1)
      fixit = true
    else
      p = r
      r -= 1
      q = 2
    end
  end

  p <= h1 && (h1 = p - 1)
  if p <= r
    need_r = true
  elseif p <= h2
    need_h2 = true
  elseif levels[h2] == levels[h1] - 1 && n - h2 == r - h1
    p <= center && (need_center = true)
  else
    center = typemax(Int)
  end

  old_p = p
  delta = q - p
  old_level_q = levels[q]
  old_width_q = widths[q]
  p = typemax(Int)
  # Copy the successor suffix from an earlier portion of the level sequence.
  # Crossing a branch boundary invalidates the corresponding extremal indices;
  # the flags below postpone their repair until that boundary is reached.
  for i in old_p:n
    levels[i] = levels[i + delta]
    if levels[i] == 2
      widths[i] = 1
    else
      p = i
      q = levels[i] == old_level_q ? old_width_q : widths[i + delta] - delta
      widths[i] = q
    end
    if need_r && levels[i] == 2
      need_r = false
      need_h2 = true
      r = i - 1
    end
    if need_h2 && levels[i] <= levels[i - 1] && i > r + 1
      need_h2 = false
      h2 = i - 1
      if levels[h2] == levels[h1] - 1 && n - h2 == r - h1
        need_center = true
      else
        center = typemax(Int)
      end
    end
    if need_center
      if levels[i] != levels[h1 - h2 + i] - 1
        need_center = false
        center = i
      else
        center = i + 1
      end
    end
  end

  if fixit
    r = n - h1 + 1
    for i in (r + 1):n
      levels[i] = i - r + 1
      widths[i] = i - 1
    end
    widths[r + 1] = 1
    h2 = n
    p = n
    q = p - 1
    center = typemax(Int)
  else
    if p == typemax(Int)
      p = levels[old_p - 1] != 2 ? old_p - 1 : old_p - 2
      q = widths[p]
    end
    if need_h2
      h2 = n
      center = levels[h2] == levels[h1] - 1 && h1 == r ? n + 1 : typemax(Int)
    end
  end
  state.q, state.p, state.h1, state.h2, state.center, state.r = q, p, h1, h2, center, r
  copy(levels), state
end
