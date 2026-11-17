"""Checked modular scalar; the prime is validated once by `_cert_prime_type`."""
struct _CertPrime{P} <: Number
  value::Int64
  _CertPrime{P}(x::Integer) where {P} = new{P}(Int64(mod(x, P)))
end
_CertPrime{P}(x::Rational) where {P} =
  _CertPrime{P}(numerator(x)) / _CertPrime{P}(denominator(x))
_CertPrime{P}(x::Nemo.QQFieldElem) where {P} = _CertPrime{P}(Rational{BigInt}(x))
_CertPrime{P}(x::_CertPrime{P}) where {P} = x
Base.zero(::Type{_CertPrime{P}}) where {P} = _CertPrime{P}(0)
Base.one(::Type{_CertPrime{P}}) where {P} = _CertPrime{P}(1)
Base.zero(x::_CertPrime) = zero(typeof(x))
Base.one(x::_CertPrime) = one(typeof(x))
Base.iszero(x::_CertPrime) = iszero(x.value)
Base.isfinite(::_CertPrime) = true
Base.:+(a::_CertPrime{P}, b::_CertPrime{P}) where {P} = _CertPrime{P}(a.value + b.value)
Base.:-(a::_CertPrime{P}, b::_CertPrime{P}) where {P} = _CertPrime{P}(a.value - b.value)
Base.:-(a::_CertPrime{P}) where {P} = _CertPrime{P}(-a.value)
Base.:*(a::_CertPrime{P}, b::_CertPrime{P}) where {P} = _CertPrime{P}(a.value * b.value)
Base.inv(a::_CertPrime{P}) where {P} =
  iszero(a) ? throw(_LocalizationPole()) : _CertPrime{P}(invmod(a.value, P))
Base.:/(a::_CertPrime{P}, b::_CertPrime{P}) where {P} = a * inv(b)
Base.:(==)(a::_CertPrime{P}, b::_CertPrime{P}) where {P} = a.value == b.value
Base.convert(::Type{_CertPrime{P}}, x::Number) where {P} = _CertPrime{P}(x)
Base.promote_rule(::Type{_CertPrime{P}}, ::Type{I}) where {P,I<:Integer} = _CertPrime{P}
Base.show(io::IO, a::_CertPrime) = print(io, a.value)

"""Convert exact specialization constants without transporting FLINT objects to workers."""
function _cert_convert_evaluation(evaluation, ::Type{T}) where {T}
  converted = Any[]
  for field in fieldnames(typeof(evaluation))
    x = getfield(evaluation, field)
    value = if field == :word_values
      [T.(v) for v in x]
    elseif field == :billey
      Dict{Tuple{Int,Int},T}()
    elseif field == :edge_factor
      Dict{NTuple{3,Int},T}()
    elseif field == :inverse_flags
      Dict{NTuple{3,Int},NTuple{2,T}}()
    elseif field == :vertex_factor
      Dict{NTuple{5,Int},NTuple{2,T}}()
    elseif x isa AbstractArray{Nemo.QQFieldElem}
      T.(x)
    else
      x
    end
    push!(converted, value)
  end
  _GWEvaluation(converted...)
end

"""Reduce preflighted geometry constants, leaving unevaluated kernels initialized to zero."""
function _cert_convert_geometry(geometry, ::Type{T}) where {T}
  converted = Any[]
  for field in fieldnames(typeof(geometry))
    x = getfield(geometry, field)
    value = if field == :kernels
      result = zeros(T, size(x))
      for index in findall(geometry.kernel_known)
        result[index] = T(x[index])
      end
      result
    elseif x isa AbstractArray{Nemo.QQFieldElem}
      T.(x)
    else
      x
    end
    push!(converted, value)
  end
  _BranchGeometry(converted...)
end

"""A fixed-size rigorous absolute upper bound and optional reciprocal of an exact constant."""
struct _CertBound <: Number
  mantissa::UInt64
  exponent::Int64
  inverse_mantissa::UInt64
  inverse_exponent::Int64
end

"""Normalize a positive binary mantissa with upward integer rounding."""
function _cert_bound_normalize(m::UInt64, e::Int64)
  m == 0 && return UInt64(0), Int64(0)
  bits = 64 - leading_zeros(m)
  if bits > 32
    shift = bits - 32
    m = (m >> shift) + UInt64((m & ((UInt64(1) << shift) - 1)) != 0)
    e = Base.Checked.checked_add(e, shift)
    if m == (UInt64(1) << 32)
      m >>= 1
      e = Base.Checked.checked_add(e, 1)
    end
  elseif bits < 32
    shift = 32 - bits
    m <<= shift
    e = Base.Checked.checked_sub(e, shift)
  end
  m, e
end

"""Enclose an exact rational in a 32-bit-mantissa binary upper bound."""
function _cert_bound_rational_upper(value::Rational)
  value = abs(value)
  iszero(value) && return UInt64(0), Int64(0)
  n, d = numerator(value), denominator(value)
  e = Int64(ndigits(n; base=2) - ndigits(d; base=2) - 31)
  m = e >= 0 ? cld(n, d << e) : cld(n << (-e), d)
  _cert_bound_normalize(UInt64(m), e)
end

function _CertBound(value::Rational)
  m, e = _cert_bound_rational_upper(value)
  mi, ei = iszero(value) ? (UInt64(0), Int64(0)) : _cert_bound_rational_upper(inv(value))
  _CertBound(m, e, mi, ei)
end
_CertBound(value::Integer) = _CertBound(big(value)//1)
_CertBound(value::Nemo.QQFieldElem) = _CertBound(Rational{BigInt}(value))
_CertBound(value::_CertBound) = value
Base.zero(::Type{_CertBound}) = _CertBound(UInt64(0), Int64(0), UInt64(0), Int64(0))
Base.one(::Type{_CertBound}) =
  _CertBound(UInt64(1) << 31, Int64(-31), UInt64(1) << 31, Int64(-31))
Base.zero(::_CertBound) = zero(_CertBound)
Base.one(::_CertBound) = one(_CertBound)
Base.iszero(x::_CertBound) = x.mantissa == 0
Base.isfinite(::_CertBound) = true
Base.:-(x::_CertBound) = x

function Base.:+(a::_CertBound, b::_CertBound)
  iszero(a) && return b
  iszero(b) && return a
  a.exponent < b.exponent && ((a, b) = (b, a))
  shift = Base.Checked.checked_sub(a.exponent, b.exponent)
  smaller = if shift >= 64
    UInt64(1)
  elseif shift == 0
    b.mantissa
  else
    (b.mantissa >> shift) + UInt64((b.mantissa & ((UInt64(1) << shift) - 1)) != 0)
  end
  m, e = _cert_bound_normalize(a.mantissa + smaller, a.exponent)
  _CertBound(m, e, UInt64(0), Int64(0))
end
Base.:-(a::_CertBound, b::_CertBound) = a + b

function Base.:*(a::_CertBound, b::_CertBound)
  m, e = _cert_bound_normalize(a.mantissa * b.mantissa,
    Base.Checked.checked_add(a.exponent, b.exponent))
  _CertBound(m, e, UInt64(0), Int64(0))
end
function Base.:/(a::_CertBound, b::_CertBound)
  b.inverse_mantissa != 0 || error("majorant division requires a fixed nonzero constant")
  a * _CertBound(b.inverse_mantissa, b.inverse_exponent, UInt64(0), Int64(0))
end
Base.convert(::Type{_CertBound}, x::Number) = _CertBound(x)
Base.promote_rule(::Type{_CertBound}, ::Type{I}) where {I<:Integer} = _CertBound

"""Precomputed rigorous reciprocal bounds for the common degree divisors."""
const _CERT_RECIPROCALS = [_CertBound(big(1)//i) for i in 1:20]
function Base.:*(a::_CertBound, b::Int)
  0 <= b <= 20 || return a * _CertBound(big(b))
  m, e = _cert_bound_normalize(a.mantissa * UInt64(b), a.exponent)
  _CertBound(m, e, UInt64(0), Int64(0))
end
Base.:*(a::Int, b::_CertBound) = b * a
Base.:/(a::_CertBound, b::Int) =
  1 <= b <= 20 ? a * _CERT_RECIPROCALS[b] : a / _CertBound(big(b))

"""Convert a binary upper bound to its exact integer ceiling."""
function _cert_bound_integer_ceiling(value::_CertBound)
  if value.exponent >= 0
    big(value.mantissa) << value.exponent
  else
    cld(big(value.mantissa), big(1) << (-value.exponent))
  end
end
