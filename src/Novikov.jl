"""
    NovikovPolynomial{R}()
    NovikovPolynomial{R}(constant)
    NovikovPolynomial{R}(terms::AbstractDict)

An exact sparse polynomial in `R` Novikov variables, with `Rational{BigInt}`
coefficients and nonnegative integer exponents. Construct from an integer or
finite rational constant, or a dictionary mapping exponent tuples/vectors to such
coefficients. Construction copies and normalizes the dictionary: equivalent
exponents are combined and zero coefficients are removed.

Arithmetic supports addition, subtraction, polynomial multiplication, and
exact scalar multiplication. Different variable counts belong to different
polynomial rings. Treat the stored `terms` dictionary as read-only; construct
a new polynomial to normalize an edited dictionary. Matrix entries returned
by [`small_quantum_multiplication`](@ref) have independent dictionaries.

# Examples
```jldoctest
julia> P = NovikovPolynomial{2};

julia> p = P(Dict((1, 0) => 2, (0, 1) => -1));

julia> evaluate(p, (2, 3))
1//1

julia> iszero(p - p)
true
```
"""
struct NovikovPolynomial{R}
  terms::Dict{NTuple{R,Int},_BigRational}

  function NovikovPolynomial{R}(source::AbstractDict) where {R}
    _check_novikov_rank(Val(R))
    terms = Dict{NTuple{R,Int},_BigRational}()
    for (exponent, coefficient) in source
      exponent isa Union{Tuple,AbstractVector} && length(exponent) == R ||
        throw(ArgumentError("expected $R exponents in a tuple or vector"))
      all(e -> e isa Integer && 0 <= e <= typemax(Int), exponent) ||
        throw(ArgumentError("Novikov exponents must be nonnegative machine integers"))
      coefficient isa Union{Integer,Rational} && isfinite(coefficient) ||
        throw(ArgumentError("Novikov coefficients must be finite integers or rationals"))
      coordinates = Tuple(exponent)
      key = ntuple(i -> Int(coordinates[i]), R)
      _accumulate_term!(terms, key, _BigRational(coefficient))
    end
    new{R}(terms)
  end

  # Internal ownership constructor. Callers must supply a fresh, normalized
  # dictionary; bypassing a second validation/copy is safe only under that
  # invariant. Public construction always goes through the method above.
  NovikovPolynomial{R}(
    terms::Dict{NTuple{R,Int},_BigRational}, ::Val{:canonical}
  ) where {R} = new{R}(terms)
end

"""Validate the variable count before constructing exponent tuples or typed dictionaries."""
function _check_novikov_rank(::Val{R}) where {R}
  R isa Int && R >= 0 ||
    throw(ArgumentError("the number of variables must be a nonnegative Int"))
  nothing
end

function NovikovPolynomial{R}() where {R}
  _check_novikov_rank(Val(R))
  NovikovPolynomial{R}(Dict{NTuple{R,Int},_BigRational}(), Val(:canonical))
end

function NovikovPolynomial{R}(constant::Union{Integer,Rational}) where {R}
  _check_novikov_rank(Val(R))
  NovikovPolynomial{R}(Dict(ntuple(_ -> 0, R) => constant))
end

"""
    _accumulate_term!(terms, exponent, coefficient)

Add a canonical exact coefficient to an owned sparse dictionary, deleting
cancellation zeros. Missing terms are assigned directly: allocating a fresh
`0//1` and adding it to every new coefficient is unnecessary GMP work.
"""
function _accumulate_term!(terms, exponent, coefficient::_BigRational)
  iszero(coefficient) && return terms
  previous = get(terms, exponent, nothing)
  value = previous === nothing ? coefficient : previous + coefficient
  iszero(value) ? delete!(terms, exponent) : (terms[exponent] = value)
  terms
end

"""
    _add_novikov_term!(polynomial, beta, coefficient)

Accumulate one validated quantum coefficient into an owned polynomial buffer.
The caller supplies compact nonnegative integer coordinates of the correct
length. This mutation is private to matrix assembly, before results escape.
"""
function _add_novikov_term!(
  p::NovikovPolynomial{R}, beta, coefficient::_BigRational
) where {R}
  _accumulate_term!(p.terms, ntuple(i -> beta[i], R), coefficient)
  p
end

Base.zero(::Type{NovikovPolynomial{R}}) where {R} = NovikovPolynomial{R}()
Base.zero(::NovikovPolynomial{R}) where {R} = zero(NovikovPolynomial{R})
Base.one(::Type{NovikovPolynomial{R}}) where {R} = NovikovPolynomial{R}(1)
Base.one(::NovikovPolynomial{R}) where {R} = one(NovikovPolynomial{R})
Base.iszero(p::NovikovPolynomial) = isempty(p.terms)
Base.broadcastable(p::NovikovPolynomial) = Ref(p)
Base.:(==)(a::NovikovPolynomial{R}, b::NovikovPolynomial{R}) where {R} = a.terms == b.terms
Base.:(==)(::NovikovPolynomial, ::NovikovPolynomial) = false
Base.isequal(a::NovikovPolynomial{R}, b::NovikovPolynomial{R}) where {R} =
  isequal(a.terms, b.terms)
Base.isequal(::NovikovPolynomial, ::NovikovPolynomial) = false
Base.hash(p::NovikovPolynomial{R}, h::UInt) where {R} = hash(p.terms, hash(R, h))
Base.copy(p::NovikovPolynomial{R}) where {R} =
  NovikovPolynomial{R}(copy(p.terms), Val(:canonical))

"""Merge two normalized polynomials with sign ±1, copying each dictionary at most once."""
function _combine_novikov(
  a::NovikovPolynomial{R}, b::NovikovPolynomial{R}, sign::Int
) where {R}
  terms = copy(a.terms)
  for (exponent, coefficient) in b.terms
    _accumulate_term!(terms, exponent, sign == 1 ? coefficient : -coefficient)
  end
  NovikovPolynomial{R}(terms, Val(:canonical))
end

Base.:+(a::NovikovPolynomial{R}, b::NovikovPolynomial{R}) where {R} =
  _combine_novikov(a, b, 1)
Base.:-(a::NovikovPolynomial{R}, b::NovikovPolynomial{R}) where {R} =
  _combine_novikov(a, b, -1)
Base.:-(p::NovikovPolynomial) = -1 * p

function Base.:*(coefficient::Union{Integer,Rational}, p::NovikovPolynomial{R}) where {R}
  isfinite(coefficient) || throw(ArgumentError("Novikov coefficients must be finite"))
  iszero(coefficient) && return zero(p)
  value = _BigRational(coefficient)
  terms = Dict(exponent => value * c for (exponent, c) in p.terms)
  NovikovPolynomial{R}(terms, Val(:canonical))
end
Base.:*(p::NovikovPolynomial, coefficient::Union{Integer,Rational}) = coefficient * p

function Base.:*(a::NovikovPolynomial{R}, b::NovikovPolynomial{R}) where {R}
  terms = Dict{NTuple{R,Int},_BigRational}()
  for (left, x) in a.terms, (right, y) in b.terms
    exponent = ntuple(i -> Base.Checked.checked_add(left[i], right[i]), R)
    _accumulate_term!(terms, exponent, x * y)
  end
  NovikovPolynomial{R}(terms, Val(:canonical))
end

for operator in (:+, :-, :*)
  @eval Base.$operator(::NovikovPolynomial, ::NovikovPolynomial) =
    throw(DimensionMismatch("Novikov polynomials must have the same number of variables"))
end

function Base.show(io::IO, p::NovikovPolynomial{R}) where {R}
  isempty(p.terms) && return print(io, "0")
  first_term = true
  for (exponent, coefficient) in sort!(collect(p.terms); by=first)
    negative = coefficient < 0
    print(io, first_term ? (negative ? "-" : "") : (negative ? " - " : " + "))
    first_term = false
    constant = all(iszero, exponent)
    magnitude = abs(coefficient)
    if constant || magnitude != 1
      print(io, magnitude)
      constant || print(io, '*')
    end
    first_variable = true
    for i in 1:R
      power = exponent[i]
      iszero(power) && continue
      first_variable || print(io, '*')
      first_variable = false
      print(io, 'q', i)
      power == 1 || print(io, '^', power)
    end
  end
end

"""Convert a finite rationally representable value to the exact evaluation scalar."""
function _novikov_value(value)
  value isa Union{Integer,Rational,AbstractFloat} && isfinite(value) ||
    throw(ArgumentError("Novikov values must be finite integers, rationals, or floats"))
  _BigRational(value)
end

"""
    evaluate(p::NovikovPolynomial)
    evaluate(p::NovikovPolynomial, q)

Evaluate exactly at the values of a tuple or vector `q`. The default specializes
every variable to one. Integers and rationals retain their exact values;
finite floating-point inputs are converted to their exact binary rational
values. Vector iteration order determines variable order, including for
vectors with nonstandard axes. The result is a `Rational{BigInt}`.
The values must be finite and their count must match the polynomial's number
of variables; malformed inputs raise `ArgumentError`.

# Examples
```jldoctest
julia> p = NovikovPolynomial{1}(Dict((2,) => 3));

julia> evaluate(p)
3//1

julia> evaluate(p, (1//2,))
3//4
```
"""
evaluate(p::NovikovPolynomial) = sum(values(p.terms); init=zero(_BigRational))

function evaluate(p::NovikovPolynomial{R}, q) where {R}
  q isa Union{Tuple,AbstractVector} && length(q) == R ||
    throw(ArgumentError("expected $R Novikov values in a tuple or vector"))
  values = map(_novikov_value, Tuple(q))
  total = zero(_BigRational)
  for (exponent, coefficient) in p.terms
    term = coefficient
    for i in 1:R
      iszero(exponent[i]) || (term *= values[i]^exponent[i])
    end
    total += term
  end
  total
end
