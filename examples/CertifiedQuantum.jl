using PartialFlagVarieties, GromovWittenFlags

X = projective_space(2)
reference = small_quantum_multiplication(X; algorithm=:branches)
mktempdir() do checkpoints
  result = small_quantum_multiplication(X; arithmetic=:certified, checkpoint=checkpoints)
  @assert result.matrix == reference.matrix
  replay = small_quantum_multiplication(X; arithmetic=:certified, checkpoint=checkpoints)
  @assert replay.matrix == result.matrix
  @assert all(c -> c["new_prime_jobs"] == 0, replay.certification)
end
println("Certified quantum multiplication and checkpoint replay verified.")
