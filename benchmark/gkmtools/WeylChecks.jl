"""Check word transport and simple-root actions across the two independent Weyl implementations."""
function weyl_checks()
  cases = ((:A, 3, SS.TypeA{3}), (:B, 2, SS.TypeB{2}),
    (:C, 3, SS.TypeC{3}), (:G, 2, SS.TypeG2))
  rows = Dict{String,Any}[]
  for (kind, rank, DT) in cases
    root_system = OS.root_system(kind, rank)
    W = OS.weyl_group(root_system)
    reflections = OS.reflection.(OS.simple_roots(root_system))
    semisimple_roots = SS.RootSystem(DT)
    elements = SS.right_coset_reps(SS.weyl_group(DT), Int[])
    images = typeof(one(W))[]
    action_mismatches = 0
    for w in elements
      image = word_to_weyl(SS.word(w), W, reflections)
      length(w) == length(image) || error("word transport changed length")
      push!(images, image)
      for i in 1:rank
        left = collect(SS.coefficients(SS.simple_root(semisimple_roots, i) * w))
        right = [Int((OS.simple_root(root_system, i) * image)[j]) for j in 1:rank]
        left == right || (action_mismatches += 1)
      end
    end
    length(unique(images)) == length(elements) == Int(OS.order(W)) ||
      error("the transported words do not give a bijection")
    A, B = SS.cartan_matrix(DT), OS.cartan_matrix(root_system)
    cartan_equal = all(A[i, j] == Int(B[i, j]) for i in 1:rank, j in 1:rank)
    cartan_equal || error("the Cartan conventions disagree for $kind$rank")
    iszero(action_mismatches) || error("the root actions disagree for $kind$rank")
    push!(
      rows,
      Dict("type" => "$kind$rank", "order" => length(elements),
        "cartan_equal" => cartan_equal,
        "simple_root_action_checks" => rank * length(elements),
        "simple_root_action_mismatches" => action_mismatches,
        "Semisimple_type" => string(typeof(first(elements))),
        "Oscar_type" => string(typeof(first(images)))),
    )
  end
  Dict("julia" => string(VERSION), "cases" => rows)
end
