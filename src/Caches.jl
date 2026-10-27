"""
    clear_localization_caches!()

Release cached flag root data, Schubert recurrence plans, tree shapes,
rooted depth sequences, and free-tree automorphisms.
Per-evaluation localization caches are local to each computation.
This does not clear caches owned by PartialFlagVarieties or Semisimple.
"""
function clear_localization_caches!()
  lock(_GW_DATA_LOCK) do
    empty!(_GW_DATA_CACHE)
  end
  lock(_TREE_AUTOMORPHISM_LOCK) do
    empty!(_TREE_AUTOMORPHISM_CACHE)
    empty!(_TREE_SHAPE_CACHE)
    empty!(_ROOTED_LEVEL_CACHE)
  end
  nothing
end
