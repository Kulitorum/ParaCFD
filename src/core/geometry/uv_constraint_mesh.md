# Deterministic face-UV constraint insertion

`uv_constraint_mesh.{h,cpp}` is an OCC-free geometry component for converting one existing
face-local triangle mesh into a planar straight-line graph whose contact edges use a
single, shared sample chain. Its STEP-import call site is deliberately deferred until
contact canonicalization supplies the final common chains.

## Required preprocessing

For every physical CAD contact curve, form one ordered sample chain before touching any
face mesh. The chain must be the union of:

- the curve's common geometric samples;
- every existing trim-boundary vertex from every participating face that lies on the
  curve; and
- every intentional constraint intersection.

Every sample receives one canonical `uint32_t` topology ID. The same ID is passed with
the corresponding face-local UV point for every face use. An existing fine-side midpoint
that is omitted from the common chain is rejected rather than silently producing a
nonconforming one-to-two edge match.

## Insertion algorithm

Insertion is transactional per polyline.

1. Locate every prescribed sample in deterministic vertex/edge/triangle order. Reuse an
   existing vertex, split an existing edge, or split one containing triangle. Splitting a
   previously constrained edge propagates all of its constraint IDs to both children,
   which permits a later curve to end at an interior T-junction.
2. For each adjacent sample pair, audit for an unsampled collinear vertex. Such a vertex
   must be added to the common chain first.
3. Trace the open segment through the existing triangles by sorting its proper crossings
   with mesh edges. Crossing a trim-boundary or an earlier constraint is rejected unless
   their intersection was supplied as a sample.
4. Remove the crossed triangle strip, recover its directed CCW boundary, split that loop
   at the two prescribed endpoints, and deterministically ear-clip the two resulting
   polygons. Their common closing edge is exactly the requested segment.
5. Tag that edge with the physical constraint ID. No crossing-dependent Steiner vertex is
   introduced on the shared contact chain.

Crossing locations are ordered by their dimensionless segment parameter. Comparisons use
`coordinate_tolerance / segment_length`, so the classification is invariant under a uniform
change in the face's native UV scale; a UV length is never compared directly with a unitless
parameter.

One-sided mesh edges define both outer and inner trim loops. The algorithm never crosses
those edges, so a successful interior insertion cannot fill a hole or change its boundary.
Samples may split a boundary edge, preserving the same loop as two child edges.

Before insertion, validation proves the indexed two-manifold/winding invariants and audits the
geometry itself. Distinct indexed edges may meet only at a shared indexed endpoint; proper
crossings, collinear overlap, and unindexed T-junctions are rejected. A separating-axis audit of
every triangle pair also rejects positive-area overlap, including one triangle wholly contained
inside another without an edge crossing. These checks are deterministic and insertion remains
transactional when any check fails.

Determinism comes from ordered edge maps, stable triangle scans, sorted crossing
parameters, and choosing the valid ear with the smallest vertex index. Replaying identical
indexed input produces bitwise-identical vertices, triangles, and constraint tags.

## Production handoff

For each source face, construct this mesh from `Poly_Triangulation` UV nodes and CCW
triangles, insert every `StepContactFaceUse::sample_uv` chain, then emit face-local render
vertices. Shared samples retain their supplied topology IDs; a caller-owned global
allocator assigns unique IDs to all remaining face-local vertices. Boundary/constraint
edge tags can then populate CAD/contact half-edge provenance in the final `TriMesh`.

The prototype uses scan-based adjacency and point location to keep its invariants visible.
A production version should replace those scans with half-edge adjacency and a spatial
locator without changing ordering or transaction semantics.

## Current limits

- Predicates use `long double` plus a caller-selected absolute UV tolerance, not adaptive
  exact arithmetic. The importer must derive that tolerance from its face/pcurve audit.
- Vertex indices and canonical topology IDs are `uint32_t`. `UINT32_MAX` is permanently
  reserved as the unassigned topology-ID sentinel; validation and every vertex allocation fail
  deterministically before an index cast or sentinel collision can occur.
- Intersecting constraints require their intersection in the canonical sample chains.
- Periodic UV seams need to be unwrapped into one consistent face-local chart first.
- Ear clipping preserves topology but does not optimize triangle angles. A deterministic
  constrained-Delaunay cleanup may flip only unconstrained, non-boundary edges afterward.
- Complexity is intentionally quadratic in mesh size for the research probe.

## Production test

From the repository root:

```powershell
cmake --build build-paraglider-ui --config Release --target uv_constraint_mesh_probe
ctest --test-dir build-paraglider-ui -C Release -R '^uv_constraint_mesh$' --output-on-failure
```

The probe covers an interior T-junction, one-to-two boundary sampling normalization,
canonical/local topology ID assignment, annular holes and boundary splits, transactional
rejection, scale-aware crossing-parameter classification, geometric edge-crossing and contained
triangle-overlap rejection, malformed winding/non-manifold rejection, allocator exhaustion,
area preservation, and bitwise deterministic replay.
