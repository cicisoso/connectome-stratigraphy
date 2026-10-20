Per-neuron measurements behind every figure (one CSV per connectome and analysis).
Coordinates in micrometres in each connectome's own space; nid = the connectome's neuron identifier.
kenyon_cells: kctype, rank (0 gamma-d, 1 gamma-m, 2 alpha'/beta', 3-6 alpha/beta p, s, m, c), clone, ped_radial and ped_burial (fraction of
  peduncle radius), lobe_radial (alpha/beta lobe), soma_tang_clone (um from the clone axis).
lineages_*: lineage, birth rank (see ESM table S11), epoch (adPN), frame, in_section, soma_dist, soma_depth, soma_tang (um), radial and burial
  (fraction of tract radius), geo_age (mean of the z-scored, sign-aligned soma depth, soma distance and fibre radial position).
hemilineage_units: unit (hemilineage|side|neuromere), birth class (primary, early_secondary, secondary), the same measures, geo_age_soma
  (soma measures only).
connectivity: synapse counts, size (MaleCNS and BANC: cable length in um; FlyWire: skeleton length; MANC: volume), synapses per size and
  the fraction of input and output synapses from each partner class.
