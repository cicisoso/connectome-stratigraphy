"""Select the neurons whose skeletons are needed and extract the FlyWire ones into a compact cache."""
from common import *
fw = load_flywire()
prim_hl = set(fw.loc[(fw.birth == 'primary') & fw.hl.notna(), 'hl'])
need = (fw.cls == 'Kenyon_Cell') | fw.hl.isin(prim_hl | {'ALad1', 'ALl1_ventral', 'ALv1', 'CREa1_dorsal', 'CREa2_medial'})
ids = fw.loc[need, 'nid'].astype('int64').tolist(); print('FlyWire neurons for the cache', len(ids))
n = build_flywire_cache(ids); print('cached nodes', n)
