#!/bin/zsh
# Regenerates every result, table and figure of the Biology Letters package from the public data.
# Data locations: see code/common.py (FLYBRAIN_DATA). Run from biolett/code.
set -e
python a0_prepare.py
for ds in MaleCNS FlyWire BANC; do python -W ignore a1_kc.py $ds; done
for ds in MaleCNS FlyWire BANC; do python -W ignore a2_lineages.py $ds; done
for ds in MaleCNS FlyWire BANC; do python -W ignore a3_hemilineages.py $ds; done
for ds in MaleCNS FlyWire BANC MANC; do python -W ignore a4_connectivity.py $ds; done
for ds in MaleCNS FlyWire BANC; do python -W ignore a6_render_data.py $ds; done
python -W ignore a5_summary.py
for f in fig1.py fig2.py fig3.py figS_renderings.py figS_more.py; do python -W ignore $f; done
python -W ignore make_tables.py
