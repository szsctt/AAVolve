#!/usr/bin/env python3
import mappy as mp
import inspect

print('mappy version:', getattr(mp, '__version__', None))
try:
    a = mp.Aligner(seq='ACGT', preset='map-ont')
except Exception as e:
    print('Failed to create Aligner:', e)
    raise

hits = list(a.map('ACGT'))
print('n_hits:', len(hits))
print('Aligner.map signature:', inspect.signature(mp.Aligner.map))
if hits:
    h = hits[0]
    print('Sample hit repr:', repr(h))
    attrs = [x for x in dir(h) if not x.startswith('_')]
    print('hit attributes:', attrs)
    for attr in ('cigar', 'cigar_str', 'cigarstring', 'cigar_strs'):
        val = getattr(h, attr, None)
        print(f"{attr}:", val, 'type=', type(val))
    # also inspect types of elements if cigar is a list
    cig = getattr(h, 'cigar', None)
    if isinstance(cig, (list, tuple)) and len(cig):
        print('First cigar element repr:', repr(cig[0]))
        print('First cigar element type:', type(cig[0]))
        try:
            print('First cigar element items:', tuple(cig[0]))
        except Exception as e:
            print('Could not unpack first cigar element:', e)
else:
    print('No hits to inspect (tiny reference may have no hits)')
