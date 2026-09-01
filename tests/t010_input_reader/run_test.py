#!/usr/bin/env python3
"""Input-reader robustness tests (see src/input_reader.h).

1. Section-order independence: the same physical input must produce a
   bit-identical kinematics fingerprint for every permutation of the
   [section] order. Includes the orders that used to fail or silently
   misbehave: [lsgrid] before [wave input irregular], [wave type] last,
   and [wave input regular] before [general input data] (which used to
   pick up depth = 0 silently).
2. Keyword-order freedom inside [wave input irregular] (options after
   the component table).
3. A misspelled section name and an unknown keyword raise warnings.
4. A malformed number in the component table is a hard error (it used
   to become a silent 0.0 component).

Run from this directory (probe must be built first, see probe.cpp):
    python3 run_test.py
"""
import itertools
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROBE = os.path.join(HERE, 'probe')

GENERAL = """[general input data]
depth 50.0
mtheta 15.0
gravity 9.81
"""
REFPOINT = """[wave reference point]
time 0.0
x 5.0
y -3.0
"""
RAMPS = """[ramps]
time_rampup 0 8.0
"""
IRREGULAR = """[wave input irregular]
nfreq 3
ndir 0
0.60000000 1.00000000e+00 0.0367923842 0.50000000 0.10000000
0.90000000 5.00000000e-01 0.0827447770 2.10000000 -0.20000000
1.30000000 2.50000000e-01 0.1723155997 4.00000000 0.30000000
"""
LSGRID = """[lsgrid]
bounds 0.0 200.0 -40.0 40.0
nx 48
ny 12
nl 8
t0 0.0
dt 0.5
stretch_beta 3.0
"""
WAVETYPE = """[wave type]
irregular
"""

REGULAR_GENERAL = """[general input data]
depth 30.0
mtheta 0.0
gravity 9.81
"""
REGULAR_INPUT = """[wave input regular]
wave_height 4.0
wave_period 9.0
"""
REGULAR_TYPE = """[wave type]
regular
"""


def run_case(sections, expect_fail=False):
    """Write the sections in the given order and run the probe."""
    with tempfile.TemporaryDirectory() as td:
        with open(os.path.join(td, 'waveinput.dat'), 'w') as fh:
            fh.write('@v301\n' + '\n'.join(sections) + '\n')
        r = subprocess.run([PROBE], cwd=td, capture_output=True, text=True,
                           timeout=300)
        probe = [ln for ln in r.stdout.splitlines()
                 if ln.startswith('PROBE')]
        if expect_fail:
            return r.returncode, r.stdout + r.stderr
        if r.returncode != 0 or not probe:
            print(r.stdout[-3000:])
            print(r.stderr[-2000:])
            raise RuntimeError(f'probe failed (exit {r.returncode})')
        return probe, r.stdout + r.stderr


def main():
    nfail = 0

    # ---- 1. irregular + lsgrid: all 120 section orders ----
    base = [WAVETYPE, GENERAL, REFPOINT, RAMPS, IRREGULAR, LSGRID]
    ref, _ = run_case(base)
    orders = list(itertools.permutations(base))
    bad = 0
    for perm in orders:
        got, _ = run_case(list(perm))
        if got != ref:
            bad += 1
    print(f'[1] section-order independence (irregular+lsgrid): '
          f'{len(orders)} permutations, {bad} mismatches'
          + (' - OK' if bad == 0 else ' - FAIL'))
    nfail += (bad != 0)

    # ---- 2. regular wave: [wave input regular] before [general] ----
    a, _ = run_case([REGULAR_TYPE, REGULAR_GENERAL, REGULAR_INPUT])
    b, _ = run_case([REGULAR_INPUT, REGULAR_TYPE, REGULAR_GENERAL])
    ok = (a == b)
    print(f'[2] regular wave, [wave input regular] first: '
          + ('identical - OK' if ok else 'DIFFERS - FAIL (silent depth=0?)'))
    nfail += (not ok)

    # ---- 3. options after the component table ----
    irregular_reordered = """[wave input irregular]
0.60000000 1.00000000e+00 0.0367923842 0.50000000 0.10000000
0.90000000 5.00000000e-01 0.0827447770 2.10000000 -0.20000000
1.30000000 2.50000000e-01 0.1723155997 4.00000000 0.30000000
nfreq 3
ndir 0
"""
    got, _ = run_case([WAVETYPE, GENERAL, REFPOINT, RAMPS,
                       irregular_reordered, LSGRID])
    ok = (got == ref)
    print('[3] nfreq/ndir after the component table: '
          + ('identical - OK' if ok else 'DIFFERS - FAIL'))
    nfail += (not ok)

    # ---- 4. misspelled section + unknown keyword warn ----
    gen_typo = GENERAL + 'mtetha 5.0\n'
    _, err = run_case([WAVETYPE, gen_typo, REFPOINT, RAMPS, IRREGULAR,
                       LSGRID, '[seabeed]\nfile nosuch.stl\n'])
    ok = '[seabeed]' in err and 'mtetha' in err
    print('[4] misspelled section/keyword warns: '
          + ('OK' if ok else 'FAIL'))
    nfail += (not ok)

    # ---- 5. malformed number is a hard error ----
    broken = IRREGULAR.replace('0.0827447770', '0.08x7447770')
    rc, out = run_case([WAVETYPE, GENERAL, IRREGULAR.replace(
        IRREGULAR, broken), LSGRID], expect_fail=True)
    ok = rc != 0 and 'malformed' in out
    print('[5] malformed table number is a hard error: '
          + ('OK' if ok else f'FAIL (exit {rc})'))
    nfail += (not ok)

    # ---- 6. wave_length AND wave_period together is an error ----
    both = REGULAR_INPUT.rstrip() + '\nwave_length 120.0\n'
    rc, out = run_case([REGULAR_TYPE, REGULAR_GENERAL, both],
                       expect_fail=True)
    ok = rc != 0 and 'either' in out
    print('[6] wave_length + wave_period together is a hard error: '
          + ('OK' if ok else f'FAIL (exit {rc})'))
    nfail += (not ok)

    # ---- 7. stretch_beta auto resolves and is order-independent ----
    lsgrid_auto = LSGRID.replace('stretch_beta 3.0', 'stretch_beta auto')
    base_auto = [WAVETYPE, GENERAL, REFPOINT, RAMPS, IRREGULAR, lsgrid_auto]
    ref_auto, err = run_case(base_auto)
    resolved = 'auto -> beta=' in err
    # omitting the keyword entirely must give the same result (auto is
    # the default), and the section order must not matter
    lsgrid_missing = LSGRID.replace('stretch_beta 3.0\n', '')
    got_missing, _ = run_case([WAVETYPE, GENERAL, REFPOINT, RAMPS,
                               IRREGULAR, lsgrid_missing])
    got_perm, _ = run_case([lsgrid_auto, IRREGULAR, RAMPS, REFPOINT,
                            GENERAL, WAVETYPE])
    differs_from_manual = (ref_auto != ref)
    ok = (resolved and got_missing == ref_auto and got_perm == ref_auto
          and differs_from_manual)
    print('[7] stretch_beta auto: '
          + ('resolved, default, order-independent, != beta=3 - OK'
             if ok else f'FAIL (resolved={resolved}, '
             f'default={got_missing == ref_auto}, '
             f'perm={got_perm == ref_auto}, '
             f'distinct={differs_from_manual})'))
    nfail += (not ok)

    # ---- 8. malformed stretch_beta is a hard error ----
    rc, out = run_case([WAVETYPE, GENERAL, REFPOINT, RAMPS, IRREGULAR,
                        LSGRID.replace('stretch_beta 3.0',
                                       'stretch_beta fast')],
                       expect_fail=True)
    ok = rc != 0 and 'stretch_beta' in out
    print('[8] malformed stretch_beta value is a hard error: '
          + ('OK' if ok else f'FAIL (exit {rc})'))
    nfail += (not ok)

    print('all reader tests passed' if nfail == 0
          else f'{nfail} READER TEST(S) FAILED')
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
