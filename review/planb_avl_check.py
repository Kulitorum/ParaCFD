"""Independent thin-surface AVL check of the unchanged production CAD.

AVL is an external reference executable, not a replacement load model in ParaCFD.
Section coordinates retain their incidence to horizontal flow. Their vertical
camber is projected on the local normal of the arched lifting surface; AVL uses
camber slopes in its linearized tangency condition. Thickness, viscosity and
separation are not represented by this comparison.
"""
from pathlib import Path
import json
import re
import subprocess
import numpy as np
from planb_section_audit import load_mesh, section_contour
from section_panel_check import resample, naca

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / 'build-paraglider-ui/reference-tools'


def run_avl(name, geometry, angle=0):
    (FOLDER / (name+'.avl')).write_text(geometry, encoding='utf-8')
    forces = FOLDER / (name+'.forces')
    strips = FOLDER / (name+'.strips')
    # Explicit overwrite answer only when the output already exists.
    commands = (f'PLOP\nG\n\nLOAD {name}.avl\nOPER\nA A {angle}\nX\n'
                f'FT\n{name}.forces\n'+('Y\n' if forces.exists() else '')+
                f'FS\n{name}.strips\n'+('Y\n' if strips.exists() else '')+
                '\nQUIT\n')
    (FOLDER / (name+'.commands')).write_text(commands, encoding='utf-8')
    run = subprocess.run([str(FOLDER / 'avl352.exe')], cwd=FOLDER,
                         input=commands, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, timeout=180,
                         creationflags=subprocess.CREATE_NO_WINDOW)
    (ROOT / 'review' / (name+'.log')).write_text(run.stdout, encoding='utf-8')
    if run.returncode or not forces.exists():
        raise RuntimeError(f'AVL failed for {name}: {run.stdout[-1800:]}')
    text = forces.read_text(encoding='utf-8')
    result = {key: float(re.search(r'\b'+key+r'\s*=\s*([-+.\dEe]+)', text)[1])
              for key in ('CLtot', 'CDtot', 'CDind', 'CLff', 'CDff', 'CYtot',
                          'Sref', 'Cref', 'Bref', 'Alpha')}
    result['name'] = name
    reasons = []
    if result['CDind'] < -1e-5:
        reasons.append('negative induced drag')
    if abs(result['CLtot']-result['CLff']) > .02:
        reasons.append('near-field and wake lift disagree')
    if abs(result['CYtot']) > .01:
        reasons.append('spurious sideforce on the symmetric input')
    result['accepted'] = not reasons
    result['rejection_reasons'] = reasons
    print(json.dumps(result), flush=True)
    return result


def rectangular_reference(code, span=20, angle=0):
    name = f'final-avl-naca{code}-ar{span}-alpha{angle}'
    coordinates = resample(naca(code, 400), 240)
    np.savetxt(FOLDER / (name+'.dat'), coordinates, header=name, comments='')
    geometry = (f'{name}\n0\n0 0 0\n{span} 1 {span}\n0 0 0\n'
                'SURFACE\nReference wing\n24 1.0 96 1.0\n'
                f'SECTION\n0 {-span/2} 0 1 0\nAFILE\n{name}.dat\n'
                f'SECTION\n0 {span/2} 0 1 0\nAFILE\n{name}.dat\n')
    return run_avl(name, geometry, angle)


def planb_geometry(chord_panels=20, span_panels=96, sections=49):
    vertices, faces = load_mesh(ROOT / 'review/final-planb-placed.obj')
    # Stay just inside the CAD end caps; retain and report the omitted span.
    fractions = .001 + .998 * .5 * (1-np.cos(np.linspace(0, np.pi, sections)))
    span = np.ptp(vertices[:, 1])
    ys = vertices[:, 1].min() + span * fractions
    contours = [section_contour(vertices, faces, y) for y in ys]
    leading = np.array([c[c[:, 0].argmin()] for c in contours])
    chords = np.array([np.ptp(c[:, 0]) for c in contours])
    cos_dihedral = 1 / np.sqrt(1+np.gradient(leading[:, 1], ys)**2)
    area = float(np.trapezoid(chords, ys))
    name = f'final-avl-planb-c{chord_panels}-s{span_panels}'
    geometry = (f'{name}; unchanged pose, horizontal flow\n0\n0 0 0\n'
                f'{area:.12g} 2.22899 {span:.12g}\n0 0 0\n'
                f'SURFACE\nPlanB\n{chord_panels} 1.0 {span_panels} 1.0\n')
    for index, (y, contour, le, chord, cosine) in enumerate(
            zip(ys, contours, leading, chords, cos_dihedral)):
        normalized = (contour-le) / chord
        normalized[:, 1] *= cosine
        dat = f'final-avl-planb-section-{index}.dat'
        np.savetxt(FOLDER / dat, resample(normalized, 240), header=dat, comments='')
        geometry += (f'SECTION\n{le[0]:.12g} {y:.12g} {le[1]:.12g} {chord:.12g} 0\n'
                     f'AFILE\n{dat}\n')
    return name, geometry, dict(projected_chord_integral_m2=area,
                               omitted_span_fraction=.002, sections=sections,
                               vertical_camber_projection='local surface normal',
                               thickness_correction=False)


if __name__ == '__main__':
    results = [rectangular_reference('0012', angle=0),
               rectangular_reference('0012', angle=4),
               rectangular_reference('2412', angle=0)]
    assert abs(results[0]['CLtot']) < 1e-6
    assert .35 < results[1]['CLtot'] < .48
    assert .15 < results[2]['CLtot'] < .27
    for nc, ns in ((12, 64), (20, 96), (28, 128)):
        name, geometry, metadata = planb_geometry(nc, ns)
        result = run_avl(name, geometry)
        result.update(metadata)
        result['lift_N_at_10ms_rho_1_225'] = result['CLtot'] * result['Sref'] * .5*1.225*100
        results.append(result)
    (ROOT / 'review/final-avl-results.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    if any(not result['accepted'] for result in results):
        raise SystemExit('REJECTED: the full CAD section-to-AVL approximation did not pass consistency checks.')
