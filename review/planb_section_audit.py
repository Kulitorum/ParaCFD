"""Audit unrotated span stations from the exported production triangle mesh."""
from pathlib import Path
from collections import defaultdict, Counter
import json
import numpy as np
from section_panel_check import solve, resample


def load_mesh(path):
    vertices, faces = [], []
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if line.startswith('v '):
            vertices.append([float(x) for x in line.split()[1:]])
        elif line.startswith('f '):
            faces.append([int(x)-1 for x in line.split()[1:]])
    return np.array(vertices), np.array(faces)


def section_contour(vertices, faces, y):
    tri = vertices[faces]
    points = defaultdict(list)
    for edge in range(3):
        a, b = tri[:, edge], tri[:, (edge+1) % 3]
        selected = ((a[:, 1] < y) & (b[:, 1] >= y)) | ((b[:, 1] < y) & (a[:, 1] >= y))
        indices = np.flatnonzero(selected)
        d = b[selected]-a[selected]
        crossings = a[selected]+d*((y-a[selected, 1])/d[:, 1])[:, None]
        for index, point in zip(indices, crossings):
            points[index].append(point[[0, 2]])
    adjacency = defaultdict(list)
    for pair in points.values():
        assert len(pair) == 2
        a, b = [tuple(np.round(p, 8)) for p in pair]
        if a != b:
            adjacency[a].append(b)
            adjacency[b].append(a)
    assert set(map(len, adjacency.values())) == {2}, Counter(map(len, adjacency.values()))
    start = max(adjacency)
    p, previous, contour = start, None, []
    while True:
        contour.append(p)
        following = [q for q in adjacency[p] if q != previous][0]
        previous, p = p, following
        if p == start:
            break
        assert len(contour) <= len(adjacency)
    assert len(contour) == len(adjacency), 'multiple section loops'
    return np.array(contour)


if __name__ == '__main__':
    vertices, faces = load_mesh('review/final-planb-placed.obj')
    span = np.ptp(vertices[:, 1])
    results = []
    for fraction in (.1, .25, .5001, .75, .9):
        y = vertices[:, 1].min()+span*fraction
        contour = section_contour(vertices, faces, y)
        chord = np.ptp(contour[:, 0])
        leading, trailing = contour[contour[:, 0].argmin()], contour[contour[:, 0].argmax()]
        normalized = (contour-leading)/chord
        result, samples = solve(resample(normalized, 800))
        result.update(span_fraction=fraction, source_y=y, chord_m=chord,
                      nose_up_chord_degrees=float(np.rad2deg(np.arctan2(leading[1]-trailing[1], chord))))
        results.append(result)
        print(json.dumps(result), flush=True)
        np.savetxt(f'review/final-planb-section-{fraction}.csv', contour,
                   header='x,z', delimiter=',', comments='')
    Path('review/final-span-section-results.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
