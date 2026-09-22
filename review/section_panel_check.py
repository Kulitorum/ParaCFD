"""Independent 2-D source/vortex panel check; not a production lift substitute.

Constant source strength per straight panel and one constant vortex-sheet
strength; no penetration at panel midpoints and equal TE surface speeds.
Reference formulation: https://web.itu.edu.tr/~atares/courses/CA/4.2_HessSmith.html
The code shares no pressure/transport implementation with ParaCFD or XFOIL.
"""
from pathlib import Path
import json
import numpy as np


def naca(code, n=800):
    m, p, thickness = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = .5 * (1 - np.cos(np.linspace(0, np.pi, n // 2 + 1)))
    yt = 5 * thickness * (.2969 * np.sqrt(x) - .126 * x - .3516 * x*x
                          + .2843 * x**3 - .1036 * x**4)
    yc, dy = np.zeros_like(x), np.zeros_like(x)
    if m:
        front = x < p
        yc[front] = m / p**2 * (2*p*x[front] - x[front]**2)
        dy[front] = 2*m / p**2 * (p-x[front])
        yc[~front] = m / (1-p)**2 * (1-2*p+2*p*x[~front]-x[~front]**2)
        dy[~front] = 2*m / (1-p)**2 * (p-x[~front])
    theta = np.arctan(dy)
    upper = np.column_stack((x-yt*np.sin(theta), yc+yt*np.cos(theta)))
    lower = np.column_stack((x+yt*np.sin(theta), yc-yt*np.cos(theta)))
    return np.vstack((upper[::-1], lower[1:]))


def resample(points, count):
    """Cosine spacing separately on each branch, on the original polyline."""
    p = points.copy()
    if np.linalg.norm(p[0]-p[-1]) > 1e-12:
        p = np.vstack((p, p[0]))
    le = np.argmin(p[:, 0])
    result = []
    for branch in (p[:le+1], p[le:]):
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(branch, axis=0), axis=1))]
        targets = arc[-1] * .5 * (1-np.cos(np.linspace(0, np.pi, count//2+1)))
        result.append(np.column_stack([np.interp(targets, arc, branch[:, a]) for a in range(2)]))
    return np.vstack((result[0], result[1][1:]))


def solve(points, flow_degrees=0):
    p = np.asarray(points, dtype=float)
    p = p[np.r_[True, np.linalg.norm(np.diff(p, axis=0), axis=1) > 1e-12]]
    if np.linalg.norm(p[0]-p[-1]) > 1e-12:
        p = np.vstack((p, p[0]))
    if np.sum(p[:-1, 0]*p[1:, 1]-p[1:, 0]*p[:-1, 1]) < 0:
        p = p[::-1]
    start, end = p[:-1], p[1:]
    centre = (start+end)/2
    length = np.linalg.norm(end-start, axis=1)
    tangent = (end-start)/length[:, None]
    normal = np.column_stack((tangent[:, 1], -tangent[:, 0]))
    left = -normal
    r = centre[:, None, :]-start[None, :, :]
    x = np.einsum('ijk,jk->ij', r, tangent)
    y = np.einsum('ijk,jk->ij', r, left)
    du = np.log((x*x+y*y)/((x-length)**2+y*y))/(4*np.pi)
    angle = np.arctan2(y, x-length)-np.arctan2(y, x)
    angle = (angle+np.pi) % (2*np.pi)-np.pi
    dv = angle/(2*np.pi)
    # Exterior limit for a CCW contour: source velocity points to the right.
    np.fill_diagonal(du, 0)
    np.fill_diagonal(dv, -.5)
    source = du[:, :, None]*tangent[None, :, :]+dv[:, :, None]*left[None, :, :]
    vortex = np.stack((-source[:, :, 1], source[:, :, 0]), axis=2)
    influence_n = np.einsum('ijk,ik->ij', source, normal)
    influence_t = np.einsum('ijk,ik->ij', source, tangent)
    vortex_n = np.einsum('ijk,ik->ij', vortex, normal).sum(axis=1)
    vortex_t = np.einsum('ijk,ik->ij', vortex, tangent).sum(axis=1)
    a = np.deg2rad(flow_degrees)
    freestream = np.array([np.cos(a), np.sin(a)])
    system = np.zeros((len(length)+1, len(length)+1))
    system[:-1, :-1] = influence_n
    system[:-1, -1] = vortex_n
    system[-1, :-1] = influence_t[0]+influence_t[-1]
    system[-1, -1] = vortex_t[0]+vortex_t[-1]
    rhs = np.r_[-normal@freestream, -(tangent[0]+tangent[-1])@freestream]
    strength = np.linalg.solve(system, rhs)
    velocity = tangent@freestream+influence_t@strength[:-1]+vortex_t*strength[-1]
    cp = 1-velocity**2
    chord = np.ptp(p[:, 0])
    force = -np.sum(cp[:, None]*normal*length[:, None], axis=0)/chord
    lift_direction = np.array([-freestream[1], freestream[0]])
    result = dict(panels=len(length), flow_degrees=flow_degrees,
                  cl=float(force@lift_direction), cd_pressure=float(force@freestream),
                  cl_circulation=float(-2*strength[-1]*sum(length)/chord),
                  source_flux=float(strength[:-1]@length),
                  residual=float(np.max(np.abs(system@strength-rhs))),
                  cp_min=float(cp.min()), cp_max=float(cp.max()))
    return result, np.column_stack((centre, cp, velocity))


if __name__ == '__main__':
    results = []
    for code in ('0012', '2412'):
        for angle in ((-4, 0, 4) if code == '0012' else (0, 4)):
            result, _ = solve(naca(code), angle)
            result['case'] = 'NACA'+code
            print(json.dumps(result), flush=True)
            results.append(result)
            if code == '0012' and angle == 0:
                assert abs(result['cl']) < 1e-6
            if code == '2412' and angle == 4:
                assert abs(result['cl']-.7414) < .005
            assert abs(result['cd_pressure']) < .003
            assert abs(result['cl']-result['cl_circulation']) < .01
    planb = np.loadtxt('review/planb-neutral-section-outline.csv', delimiter=',', skiprows=1)
    planb = (planb-planb[planb[:, 0].argmin()])/np.ptp(planb[:, 0])
    for count in (200, 400, 800, 1600):
        result, samples = solve(resample(planb, count))
        result['case'] = 'PlanB unchanged midspan, horizontal flow'
        print(json.dumps(result), flush=True)
        results.append(result)
        np.savetxt(f'review/final-panel-planb-{count}.csv', samples, delimiter=',',
                   header='x,z,cp,tangential_velocity', comments='')
    assert abs(results[-1]['cl']-results[-2]['cl']) < .001
    assert abs(results[-1]['cl']-results[-1]['cl_circulation']) < .001
    Path('review/final-panel-results.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
