#!/usr/bin/env python3
"""validate_vti.py - independent validator for the hand-rolled VTI writer.

Parses a VTK ImageData (.vti) file written with appended raw binary (header_type
UInt64, LittleEndian, Float32 point scalars) WITHOUT any VTK dependency, checks
structural consistency (WholeExtent vs each array's declared byte count), and,
with --check-synthetic, re-derives the M0 analytic field and compares.

The synthetic formula here MUST stay identical to src/core/vti_synthetic.h.

Exit code 0 => valid; non-zero => invalid (message on stderr).
"""
import argparse
import math
import re
import struct
import sys


def _attr(text, name):
    m = re.search(name + r'="([^"]*)"', text)
    return m.group(1) if m else None


def synthetic_value(i, j, k, nx, ny, nz, sx, sy, sz):
    x, y, z = i * sx, j * sy, k * sz
    Lx, Ly, Lz = nx * sx, ny * sy, nz * sz
    fx = math.sin(2.0 * math.pi * x / Lx) if Lx > 0 else 0.0
    fy = math.cos(2.0 * math.pi * y / Ly) if Ly > 0 else 0.0
    fz = (z / Lz) if Lz > 0 else 0.0
    return fx * fy + 0.5 * fz


def parse_vti(path):
    with open(path, "rb") as f:
        raw = f.read()

    # Header text ends at the appended underscore; decode it as latin-1 so byte
    # offsets line up 1:1 with the raw buffer.
    head = raw.decode("latin-1", errors="replace")

    if 'type="ImageData"' not in head:
        raise ValueError("not a VTK ImageData file")
    if _attr(head, "header_type") != "UInt64":
        raise ValueError("expected header_type=UInt64")
    if _attr(head, "byte_order") != "LittleEndian":
        raise ValueError("expected byte_order=LittleEndian")

    ext = _attr(head, "WholeExtent")
    if ext is None:
        raise ValueError("missing WholeExtent")
    e = [int(v) for v in ext.split()]
    if len(e) != 6:
        raise ValueError("bad WholeExtent")
    nx, ny, nz = e[1] - e[0] + 1, e[3] - e[2] + 1, e[5] - e[4] + 1
    if min(nx, ny, nz) < 1:
        raise ValueError("non-positive dimensions")

    org = [float(v) for v in _attr(head, "Origin").split()]
    spc = [float(v) for v in _attr(head, "Spacing").split()]
    if len(org) != 3 or len(spc) != 3:
        raise ValueError("bad Origin/Spacing")

    npoints = nx * ny * nz

    ap = head.find("<AppendedData")
    if ap < 0:
        raise ValueError("missing AppendedData")
    us = raw.find(b"_", ap)
    if us < 0:
        raise ValueError("missing appended-data underscore")
    base = us + 1

    # Collect DataArrays in PointData.
    pd = head.find("<PointData")
    pd_end = head.find("</PointData>")
    if pd < 0 or pd_end < 0:
        raise ValueError("missing PointData")
    fields = {}
    for m in re.finditer(r"<DataArray\b[^>]*?/?>", head[pd:pd_end]):
        tag = m.group(0)
        name = _attr(tag, "Name")
        typ = _attr(tag, "type")
        off = _attr(tag, "offset")
        if typ != "Float32":
            raise ValueError("DataArray %r not Float32" % name)
        offset = int(off)
        hpos = base + offset
        (nbytes,) = struct.unpack_from("<Q", raw, hpos)
        expect = npoints * 4
        if nbytes != expect:
            raise ValueError("array %r byte-count %d != expected %d" % (name, nbytes, expect))
        dpos = hpos + 8
        if dpos + nbytes > len(raw):
            raise ValueError("array %r data past EOF" % name)
        vals = struct.unpack_from("<%df" % npoints, raw, dpos)
        fields[name] = vals

    if not fields:
        raise ValueError("no DataArray found")

    return {
        "dims": (nx, ny, nz),
        "origin": org,
        "spacing": spc,
        "npoints": npoints,
        "fields": fields,
    }


def main():
    ap = argparse.ArgumentParser(description="Validate a hand-rolled .vti file")
    ap.add_argument("path")
    ap.add_argument("--check-synthetic", action="store_true",
                    help="also verify payload equals the M0 analytic field")
    ap.add_argument("--tol", type=float, default=1e-5,
                    help="absolute tolerance for --check-synthetic")
    args = ap.parse_args()

    try:
        info = parse_vti(args.path)
    except Exception as ex:  # noqa: BLE001
        print("INVALID: %s" % ex, file=sys.stderr)
        return 1

    nx, ny, nz = info["dims"]
    print("OK header: dims=%dx%dx%d origin=%s spacing=%s fields=%s"
          % (nx, ny, nz, info["origin"], info["spacing"], list(info["fields"])))

    # Sanity: all values finite.
    for name, vals in info["fields"].items():
        for v in vals:
            if not math.isfinite(v):
                print("INVALID: field %r has non-finite value" % name, file=sys.stderr)
                return 1

    if args.check_synthetic:
        if "synthetic" not in info["fields"]:
            print("INVALID: no 'synthetic' field to check", file=sys.stderr)
            return 1
        sx, sy, sz = info["spacing"]
        vals = info["fields"]["synthetic"]
        max_err = 0.0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    idx = i + nx * (j + ny * k)
                    ref = struct.unpack("<f", struct.pack("<f", synthetic_value(i, j, k, nx, ny, nz, sx, sy, sz)))[0]
                    err = abs(vals[idx] - ref)
                    if err > max_err:
                        max_err = err
        print("synthetic check: max_abs_err=%.3e tol=%.1e" % (max_err, args.tol))
        if max_err > args.tol:
            print("INVALID: synthetic field mismatch", file=sys.stderr)
            return 1

    print("VALID: %s" % args.path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
