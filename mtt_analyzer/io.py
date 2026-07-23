"""I/O helpers for parsing CSV input and well positions."""

from __future__ import annotations

import csv
import re
from typing import List, Dict, Optional, Tuple


def parse_concentration(text):
    """Parse a concentration string like "1uM", "10 nM", "5e-6"."""
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        pass
    m = re.match(r"^\s*([0-9.eE+-]+)\s*(nM|uM|mM|M|pM|nm|um|mm|m|mg/mL|ug/mL|ng/mL)\s*$", s)
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2)
    table = {
        "M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9, "pM": 1e-12,
        "m": 1e-3, "mm": 1e-3, "um": 1e-6, "nm": 1e-9,
    }
    if unit in ("mg/mL", "ug/mL", "ng/mL"):
        return val
    return val * table.get(unit, 1.0)


def parse_csv(path):
    """Read a CSV file and return a list of dict rows."""
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: (v.strip() if v else v) for k, v in row.items()})
    return rows


def expand_wells(pattern, layout=None):
    """Expand a well pattern into a list of wells."""
    if layout is None:
        rows, cols = 8, 12
    else:
        rows, cols = layout
    out = []
    for token in re.split(r"[,\s\n]+", pattern):
        token = token.strip()
        if not token:
            continue
        m = re.match(r"^([A-P])(\d+)-([A-P])(\d+)$", token)
        if m:
            r1, c1, r2, c2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
            ra, rb = min(r1, r2), max(r1, r2)
            ca, cb = min(c1, c2), max(c1, c2)
            for r in range(ord(ra), ord(rb) + 1):
                for c in range(ca, cb + 1):
                    out.append(chr(r) + str(c))
            continue
        m = re.match(r"^([A-P])(\d+)$", token)
        if m:
            out.append(m.group(1) + m.group(2))
            continue
        m = re.match(r"^\[([A-P])-([A-P])\](\d+)$", token)
        if m:
            r1, r2, c = m.group(1), m.group(2), m.group(3)
            for r in range(ord(r1), ord(min(r2, chr(ord("A") + rows - 1))) + 1):
                out.append(chr(r) + c)
            continue
    seen = set()
    deduped = []
    for w in out:
        if w not in seen:
            seen.add(w)
            deduped.append(w)
    return deduped
