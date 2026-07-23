"""Data models for plates, groups, and wells."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Dict, Optional


PLATE_KEY_SEP = "::"


def make_well_key(plate_id, well):
    return plate_id + PLATE_KEY_SEP + well


def parse_plate_from_key(key):
    return key.split(PLATE_KEY_SEP, 1)[0]


def parse_well_from_key(key):
    return key.split(PLATE_KEY_SEP, 1)[1] if PLATE_KEY_SEP in key else key


def well_id_from_row_col(row, col):
    return chr(ord("A") + row - 1) + str(col)


@dataclass
class Plate:
    id: str
    name: str
    format: int = 96
    well_data: Dict[str, Optional[float]] = field(default_factory=dict)

    @property
    def layout(self):
        if self.format == 96:
            return {"rows": 8, "cols": 12, "max": 96}
        if self.format == 24:
            return {"rows": 4, "cols": 6, "max": 24}
        return {"rows": 2, "cols": 3, "max": 6}


@dataclass
class Group:
    id: str
    name: str
    color: str = "#0891b2"
    wells: List[str] = field(default_factory=list)


@dataclass
class Session:
    plates: List[Plate] = field(default_factory=list)
    active_plate_id: Optional[str] = None
    groups: List[Group] = field(default_factory=list)
    selected_wells: List[str] = field(default_factory=list)
    control_group_id: Optional[str] = None
    stats_mode: str = "aggregate"
    chart_mode: str = "abs"
