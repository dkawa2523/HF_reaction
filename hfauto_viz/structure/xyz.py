from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path

@dataclass
class Atom:
    symbol: str
    x: float
    y: float
    z: float
Frame = list[Atom]

def parse_xyz_frames(text: str) -> list[Frame]:
    lines=[ln.rstrip() for ln in text.splitlines() if ln.strip()]
    frames=[]; i=0
    while i < len(lines):
        try: n=int(lines[i].split()[0])
        except Exception: break
        i += 1
        if i < len(lines): i += 1
        atoms=[]
        for _ in range(n):
            if i >= len(lines): break
            parts=lines[i].split(); i+=1
            if len(parts)>=4:
                try: atoms.append(Atom(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                except Exception: pass
        if atoms: frames.append(atoms)
    return frames

def read_xyz_frames(path: str | Path | None) -> list[Frame]:
    if not path: return []
    p=Path(path)
    if not p.exists(): return []
    return parse_xyz_frames(p.read_text(encoding='utf-8'))

def frame_to_xyz(frame: Frame, comment: str = "") -> str:
    return str(len(frame)) + "\n" + comment + "\n" + "\n".join(f"{a.symbol} {a.x:.8f} {a.y:.8f} {a.z:.8f}" for a in frame) + "\n"

def frames_to_xyz(frames: list[tuple[str, Frame]] | list[Frame]) -> str:
    out=[]
    for i,item in enumerate(frames):
        if isinstance(item, tuple): comment, frame = item
        else: comment, frame = f"frame {i}", item
        out.append(frame_to_xyz(frame, comment))
    return "".join(out)

def read_xyz(path: str | Path | None) -> str:
    frames = read_xyz_frames(path)
    return frames_to_xyz(frames) if frames else ""

def reaction_path_frames(reactant: Frame | None, ts: Frame | None, product: Frame | None) -> list[tuple[str, Frame]]:
    frames=[]
    if reactant: frames.append(("reactant", reactant))
    if ts: frames.append(("transition_state", ts))
    if product: frames.append(("product", product))
    return frames


def interpolate_frames(a: Frame, b: Frame, n: int = 10, prefix: str = 'frame') -> list[tuple[str, Frame]]:
    if not a or not b or len(a) != len(b):
        return [(prefix + '_0', a)] if a else []
    frames=[]
    steps=max(n,2)
    for i in range(steps):
        t=i/(steps-1)
        frame=[]
        for aa,bb in zip(a,b):
            frame.append(Atom(aa.symbol, aa.x*(1-t)+bb.x*t, aa.y*(1-t)+bb.y*t, aa.z*(1-t)+bb.z*t))
        frames.append((f'{prefix}_{i:03d}', frame))
    return frames

def concatenate_xyz_frames(frames: list[tuple[str, Frame]] | list[Frame]) -> str:
    return frames_to_xyz(frames)

def read_xyz_tuple(path: str | Path | None) -> tuple[str, Frame]:
    frames=read_xyz_frames(path)
    return ('frame 0', frames[0]) if frames else ('', [])
