#!/usr/bin/env python3
"""
Small helper to normalise JUCE's generated moduleinfo.json files after a build.
The generator occasionally emits trailing commas; the VST3 validator rejects them.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from typing import Iterable


def _fix_moduleinfo(path: pathlib.Path) -> bool:
    if not path.exists():
        return False

    try:
        text = path.read_text()
    except OSError as exc:
        print(f"[mirror7] Unable to read {path}: {exc}", file=sys.stderr)
        return False

    text = re.sub(r",\s*(\]|\})", r"\1", text)

    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        print(f"[mirror7] Invalid JSON in {path}: {exc}", file=sys.stderr)
        return False

    try:
        path.write_text(json.dumps(payload, indent=2) + "\n")
    except OSError as exc:
        print(f"[mirror7] Unable to write {path}: {exc}", file=sys.stderr)
        return False

    return True


def _collect_targets(artefact_root: pathlib.Path, vst3_name: str, extra_dirs: Iterable[pathlib.Path]) -> list[pathlib.Path]:
    rel_module = pathlib.Path(f"{vst3_name}.vst3") / "Contents" / "Resources" / "moduleinfo.json"

    locations = [
        artefact_root / "VST3" / rel_module,
        artefact_root / "Debug" / "VST3" / rel_module,
        artefact_root / "Release" / "VST3" / rel_module,
    ]

    for extra in extra_dirs:
        locations.append(extra / rel_module)

    return locations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Patch trailing commas in JUCE's moduleinfo.json export.")
    parser.add_argument("--artefact-root", type=pathlib.Path, required=True, help="Root folder where JUCE copies plug-in artefacts (Mirror7_artefacts).")
    parser.add_argument("--vst3-name", type=str, default="Mirror7", help="VST3 bundle base name (defaults to Mirror7).")
    parser.add_argument("--user-vst3-dir", type=pathlib.Path, default=None, help="Optional manual user VST3 directory.")

    args = parser.parse_args(argv)

    extra_dirs = []
    if args.user_vst3_dir is not None:
        extra_dirs.append(args.user_vst3_dir.expanduser())
    else:
        extra_dirs.append(pathlib.Path.home() / "Library" / "Audio" / "Plug-Ins" / "VST3")

    artefact_root = args.artefact_root
    if not artefact_root.exists():
        print(f"[mirror7] Artefact root {artefact_root} does not exist yet; skipping module info patch.")
        return 0

    targets = _collect_targets(artefact_root, args.vst3_name, extra_dirs)

    patched = 0
    for candidate in targets:
        if _fix_moduleinfo(candidate):
            patched += 1

    print(f"[mirror7] Patched {patched} moduleinfo.json file(s)." if patched else "[mirror7] No moduleinfo.json files required patching.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
