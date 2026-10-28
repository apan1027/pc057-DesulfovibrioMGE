# Generated analysis data

Run `python3 dev/prepare_inputs.py` from the repository root. The compressed archive is verified before extraction. All 72 extracted payloads are checked against `metadata/input-manifest.json`; existing input files are never silently overwritten. An interrupted or changed input preparation fails verification on the next invocation.

Shared inputs are placed under `00-raw/d00-resource/pc057/`. Generated outputs use the numbered step directories. Git ignores extracted inputs and bulk outputs; it tracks only this explanation and seven selected PNG previews. The input ZIP is explicitly tracked at repository `data/pc057-frozen-inputs-2026.10.28.zip`, so the supported workflow does not depend on an unfulfilled future download location.
