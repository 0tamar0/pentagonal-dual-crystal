# Pentagonal Dual Crystal (Ver.1.0)

[日本語](README.md)

An independent, not externally peer-reviewed project on recursively constructed regular-pentagon crystals A and B and their gap components a and b. For the specified coordinate model and prescribed recursive atomic decomposition, the computer-assisted proof establishes non-overlap, edge-to-edge or vertex-to-vertex tile contacts, and ten congruence classes of contacting atomic component pairs at all finite generations.

- [Illustrated guide (Japanese)](docs/Pentagonal_Dual_Crystal_Illustrated_Guide_JA_v1.0.pdf)
- [Paper (Japanese)](docs/Pentagonal_Dual_Crystal_Public_Paper_JA_v1.0.pdf)
- [Verification appendix (Japanese)](docs/Pentagonal_Dual_Crystal_Verification_Appendix_JA_v1.0.pdf)
- [Certificates](certificates/) and [source for inspection](proof/)

Run `python3 verify.py` from the repository root with Python 3 (standard library only). It checks the release archive and runs the original verification in a temporary directory. `--check-only` checks hashes without running numerical certificates. The complete reproducible environment is in the release archive; the extracted source directory is for inspection. The proof uses finite-state closure and strong induction, not exhaustive generation of 1,000 levels.

The theorem does not address arbitrary placements, infinite limits, plane tilings, contour connectivity, or minimality of the state count. No external peer review or proof-assistant verification has been performed. ChatGPT was used substantially in proof construction, coding, testing and documentation. Corrections and counterexamples are welcome. No public license has been selected.
