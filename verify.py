"""Verify the preserved Ver.1.0 archive; run the original audit in a temporary directory."""
import argparse, hashlib, json, subprocess, sys, tempfile, zipfile
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--check-only',action='store_true'); args=p.parse_args()
archive=Path(__file__).resolve().parent/'release'/'Pentagonal_Dual_Crystal_Integrated_Proof_JA_v1.0_bundle.zip'
def sha(b): return hashlib.sha256(b).hexdigest()
assert sha(archive.read_bytes())=='4bd367aec3433a996ddf70402d8e6c56ab730553b9a4176dcd339cfbf251343f', 'Archive SHA-256 mismatch'
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    manifest=json.loads(z.read('MANIFEST_SHA256.json'))
    assert all(sha(z.read(n))==h for n,h in manifest.items())
    print(f'Archive: {len(manifest)} SHA-256 checks PASS',flush=True)
    if not args.check_only:
        with tempfile.TemporaryDirectory(prefix='dual_crystal_') as tmp:
            z.extractall(tmp)
            subprocess.run([sys.executable,str(Path(tmp)/'verify_v10.py')],cwd=tmp,check=True)
