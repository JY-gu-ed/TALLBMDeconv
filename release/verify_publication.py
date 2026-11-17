from pathlib import Path
import hashlib,tarfile
root=Path(__file__).resolve().parents[1]
for line in (root/'release/INFERENCE_SHA256SUMS.txt').read_text().splitlines():
 h,p=line.split('  ',1)
 assert hashlib.sha256((root/p).read_bytes()).hexdigest()==h,p
for p in root.rglob('*'):
 if not p.is_file() or '.git' in p.parts: continue
 assert p.suffix.lower() not in {'.rds','.fcs','.xlsx','.h5','.h5ad','.bam','.fastq'},str(p)
with tarfile.open(root/'dist/TALLBMDeconv_0.1.1.tar.gz') as t:
 for m in t.getmembers():
  if m.name.lower().endswith(('.rds','.fcs','.xlsx','.h5','.h5ad')):
   # R-generated vignette index, not a reference/data matrix.
   assert m.name=='TALLBMDeconv/build/vignette.rds' and m.size==233,m.name
print('PASS_CODE_ONLY_AND_UNCHANGED_INFERENCE')
