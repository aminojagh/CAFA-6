import os
import hashlib


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def sha256_file(path):
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()