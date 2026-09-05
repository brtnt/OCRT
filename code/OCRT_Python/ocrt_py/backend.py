"""Array-backend selector: numpy (default) or CuPy (GPU).

Selection:
  - env OCRT_PY_GPU=1  -> import cupy as the array module
  - else               -> numpy

All batched LUT-mode code uses `from .backend import xp, to_np, GPU`.
Algorithms are backend-agnostic; only the array module differs, so the
numpy path remains the bit-reference and the CuPy path targets the
7-significant-digit gate (GPU FMA/reduction order differs from CPU BLAS).
"""
import os

GPU = os.environ.get('OCRT_PY_GPU', '0') == '1'

if GPU:
    import cupy as xp  # noqa: F401

    def to_np(a):
        import cupy
        return cupy.asnumpy(a)

    def to_xp(a):
        # numpy(또는 리스트)를 GPU 배열로. None은 그대로.
        return None if a is None else xp.asarray(a)

    def sync():
        import cupy
        cupy.cuda.Stream.null.synchronize()

    def erfc(a):
        from cupyx.scipy.special import erfc as _erfc
        return _erfc(a)
else:
    import numpy as xp  # noqa: F401

    def to_np(a):
        return a

    def to_xp(a):
        return a

    def sync():
        pass

    def erfc(a):
        from scipy.special import erfc as _erfc
        return _erfc(a)
