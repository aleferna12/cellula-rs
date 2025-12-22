from pathlib import Path as _Path
if len(list(_Path(__file__).parent.glob("rust.cp*"))) > 0:
    from . import rust
del _Path

from . import calculate_adh
