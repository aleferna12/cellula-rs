from pathlib import Path
import matplotlib
import colorir as cl

PROJECT_DIR = Path(__file__).resolve().parent.parent
cl.config.DEFAULT_PALETTES_DIR = PROJECT_DIR / "data" / "palettes"
CS = cl.Palette.load()
matplotlib.use("Agg")
