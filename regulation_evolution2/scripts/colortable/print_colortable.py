import argparse
from colorir import *


def get_parser():
    def run(args):
        print_colortable(args.filepath)

    parser = argparse.ArgumentParser(description="Print a color table file to the terminal")
    parser.add_argument("filepath")
    parser.set_defaults(run=run)
    return parser


def read_colortable(filepath) -> StackPalette:
    with open(filepath) as file:
        pal_raw = file.read()
    pal = StackPalette()
    color_fmt = ColorFormat(sRGB, max_rgb=255)
    for row in pal_raw.split("\n"):
        if not row:
            continue
        color = tuple(int(x) for x in row.split()[1:])
        pal.add(color_fmt.format(color))
    return pal


def print_colortable(path):
    pal = read_colortable(path)
    palstr = swatch(pal, file=None)
    to_sub = palstr.count("\n")
    n_size = len(str(to_sub - 1))
    palstr = "0 " + " " * (n_size - 1) + palstr
    for i in range(to_sub):
        palstr = palstr.replace("\n", f"@{i + 1} {' ' * (n_size - len(str(i + 1)))}", 1)
    print(palstr.replace("@", "\n"))