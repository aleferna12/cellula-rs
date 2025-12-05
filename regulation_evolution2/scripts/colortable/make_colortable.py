import argparse
from colorir import *
from scripts.colortable.print_colortable import read_colortable, print_colortable


# Make changes to this function to create different palettes
def get_parser():
    def run(args):
        pal = StackPalette.load(args.palettes)
        write_colortable(pal, args.outputfile)
        if args.print:
            print_colortable(args.outputfile)

    parser = argparse.ArgumentParser(
        description="Create a color table file from a set of predefined palettes"
    )
    parser.add_argument("outputfile", help="Output file for the color table")
    parser.add_argument(
        "palettes",
        default=["categ8", "chemgrad64_2", "miggrad8", "divgrad8", "tab20"],
        nargs="*",
        help="Palettes to add to the color table (in order). The palettes files "
             "must be located in 'data/palettes'. To generate palettes use the "
             "'colorir' python package (default: categ8 chemgrad_2 miggrad8 divgrad8)"
    )
    parser.add_argument("-p",
                        "--print",
                        help="Print the palette after writing it",
                        action="store_true")
    parser.set_defaults(run=run)
    return parser


def write_colortable(spalette: StackPalette, outfile):
    spalette.color_format = ColorFormat(sRGB, max_rgb=255, round_to=0)
    out_str = ""
    for i, color in enumerate(spalette):
        out_str += f"{i} {color}\n"
    with open(outfile, 'w') as file:
        file.write(out_str.replace('(', '').replace(')', '').replace(',', ''))
