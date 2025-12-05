import argparse
import logging
import skimage
import numpy as np
import pandas as pd
from enlighten import Counter
from scripts.fileio import *

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        templatedf = parse_cell_data(args.templatefile)
        celldf, latdf = make_competition(templatedf, args.imgfile, args.size, args.cell_length)
        logger.info("Writing output competition files")
        write_cell_data(celldf, args.outcellfile)
        write_lattice(latdf, args.latticefile)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Create cell and lattice CSV files to start competition experiments."
    )
    input_ = parser.add_argument_group("input")
    input_.add_argument(
        "templatefile",
        help="CSV file containing the template cells. The number of different groups must match "
             "the number of unique pixel colors in 'imgfile'. If more than one template has the "
             "same group attribute, they will be sampled sequentially when projected onto the "
             "lattice. These files can be generated with the 'make_templates' script"
    )
    input_.add_argument(
        "imgfile",
        help="PNG file containing a drawing of the initial setting of the simulation."
             "Colors represent entries on the 'templatefile' table (ordered in reverse by RGB value)."
             "The background must be white and will be ignored. Be careful not to have any fading "
             "colors, as these will be interpreted as additional cell templates. The image "
             "dimensions should match those of the simulation lattice"
    )
    output = parser.add_argument_group("output")
    output.add_argument("outcellfile", help="CSV output file containing cell data")
    output.add_argument("latticefile", help="CSV output lattice file")
    parser.add_argument(
        "-s",
        "--size",
        help="Size of the axis of the final lattice (default: size of the provided image)",
        default=None,
        type=int
    )
    parser.add_argument(
        "-c",
        "--cell_length",
        help="Diameter of the initialized cells (default: %(default)s)",
        default=7,
        type=int
    )
    parser.set_defaults(run=run)
    return parser


def make_competition(templatedf, imgfile, lattsize=None, cell_length=7):
    logger.info("Making competition file from image and templates")

    if imgfile[-3:].lower() != "png":
        raise ValueError("input image must be a PNG file")
    img = skimage.io.imread(imgfile)
    if img.shape[0] != img.shape[1]:
        raise ValueError("input image is not a square")

    if lattsize is None:
        lattsize = img.shape[0]
    img = skimage.transform.resize(
        img,
        (lattsize, lattsize),
        order=0,
        anti_aliasing=False
    )

    templatedf = templatedf.reset_index(drop=True)
    templatedf["time"] = 0
    gdf = templatedf.groupby("group")

    colors, indexes = np.unique(np.reshape(img, (-1, img.shape[2])), axis=0, return_inverse=True)
    if len(colors) - 1 != len(gdf):
        raise ValueError("number of non-white pixel colors in 'imgfile' and unique 'group' "
                         f"attributes in 'templatefile' must match (got {len(colors) - 1} and "
                         f"{len(gdf)} respectively)")

    # Reverse colors and indexes so they are ordered by reverse rgb when matching the cells
    colors = colors[::-1]
    indexes = np.reshape(indexes, (lattsize, lattsize))
    indexes = len(colors) - 1 - indexes

    # Basically we sample cell_length-interspaced points from image
    group_map = indexes[::cell_length, ::cell_length]
    # Then create a unique sigma for each point where a non-zero value was found in the image
    sigma_lat = np.arange(1, group_map.shape[0] ** 2 + 1).reshape(group_map.shape)
    sigma_count = np.unique(np.where(group_map, sigma_lat, 0), return_inverse=True)
    sigma_map = sigma_count[1].reshape(sigma_lat.shape)
    # Then tile the pattern to get back to lattice size
    lat = sigma_map.repeat(cell_length, axis=0).repeat(cell_length, axis=1)
    if lattsize % cell_length:
        margin = cell_length - (lattsize % cell_length)
        lat = lat[:-margin, :-margin]

    cells = []
    pbar = Counter(desc="Cells created", total=len(sigma_count[0]) - 1)
    for sigma, group in np.stack([sigma_map, group_map], axis=2).reshape((-1, 2)):
        if sigma != 0:
            groupdf = gdf.get_group(group - 1)
            cell_attrs = groupdf.iloc[sigma % len(groupdf)].copy()
            cell_attrs["sigma"] = sigma
            cell_attrs["ancestor"] = sigma
            cell_attrs["group"] = group - 1
            cells.append(cell_attrs)
            pbar.update()

    celldf = pd.DataFrame(cells).set_index("sigma", drop=False).sort_index()
    celldf.index.name = "sigma_i"

    # Remove borders
    trimmed_lat = lat[1:-1, 1:-1]
    latdf = pd.DataFrame(trimmed_lat)

    return celldf, latdf
