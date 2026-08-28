"""Scale the physical grid coordinates in an HDF5 magnetic-field file.

The magnetic-field values are not changed. OpenPMD grid metadata such as
``gridSpacing`` and ``gridGlobalOffset`` is multiplied by ``scale_factor``.
Explicit coordinate datasets named ``x``, ``y``, or ``z`` are scaled too.

Install the dependency with:
    python -m pip install h5py numpy

Example:
    python scale_hdf5_grid.py w7x_Bfield.h5 w7x_Bfield_scaled.h5 --scale-factor 0.001
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import h5py
import numpy as np


GRID_ATTRIBUTES = {"gridspacing", "gridglobaloffset"}
COORDINATE_DATASETS = {"x", "y", "z"}


def _scaled_attribute(value: object, scale_factor: float) -> object:
    """Return an HDF5 attribute value multiplied by the scale factor."""
    if isinstance(value, np.ndarray):
        return value * scale_factor
    if isinstance(value, np.number):
        return value * scale_factor
    if isinstance(value, (int, float, complex)):
        return value * scale_factor
    return value


def scale_hdf5_grid(
    input_path: Path,
    output_path: Path,
    scale_factor: float,
) -> None:
    """Write a copy of input_path with its physical grid scaled."""
    if scale_factor <= 0:
        raise ValueError("scale_factor must be greater than zero")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("input and output paths must be different")
    if output_path.exists():
        raise FileExistsError(f"output file already exists: {output_path}")

    shutil.copy2(input_path, output_path)

    with h5py.File(output_path, "r+") as handle:
        scaled_attributes = 0
        scaled_coordinate_datasets = 0

        def scale_item(name: str, item: h5py.Dataset | h5py.Group) -> None:
            nonlocal scaled_attributes, scaled_coordinate_datasets

            for attribute_name, attribute_value in item.attrs.items():
                if attribute_name.lower() in GRID_ATTRIBUTES:
                    item.attrs[attribute_name] = _scaled_attribute(
                        attribute_value, scale_factor
                    )
                    scaled_attributes += 1

            if isinstance(item, h5py.Dataset):
                dataset_name = name.rsplit("/", maxsplit=1)[-1].lower()
                if dataset_name in COORDINATE_DATASETS and np.issubdtype(
                    item.dtype, np.number
                ):
                    item[...] = item[...] * scale_factor
                    scaled_coordinate_datasets += 1

        handle.visititems(scale_item)

    print(f"Created: {output_path}")
    print(f"scale_factor: {scale_factor}")
    print(f"Scaled grid attributes: {scaled_attributes}")
    print(f"Scaled coordinate datasets: {scaled_coordinate_datasets}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path, help="Original HDF5 file")
    parser.add_argument("output_file", type=Path, help="Scaled HDF5 file to create")
    parser.add_argument(
        "--scale-factor",
        dest="scale_factor",
        type=float,
        required=True,
        help="Factor applied to physical grid coordinates and spacing",
    )
    args = parser.parse_args()

    if not args.input_file.is_file():
        parser.error(f"input file does not exist: {args.input_file}")

    try:
        scale_hdf5_grid(args.input_file, args.output_file, args.scale_factor)
    except (FileExistsError, OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
