"""Scale and translate an HDF5 file's grid and field data."""

from pathlib import Path
import shutil
import h5py
import numpy as np

### USER-DEFINED VARIABLES ###
INPUT_PATH = Path('example-femm-3d.h5').resolve()
OUTPUT_PATH = Path('copy.h5').resolve()
GRID_SCALE = 1.
GRID_TRANSLATION = np.array([0., 0., -2.5])
B_SCALE = 1.
E_SCALE = 1.


def change_item(name: str, item: h5py.Dataset | h5py.Group) -> None:
    '''Scale an HDF5 item's attributes and datasets.'''

    grid_attributes = {'gridSpacing', 'gridGlobalOffset'}

    for attribute_name, attribute_value in item.attrs.items():
        if attribute_name in grid_attributes and np.issubdtype(attribute_value.dtype, np.number):
            item.attrs[attribute_name] *= GRID_SCALE
            if attribute_name == 'gridGlobalOffset':
                item.attrs[attribute_name] += GRID_TRANSLATION

    if isinstance(item, h5py.Dataset):
        dataset_name = name.rsplit("/", maxsplit=2)
        if len(dataset_name) > 1:
            if dataset_name[-2] == 'B' and np.issubdtype(item.dtype, np.number):
                item[...] *= B_SCALE
            elif dataset_name[-2] == 'E' and np.issubdtype(item.dtype, np.number):
                item[...] *= E_SCALE


def main() -> None:

    # Throw errors
    if GRID_SCALE <= 0:
        raise ValueError("GRID_SCALE must be greater than zero")
    if B_SCALE <= 0:
        raise ValueError("B_SCALE must be greater than zero")
    if E_SCALE <= 0:
        raise ValueError("E_SCALE must be greater than zero")
    if INPUT_PATH == OUTPUT_PATH:
        raise ValueError("input and output paths must be different")
    if OUTPUT_PATH.exists():
        raise FileExistsError(f"output file already exists: {OUTPUT_PATH}")

    shutil.copy2(INPUT_PATH, OUTPUT_PATH)

    with h5py.File(OUTPUT_PATH, "r+") as file:
        file.visititems(change_item)

    print(f"Created: {OUTPUT_PATH}")
    print(f"Grid Scale Factor: {GRID_SCALE}")
    print(f"Grid Translation: {GRID_TRANSLATION}")
    print(f"B-Field Scale Factor: {B_SCALE}")
    print(f"E-Field Scale Factor: {E_SCALE}")

if __name__ == "__main__":
    main()