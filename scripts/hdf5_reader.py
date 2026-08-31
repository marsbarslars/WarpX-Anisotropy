"""Print the structure of an HDF5 file, including groups, datasets, and optionally attributes."""

import h5py
from pathlib import Path

### USER-DEFINED VARIABLES ###
INPUT_PATH = Path('copy.h5').resolve()
PRINT_ATTRIBUTES = True


def print_hdf5_structure(path: Path):
    '''Print the structure of an HDF5 file, including groups, datasets, and optionally attributes.'''

    with h5py.File(path, "r") as file:

        def print_item(name, item):
            indent = "  " * name.count("/")

            if isinstance(item, h5py.Group):
                print(f"{indent}[Group]   {name}")

            elif isinstance(item, h5py.Dataset):
                print(
                    f"{indent}[Dataset] {name} "
                    f"shape={item.shape} "
                    f"dtype={item.dtype}"
                )

                if item.size <= 5:
                    print(f"{indent}           info: {item[...]}")
                else:
                    print(f"{indent}           preview: {item[...].flatten()[:5]}...")

            # Print attributes
            if PRINT_ATTRIBUTES == True:
                if item.attrs:
                    print(f"{indent}           @attributes:")
                    for key, value in item.attrs.items():
                        print(f"{indent}             {key} = {value}")

        file.visititems(print_item)


if __name__ == "__main__":

    print_hdf5_structure(INPUT_PATH)