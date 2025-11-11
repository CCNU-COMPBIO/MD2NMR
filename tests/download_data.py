import os
import requests
from tqdm import tqdm
# Download the dataset from Zenodo
def download_dataset(url, filename):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size = int(response.headers.get('content-length', 0))  # Total file size in bytes
        block_size = 1024  # Bytes per chunk
        os.makedirs(os.path.dirname(filename), exist_ok=True)  # Ensure target directory exists

        # Display progress bar using tqdm
        with open(filename, 'wb') as f, tqdm(
            total=total_size,
            unit='B',
            unit_scale=True,
            desc=os.path.basename(filename),
            ncols=80
        ) as progress_bar:
            for data in response.iter_content(block_size):
                f.write(data)
                progress_bar.update(len(data))

    print(f"\n✅ Downloaded {filename} successfully.")
"""
# Define the dataset URL and filenames

"""
# Define the dataset URL and filenames
topology_url = "https://zenodo.org/records/16753027/files/H3.pdb"
trajectory_url = "https://zenodo.org/records/16753027/files/H3_1.xtc"
topology_filename = "./data/H3.pdb"
trajectory_filename = "./data/H3_1.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)


topology_url = "https://zenodo.org/records/16753027/files/WT_rw_run1.pdb"
trajectory_url = "https://zenodo.org/records/16753027/files/WT_rw_run1_2000ns_40ps.xtc"
topology_filename = "./data/WT_rw_run1.pdb"
trajectory_filename = "./data/WT_rw_run1_2000ns_40ps.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)
