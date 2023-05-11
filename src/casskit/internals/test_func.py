import requests

from . import _cache_manager, _logger


def fetch_dataset(dataset_url):    
    # Assuming that the URL is the unique identifier for the dataset
    dataset_name = dataset_url.split("/")[-1]
    cache_path = _cache_manager.cache_dir / dataset_name

    if cache_path.exists():
        with open(cache_path, "r") as f:
            data = f.read()
    else:
        # Fetch the data
        response = requests.get(dataset_url)
        response.raise_for_status()

        data = response.text

        # Cache the data
        with open(cache_path, "w") as f:
            f.write(data)

        _logger.info(f"Fetched dataset {dataset_name} from {dataset_url}")

    return data
