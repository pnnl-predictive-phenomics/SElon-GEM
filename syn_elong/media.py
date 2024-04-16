import pathlib
import pandas as pd

path = pathlib.Path(__file__).parent


m9_path = path.joinpath('data', 'media', 'min_media_all_wavelength.csv').__str__()
min_media = pd.read_csv(m9_path)

min_media = min_media.set_index('exchange')['uptake'].to_dict()
