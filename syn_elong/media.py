import pathlib
import pandas as pd

path = pathlib.Path(__file__).parent


m9_path = path.joinpath('data', 'media', 'min_media_all_wavelength.csv').__str__()
m9_media = pd.read_csv(m9_path)

m9_media = m9_media.set_index('exchange')['uptake'].to_dict()


syn_min_media_path = path.joinpath('data', 'media', 'syn_min_media.csv').__str__()
syn_min_media = pd.read_csv(syn_min_media_path)

syn_min_media = syn_min_media.set_index('exchange')['uptake'].to_dict()
