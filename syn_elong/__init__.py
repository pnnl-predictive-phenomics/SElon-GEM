import cobra
import pathlib
import pandas as pd


path = pathlib.Path(__file__).parent

model = cobra.io.read_sbml_model(
    path.joinpath('syn_elong.xml').__str__()
)
ijb792 = cobra.io.load_json_model(path.joinpath('iJB792.json').__str__())
ims837 =cobra.io.load_json_model(path.joinpath('iMS837.json').__str__())

exp_file_path = path.joinpath('data', 'experiments.yml').__str__()
expected_metab = pd.read_csv(
    path.joinpath('data', 'excreted', 'metabolites.csv').__str__(),
).bigg_id