
import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1"
#import os
#os.environ['JAVA_HOME'] = r'C:\Program Files\Java\jdk-21'
import straindesign as sd
import cobra

import pandas as pd
pd.set_option('display.float_format', lambda x: f'{x:.3f}')
from syn_elong.strain_design_simplified import compute_strain_designs, StrainDesign
from syn_elong import model
from syn_elong.media import m9_media

# model.solver = 'glpk'

model.medium = m9_media
consistent_model = cobra.flux_analysis.fastcc(model)
consistent_model.medium = m9_media


rxn_cost = {}
for rxn in consistent_model.reactions:
    if rxn.id.startswith('EX_') or rxn.id.startswith('BIOMASS_') or\
            rxn.id.startswith('SK_') or rxn.id.startswith('PHOA') or\
            rxn.id.startswith('DM_') or rxn.id.startswith('SK_'):
        continue
    if rxn.gene_reaction_rule == '':
        continue
    else:
        rxn_cost[rxn.id] = 1
rxn_cost.pop('BCT1_syn')
rxn_cost.pop('ASPO6')
rxn_cost.pop('ASPO5')

import logging
logging.basicConfig(level=logging.INFO)
with consistent_model as m:
    m.reactions.get_by_id('EX_sucr_e').lower_bound = 0
    module_suppress = sd.SDModule(
        m,
        sd.names.SUPPRESS,
        constraints=['EX_sucr_e - 1 BIOMASS__1 <= 0',
                     # 'BIOMASS__1 >= 0.01'
                     ]
    )
    module_protect = sd.SDModule(
        m,
        sd.names.PROTECT,
        constraints='BIOMASS__1>=.1'
    )
    module_tilted_optknock = sd.SDModule(m, sd.names.OPTKNOCK,
                                         inner_objective='BIOMASS__1 ',  # - 0.001 EX_14bdo_e
                                         outer_objective='EX_sucr_e',
                                         constraints=['BIOMASS__1  >= 0.2', 'EX_sucr_e >=3'])
    module_optcouple = sd.SDModule(m, sd.names.OPTCOUPLE,
                                   inner_objective='BIOMASS__1',
                                   prod_id='EX_sucr_e',
                                   min_gcp=0.1)
    sd_helper = StrainDesign(
        m,
        # sd_modules=[module_suppress, module_protect],
        sd_modules=[module_optcouple],
        ko_cost=rxn_cost,
        # gene_kos=True,
    )

    sols = sd_helper.run(max_solutions=10, max_cost=5, time_limit=3000, solution_approach=sd.names.ANY)


print("Completed run")
for r in  sols.reaction_sd:
    print(r)

# Print solutions
print(f"One compressed solution with cost {sols.sd_cost[0]} found and "+\
      f"expanded to {len(sols.reaction_sd)} solutions in the uncompressed netork.")
print(f"Example knockout set: {[s for s in sols.reaction_sd[0]]}")



"""
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'NTRIRfx': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'NTRARf2': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'NO3abcpp': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'GND': -1.0, 'PSP_L': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'GND': -1.0, 'PGCD': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'GND': -1.0, 'PSERT': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'G6PDH2r': -1.0, 'PSP_L': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'G6PDH2r': -1.0, 'PGCD': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'G6PDH2r': -1.0, 'PSERT': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'PGL': -1.0, 'PSP_L': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'PGL': -1.0, 'PGCD': -1.0}
{'CYTBD4cm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'PGL': -1.0, 'PSERT': -1.0}
{'NDHPQRcm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'MNHNAtpp': -1.0, 'NTRIRfx': -1.0}
{'NDHPQRcm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'MNHNAtpp': -1.0, 'NTRARf2': -1.0}
{'NDHPQRcm': -1.0, 'BCT1_syn': -1.0, 'PCXHtpp': -1.0, 'MNHNAtpp': -1.0, 'NO3abcpp': -1.0}
"""