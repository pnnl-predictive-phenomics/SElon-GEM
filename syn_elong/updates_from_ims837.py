"""
File adds reactions from the IMS837 model.
The actual IMS837 model doesn't work with memote due to sbml compliance. This file attempts to recreate
the reactions/changes from the IMS837 model, building off the current model. The code is copied from the ims837 repo.
There were 4 reactions missing that were assumed to exist. For those 4, we grabbed them from the universal model
or from the ijb792 model.


"""
from cobra import Reaction, Metabolite
from concerto.utils import universal_model
from syn_elong import ijb792

def add_reactions_not_in_model_but_ims837_thinks_exists(model):
    # reaction not in the model currently, so adding it
    rxns_to_add = []
    rxns_in_universal_model = ['FNOR', 'PC6YM']
    for rxn in rxns_in_universal_model:
        rxns_to_add.append(universal_model.reactions.get_by_id(rxn).copy())


    rxns_from_ijb792 = ['LYCBC1', 'MOGDS']
    for rxn in rxns_from_ijb792:
        rxns_to_add.append(ijb792.reactions.get_by_id(rxn).copy())

    model.add_reactions(rxns_to_add)
    return model


def update_5(model):
    model = add_reactions_not_in_model_but_ims837_thinks_exists(model)
    # add gene Synpcc7942_B2620 to CAT reaction
    model.reactions.get_by_id('CAT').gene_reaction_rule = '( Synpcc7942_1656 or Synpcc7942_B2620 )'

    # where does this logic come from?
    # add genes Synpcc7942_1088, Sypncc7942_0239, Synpcc7942_1630, Synpcc7942_2542, Synpcc7942_1054 and Synpcc7942_2158 to CBFCum reaction
    model.reactions.get_by_id('CBFCum').gene_reaction_rule = \
        '(Synpcc7942_1231 and Synpcc7942_2331 and Synpcc7942_1232 and Synpcc7942_2332 and Synpcc7942_1088 ' \
        'and Synpcc7942_1479 and ((Synpcc7942_0239) or (Synpcc7942_1630) or (Synpcc7942_2542)) and ' \
        'Synpcc7942_0113 and Synpcc7942_2426 and Synpcc7942_0475 and Synpcc7942_1053 and Synpcc7942_1048 and' \
        ' Synpcc7942_1052 and Synpcc7942_1047 and Synpcc7942_1049 and Synpcc7942_1050 and Synpcc7942_1051 and ' \
        'Synpcc7942_1054 and Synpcc7942_1055 and Synpcc7942_2030 and Synpcc7942_0327 and Synpcc7942_0326 and ' \
        'Synpcc7942_0325 and Synpcc7942_0328 and Synpcc7942_2158 and Synpcc7942_0322 and Synpcc7942_0323 and ' \
        'Synpcc7942_1478 and Synpcc7942_0511)'
    # add genes Synpcc7942_0201 and Sypncc7942_0202 to CYOOum reaction
    model.reactions.get_by_id('CYOOum').gene_reaction_rule = \
        '(Synpcc7942_2602 or Synpcc7942_0201) and (Synpcc7942_2603 or Synpcc7942_0202) and Synpcc7942_2604'
    # add gene Synpcc7942_1408 to FE3abcppreaction
    model.reactions.get_by_id('FE3abcpp').gene_reaction_rule = \
        '(Synpcc7942_1409 and Synpcc7942_1407 and Synpcc7942_1406 and Synpcc7942_1408) or (Synpcc7942_2175 ' \
        'and Synpcc7942_1407 and Synpcc7942_1406 and Synpcc7942_1408)'



    # add genes Synpcc7942_0338, Synpcc7942_0698, Synpcc7942_0898, Synpcc7942_1749, Synpcc7942_2581, Synpcc7942_1499,
    # Synpcc7942_0814 and Synpcc7942_1541 to FNOR reaction
    model.reactions.get_by_id('FNOR').gene_reaction_rule = \
        '((Synpcc7942_0338 or Synpcc7942_0698  or Synpcc7942_0898 or Synpcc7942_1749) and Synpcc7942_2581 and ' \
        'Synpcc7942_1499  and Synpcc7942_0814 and Synpcc7942_1541 and Synpcc7942_0978)'
    # add gene Synpcc7942_0508 to GGDPR reaction
    model.reactions.get_by_id('GGDPR').gene_reaction_rule = '(Synpcc7942_0385) or (Synpcc7942_0508)'
    # add gene Synpcc7942_0357 to H2CO3_NAt_syn reaction
    model.reactions.get_by_id('H2CO3_NAt_syn').gene_reaction_rule = '(Synpcc7942_1475) or (Synpcc7942_0357)'
    # add genes Synpcc7942_1388 and Synpcc7942_B2619 to HCO3E_1_cx reaction
    model.reactions.get_by_id('HCO3E_1_cx').gene_reaction_rule = \
        '(Synpcc7942_1447) or (Synpcc7942_1388) or (Synpcc7942_B2619) or (Synpcc7942_1423)'
    # add gene Synpcc7942_1378 to LPADSS2 reaction
    model.reactions.get_by_id('LPADSS2').gene_reaction_rule = '(Synpcc7942_0932) and (Synpcc7942_1378)'
    # add gene Synpcc7942_0652 to LYCBC1 reaction
    model.reactions.get_by_id('LYCBC1').gene_reaction_rule = '(Synpcc7942_2062) or (Synpcc7942_0652)'
    # add genes Synpcc7942_0888 and Synpcc7942_0498 to MAN1PT reaction
    model.reactions.get_by_id('MAN1PT').gene_reaction_rule = \
        '(Synpcc7942_1608) or (Synpcc7942_0888) or (Synpcc7942_0498) or (Synpcc7942_1973)'
    # add gene Synpcc7942_1763 to MI3PP reaction
    model.reactions.get_by_id('MI3PP').gene_reaction_rule = '(Synpcc7942_2582) or (Synpcc7942_1763)'
    # add gene Synpcc7942_1910 to NPHBDC reaction
    model.reactions.get_by_id('NPHBDC').gene_reaction_rule = \
        '(Synpcc7942_2588 and Synpcc7942_2055 and Synpcc7942_1910 and Synpcc7942_0135)'
    # add gene Synpcc7942_1266 to NTPP2 reaction
    model.reactions.get_by_id('NTPP2').gene_reaction_rule = '(Synpcc7942_1266) or (Synpcc7942_1493)'
    # add gene Synpcc7942_0706 to PC6YM reaction
    model.reactions.get_by_id('PC6YM').gene_reaction_rule = '(Synpcc7942_0706) or (Synpcc7942_1850)'
    # add gene Synpcc7942_0693 to PGLYCP reaction
    model.reactions.get_by_id('PGLYCP').gene_reaction_rule = '(Synpcc7942_2613) or (Synpcc7942_0693)'
    # add gene Synpcc7942_1450 to PIuabcpp reaction
    model.reactions.get_by_id('PIuabcpp').gene_reaction_rule = \
        '((Synpcc7942_2444) or (Synpcc7942_2445)) and (Synpcc7942_2443) and (Synpcc7942_2442) ' \
        'and (Synpcc7942_1450 or Synpcc7942_2441)'
    # add gene Synpcc7942_0974 to UDPGD reaction
    model.reactions.get_by_id('UDPGD').gene_reaction_rule = '(Synpcc7942_0973 or Synpcc7942_0974)'
    # add gene Synpcc7942_1319 to UPPRT reaction
    model.reactions.get_by_id('UPPRT').gene_reaction_rule = '( Synpcc7942_1715 or Synpcc7942_1319 )'
    # add gene Synpcc7942_B2633 to MOGDS reaction
    model.reactions.get_by_id('MOGDS').gene_reaction_rule = '( Synpcc7942_1189 or Synpcc7942_1211 or Synpcc7942_B2633 )'
    # add gene Synpcc7942_1516 to PGM reaction
    model.reactions.get_by_id('PGM').gene_reaction_rule = '( Synpcc7942_0485 or Synpcc7942_0469 or Synpcc7942_1516 or Synpcc7942_2078 )'

    m1 = model.metabolites.get_by_id('h2o_c')  # h2o
    m2 = Metabolite(
        'acmama_c',
        formula='C14H23N2O9',
        name='N-Acetyl-D-muramoyl-L-alanine',
        compartment='c')  # N-Acetyl-D-muramoyl-L-alanine
    m3 = model.metabolites.get_by_id('ala__L_c')  # L-alanine
    m4 = Metabolite(
        'acmam_c',
        formula='C11H18NO8',
        name='N-Acetyl-D-muramoate',
        compartment='c')  # N-Acetyl-D-muramoate


    # add AMAA reaction:
    reaction = Reaction('AMAA')
    reaction.name = 'N-acetylmuramoyl-L-alanine amidase'
    reaction.subsystem = 'Peptidoglycan Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: - 1.0, m3: 1.0, m4: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('h_c')  # H+
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = Metabolite(
        'csn_c',
        formula='C4H5N3O',
        name='Cytosine',
        compartment='c')  # Cytosine
    m4 = model.metabolites.get_by_id('nh4_c')  # NH4
    m5 = model.metabolites.get_by_id('ura_c')  # uracil

    # add CSND reaction:
    reaction = Reaction('CSND')
    reaction.name = 'Cytosine deaminase'
    reaction.subsystem = 'Pyrimidine metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: - 1.0, m3: -1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])
    m1 = model.metabolites.get_by_id('h2o_c')  # h2o
    m2 = Metabolite(
        'cyst-L_c',
        formula='C7H14N2O4S',
        name='L-Cystathionine',
        compartment='c')  # L-Cystathionine
    m3 = model.metabolites.get_by_id('pyr_c')  # pyruvate
    m4 = model.metabolites.get_by_id('nh4_c')  # NH4
    m5 = model.metabolites.get_by_id('hcys__L_c')  # L-Homocysteine

    # add CYSTL reaction:
    reaction = Reaction('CYSTL')
    reaction.name = 'cystathionine b-lyase'
    reaction.subsystem = 'Others'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: - 1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])
    m1 = model.metabolites.get_by_id('13dpg_c')  # 3-Phospho-D-glyceroyl phosphate
    m2 = model.metabolites.get_by_id('h_c')  # H+
    m3 = Metabolite(
        '23dpg_c',
        formula='C3H3O10P2',
        name='2,3-Disphospho-D-glycerate',
        compartment='c')  # 2,3-Disphospho-D-glycerate


    # add DPGM reaction:
    reaction = Reaction('DPGM')
    reaction.name = 'Diphosphoglyceromutase'
    reaction.subsystem = 'Glycolysis/Gluconeogenesis'
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: 1.0, m3: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('5drib_c')  # 5-Deoxy-D-ribose
    m3 = model.metabolites.get_by_id('adp_c')  # ADP
    m4 = model.metabolites.get_by_id('h_c')  # H+
    m5 = Metabolite(
        '2dr5p_c',
        formula='C5H9O7P',
        name='2-Deoxy-D-ribose 5-phosphate',
        compartment='c')  # 2-Deoxy-D-ribose 5-phosphate


    # add DRBK reaction:
    reaction = Reaction('DRBK')
    reaction.name = 'Deoxyribokinase'
    reaction.subsystem = 'Pentose phosphate pathway'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('2dr5p_c')  # 2-Deoxy-D-ribose 5-phosphate
    m2 = model.metabolites.get_by_id('g3p_c')  # D-Glyceraldehyde 3-phosphate
    m3 = model.metabolites.get_by_id('acald_c')  # Acetaldehyde


    # add DRPA reaction:
    reaction = Reaction('DRPA')
    reaction.name = 'deoxyribose-phosphate aldolase'
    reaction.subsystem = 'Pentose phosphate pathway'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: 1.0, m3: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('nad_c')  # NAD+
    m2 = model.metabolites.get_by_id('glc__D_c')  # D-glucose
    m3 = model.metabolites.get_by_id('h_c')  # H+
    m4 = model.metabolites.get_by_id('nadh_c')  # NADH
    m5 = Metabolite(
        'g15lac_c',
        formula='C6H10O6',
        name='D-Glucono-1,5-lactone',
        compartment='c')  # D-Glucono-1,5-lactone


    # add G1Dxreaction:
    reaction = Reaction('G1Dx')
    reaction.name = 'Glucose 1 dehydrogenase NAD'
    reaction.subsystem = 'Carbohydrates and related molecules'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('gtp_c')  # GTP
    m2 = model.metabolites.get_by_id('ppi_c')  # Diphosphate
    m3 = Metabolite(
        '35cgmp_c',
        formula='C10H11N5O7P',
        name='3 5 Cyclic GMP',
        compartment='c')  # 3',5'-Cyclic GMP


    # add GUACYC reaction:
    reaction = Reaction('GUACYC')
    reaction.name = 'Diguanylate cyclase'
    reaction.subsystem = 'Purine Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: 1.0, m3: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('nadh_c')  # NADH
    m2 = Metabolite(
        'rdxo_c',
        formula='Fe1SO',
        name='rubredoxin_oxidized',
        compartment='c')  # rubredoxin_oxidized
    m3 = model.metabolites.get_by_id('h_c')  # H+
    m4 = model.metabolites.get_by_id('nad_c')  # NAD+
    m5 = Metabolite(
        'rdxr_c',
        formula='Fe1SO',
        name='rubredoxin_reduced',
        compartment='c')  # rubredoxin_reduced

    # add RDXRr reaction:
    reaction = Reaction('RDXRr')
    reaction.name = 'Rubrerythrin'
    reaction.subsystem = 'Unassigned'
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -2.0, m3: 1.0, m4: 1.0, m5: 2.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('dcaACP_c')  # Decanoyl-[acyl-carrier protein]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'dca_c',
        formula='C10H19O2',
        name='Decanoate__n_C100_',
        compartment='c')  # Decanoate__n_C100_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA100ACPHi reaction:
    reaction = Reaction('FA100ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('ddcaACP_c')  # Dodecanoyl-[acyl-carrier protein]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'ddca_c',
        formula='C12H23O2',
        name='Dodecanoate__n_C120_',
        compartment='c')  # Dodecanoate__n_C120_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA120ACPHi reaction:
    reaction = Reaction('FA120ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('myrsACP_c')  # Tetradecanoyl-[acyl-carrier protein]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'ttdca_c',
        formula='C14H27O2',
        name='tetradecanoate__n_C140_',
        compartment='c')  # tetradecanoate__n_C140_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA140ACPHi reaction:
    reaction = Reaction('FA140ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('palmACP_c')  # (9Z)-Hexadecanoyl-[acp]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'hdca_c',
        formula='C16H31O2',
        name='Hexadecanoate__n_C160_',
        compartment='c')  # Hexadecanoate__n_C160_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA160ACPHi reaction:
    reaction = Reaction('FA160ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('hdeACP_c')  # Hexadecenoyl-[acyl-carrier protein]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'hdcea_c',
        formula='C16H29O2',
        name='Hexadecenoate__n_C161_',
        compartment='c')  # Hexadecenoate__n_C161_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA161ACPHi reaction:
    reaction = Reaction('FA161ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('3ooctdACP_c')  # 3-Oxostearoyl-[acp]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = model.metabolites.get_by_id('ocdca_c')  # Octadecanoate (n-C18:0)
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA180ACPHi reaction:
    reaction = Reaction('FA180ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('octe9ACP_c')  # Oleoyl-[acyl-carrier protein]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'octe9_c',
        formula='C18H32O',
        name='Oleoyl',
        compartment='c')  # Oleoyl
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA181ACPHi reaction:
    reaction = Reaction('FA181ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('ocACP_c')  # Octanoyl-[acp]
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('ACP_c')  # Acyl-carrier protein
    m4 = Metabolite(
        'octa_c',
        formula='C8H15O2',
        name='octanoate__n_C80_',
        compartment='c')  # octanoate__n_C80_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FA80ACPHi reaction:
    reaction = Reaction('FA80ACPHi')
    reaction.name = 'fatty acyl ACP hydrolase'
    reaction.subsystem = 'Fatty Acid Biosynthesis'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'dcacoa_c',
        formula='C31H50N7O17P3S',
        name='Decanoyl_CoA__n_C100CoA_',
        compartment='c')  # Decanoyl_CoA__n_C100CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('dca_c')  # Decanoate__n_C100_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE100 reaction:
    reaction = Reaction('FACOAE100')
    reaction.name = 'fatty acid CoA thioesterase decanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'ddcacoa_c',
        formula='C33H54N7O17P3S',
        name='Dodecanoyl_CoA__n_C120CoA_',
        compartment='c')  # Dodecanoyl_CoA__n_C120CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('ddca_c')  # Dodecanoate__n_C120_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE120 reaction:
    reaction = Reaction('FACOAE120')
    reaction.name = 'fatty acid CoA thioesterase dodecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'tdcoa_c',
        formula='C35H58N7O17P3S',
        name='Tetradecanoyl_CoA__n_C140CoA_',
        compartment='c')  # Tetradecanoyl_CoA__n_C140CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('ttdca_c')  # tetradecanoate__n_C140_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE140 reaction:
    reaction = Reaction('FACOAE140')
    reaction.name = 'fatty acid CoA thioesterase tetradecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'pmtcoa_c',
        formula='C37H62N7O17P3S',
        name='Palmitoyl_CoA__n_C160CoA_',
        compartment='c')  # Palmitoyl_CoA__n_C160CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('hdca_c')  # Hexadecanoate__n_C160_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE160 reaction:
    reaction = Reaction('FACOAE160')
    reaction.name = 'fatty acid CoA thioesterase hexadecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'hdcoa_c',
        formula='C37H60N7O17P3S',
        name='Hexadecenoyl_CoA__n_C161CoA_',
        compartment='c')  # Hexadecenoyl_CoA__n_C161CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('hdcea_c')  # Hexadecenoate__n_C161_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE161 reaction:
    reaction = Reaction('FACOAE161')
    reaction.name = 'fatty acid CoA thioesterase hexadecenoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'stcoa_c',
        formula='C39H66N7O17P3S',
        name='Stearoyl_CoA__n_C180CoA_',
        compartment='c')  # Stearoyl_CoA__n_C180CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('ocdca_c')  # Octadecanoate (n-C18:0)
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE180 reaction:
    reaction = Reaction('FACOAE180')
    reaction.name = 'fatty acid CoA thioesterase octadecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'octe9coa_c',
        formula='C39H64N7O17P3S',
        name='Octadecenoyl_CoA__n_C181CoA_',
        compartment='c')  # Octadecenoyl_CoA__n_C181CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('octe9_c')  # Oleoyl
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE181 reaction:
    reaction = Reaction('FACOAE181')
    reaction.name = 'fatty acid CoA thioesterase octadecenoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'occoa_c',
        formula='C29H46N7O17P3S',
        name='Octanoyl_CoA__n_C80CoA_',
        compartment='c')  # Octanoyl_CoA__n_C80CoA_
    m2 = model.metabolites.get_by_id('h2o_c')  # h2o
    m3 = model.metabolites.get_by_id('coa_c')  # CoA
    m4 = model.metabolites.get_by_id('octa_c')  # octanoate__n_C80_
    m5 = model.metabolites.get_by_id('h_c')  # H+

    # add FACOAE80 reaction:
    reaction = Reaction('FACOAE80')
    reaction.name = 'fatty acid CoA thioesterase octanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: 1.0, m4: 1.0, m5: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('dca_c')  # Decanoate__n_C100_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('dcacoa_c')  # Decanoyl_CoA__n_C100CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL100i reaction:
    reaction = Reaction('FACOAL100i')
    reaction.name = 'fatty acid CoA ligase decanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('ddca_c')  # Dodecanoate__n_C120_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('ddcacoa_c')  # Dodecanoyl_CoA__n_C120CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL120i reaction:
    reaction = Reaction('FACOAL120i')
    reaction.name = 'fatty acid CoA ligase dodecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('ttdca_c')  # tetradecanoate__n_C140_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('tdcoa_c')  # Tetradecanoyl_CoA__n_C140CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL140i reaction:
    reaction = Reaction('FACOAL140i')
    reaction.name = 'fatty acid CoA ligase tetradecanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('hdca_c')  # Hexadecanoate__n_C160_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('pmtcoa_c')  # Palmitoyl_CoA__n_C160CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL160i reaction:
    reaction = Reaction('FACOAL160i')
    reaction.name = 'C160 fatty acid activation'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('hdcea_c')  # Hexadecenoate__n_C161_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('hdcoa_c')  # Hexadecenoyl_CoA__n_C161CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL161i reaction:
    reaction = Reaction('FACOAL161i')
    reaction.name = 'fatty acid CoA ligase hexadecenoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('ocdca_c')  # Octadecanoate (n-C18:0)
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('stcoa_c')  # Stearoyl_CoA__n_C180CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL180i reaction:
    reaction = Reaction('FACOAL180i')
    reaction.name = 'C180 fatty acid activation'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('octe9_c')  # Oleoyl
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('octe9coa_c')  # Octadecenoyl_CoA__n_C181CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL181i reaction:
    reaction = Reaction('FACOAL181i')
    reaction.name = 'C181 fatty acid activation'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = model.metabolites.get_by_id('atp_c')  # ATP
    m2 = model.metabolites.get_by_id('coa_c')  # CoA
    m3 = model.metabolites.get_by_id('octa_c')  # octanoate__n_C80_
    m4 = model.metabolites.get_by_id('amp_c')  # AMP
    m5 = model.metabolites.get_by_id('occoa_c')  # Octanoyl_CoA__n_C80CoA_
    m6 = model.metabolites.get_by_id('ppi_c')  # Diphosphate

    # add FACOAL80i reaction:
    reaction = Reaction('FACOAL80i')
    reaction.name = 'fatty acid CoA ligase octanoate'
    reaction.subsystem = 'Fatty Acid Metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -1.0, m3: -1.0, m4: 1.0, m5: 1.0, m6: 1.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'bglyg4n_c',
        formula='C24H40O20',
        name='branching_glycogen__4_units_',
        compartment='c')  # branching_glycogen__4_units_
    m2 = model.metabolites.get_by_id('pi_c')  # Orthophosphate
    m3 = model.metabolites.get_by_id('g1p_c')  # D-Glucose 1-phosphate

    # add GLCP4 reaction:
    reaction = Reaction('GLCP4')
    reaction.name = 'glycogen phosphorylase 4 units'
    reaction.subsystem = 'Glycogen and sucrose metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -4.0, m3: 4.0})
    model.add_reactions([reaction])

    m1 = Metabolite(
        'glyg4n_c',
        formula='C24H40O20',
        name='glycogen__4_units__linear_glucan',
        compartment='c')  # glycogen__4_units__linear_glucan
    m2 = model.metabolites.get_by_id('pi_c')  # Orthophosphate
    m3 = model.metabolites.get_by_id('g1p_c')  # D-Glucose 1-phosphate

    # add GLCP3 reaction:
    reaction = Reaction('GLCP3')
    reaction.name = 'glycogen phosphorylase 4 units'
    reaction.subsystem = 'Glycogen and sucrose metabolism'
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({m1: - 1.0, m2: -4.0, m3: 4.0})
    model.add_reactions([reaction])
    return model