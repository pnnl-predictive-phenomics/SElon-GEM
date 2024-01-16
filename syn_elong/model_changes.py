"""
Synnoch elongatus model.


Unmodified starting model is iJB785.xml.
Living model is syn_elong.xml


Adding on biolog exchange reactions to model to accurate reflect biolog tests.

"""
from memote.suite.cli.reports import diff
import cobra
import os
import logging

from concerto.utils.biolog_help import add_biolog_exchanges, universal_model

log = logging.getLogger()

_file_path = os.path.dirname(__file__)
starting_model_f_name = 'iJB785.xml'
s_model_path = os.path.join(_file_path, starting_model_f_name)

starting_model = cobra.io.read_sbml_model(s_model_path)
# starting_model = cobra.io.load_json_model('iJB792.json')
starting_model.id = "syn_elong"

output_model_name = 'syn_elong.xml'
output_model_path = os.path.join(_file_path, output_model_name)


def write_model(model):
    cobra.io.write_sbml_model(model, output_model_path)


def update_1(model):
    # add missing biolog reactions to model
    log.info("Adding RT to prefix")
    model = add_biolog_exchanges(model)
    return model


def update_2(model):
    # Add the sucrose transporter SUCRt2 and associated gene
    log.info("Adding Sucrose Transporter")

    # Copy sucrose transport from universal model
    sucr_transport = universal_model.reactions.SUCRt2.copy()
    sucr_transport.lower_bound = -1000
    # Add the transport reaction to the model
    model.add_reactions([sucr_transport])
    # Set the lower bounds for the sucrose reactions
    model.reactions.EX_sucr_e.lower_bound = 0.011
    # Add the gene name of the sucrose transporter to the model
    gene_add = cobra.core.Gene(id='cscB', name='cscB', functional=True)
    model.genes.add(gene_add)

    model.reactions.SUCRt2.gene_reaction_rule = '( cscB )'

    # upper bound of biomass shoul dnot be limited
    model.reactions.BIOMASS__1.upper_bound = 1000

    return model


def update_3(model):
    # updates the model to add the reactions from ijb792
    rxns = ['MPTSS', 'MOADSUx', 'GTPC', 'CPMPS', 'MPTS', 'MPTAT', 'MOCOS', 'MDH']
    reactions_to_add = set()
    for reaction_id in rxns:
        reaction = universal_model.reactions.get_by_id(reaction_id)
        reactions_to_add.add(reaction.copy())
    model.add_reactions(reactions_to_add)
    model.remove_reactions([universal_model.reactions.get_by_id('ORNTA')])
    return model


def update_4(model):
    # adding reactions from gapfilling
    reaction_ids = {
        'GLCS1', 'SADT2', 'ATPSu', 'PPND2', 'aratyr1', 'CDPDAGS1819Z160',
        'ASNTRS_1', 'GLGB', 'FRDO', 'ACGAM6PS', 'PGPS1819Z160', 'PAPSSH', 'GLYOX_2',
        'APSR', 'AIRC1', 'ILETA2', 'PRFGS_1', 'AHSERL2', 'FASC161ACP', 'GLCS2', 'UAGDP_c', 'P5CCD', 'UGP_PGM_c',
    }
    reactions_to_add = set()
    for reaction_id in reaction_ids:
        reaction = universal_model.reactions.get_by_id(reaction_id)
        reactions_to_add.add(reaction.copy())
    model.add_reactions(reactions_to_add)
    return model


def process_model_steps():
    # Fix compartments
    model = update_1(starting_model)
    model = update_2(model)
    model = update_3(model)
    # don't add these reactions yet
    # model = update_4(model)
    write_model(model)


if __name__ == '__main__':
    process_model_steps()
    model_paths = [s_model_path, output_model_path]
    diff(
        [
            *model_paths,
            '--filename', os.path.join(_file_path, 'model_differences.html')
        ]
    )
