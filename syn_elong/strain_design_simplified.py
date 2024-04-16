#!/usr/bin/env python3
#
# Copyright 2022 Max Planck Insitute Magdeburg
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
#
#
"""Function: computing metabolic strain designs (compute_strain_designs)"""

from contextlib import redirect_stdout, redirect_stderr
import numpy as np
import json
import io
from copy import deepcopy
import logging
from cobra import Model
from cobra.manipulation import rename_genes
from straindesign import SDModule, SDSolutions, select_solver, fva, DisableLogger, SDMILP
from straindesign.names import *
from straindesign.networktools import   remove_ext_mets, remove_dummy_bounds, bound_blocked_or_irrevers_fva, \
                                        remove_irrelevant_genes, extend_model_gpr, extend_model_regulatory, \
                                        compress_model, compress_modules, compress_ki_ko_cost, expand_sd, filter_sd_maxcost

class StrainDesign(object):
    def __init__(self, model: Model, sd_modules: [SDModule], solver: str = None, M: int = None, compress: bool = True,
                 gene_kos: bool = False, ko_cost: dict = None, ki_cost: dict = None, gko_cost: dict = None, gki_cost: dict = None,
                 reg_cost: dict = None, solution_approach: str = 'any', advanced: bool = False, use_scenario: bool = False, ) -> object:
        """
        Computes strain designs for a user-defined strain design problem

        A number of arguments can be specified to detail the problem and influence the solution process.
        This function supports the computation of Minimal Cut Sets (MCS), OptKock, RobustKnock and OptCouple
        strain designs. It is possible to combine any of the latter ones with the MCS approach, e.g., to
        engineer growth coupled production, but also suppress the production of an undesired by-product.
        The computation can be started in two different ways. Either by specifying the computation parameters
        individually or reuse a parameters dictionary from a previous computation. CNApy stores strain design
        setup dics as JSON ".sd"-files that can be loaded in python and used as an input for this function.

        Example:
            sols = compute_strain_designs(model, sd_modules=[sd_module1, sd_module2], solution_approach = 'any')

        Args:
            model (cobra.Model):
                A metabolic model that is an instance of the cobra.Model class. The model may or may not
                contain genes/GPR-rules.

            sd_modules ([straindesign.SDModule]):
                List of strain design modules that describe the sub-problems, such as the MCS-like protection
                or suppression of flux subspaces or the OptKnock, RobustKnock or OptCouple objective and
                constraints. The list of modules determines the global objective function of the strain design
                computation. If only SUPPRESS and PROTECT modules are used, the strain design computation is
                MCS-like, such that the number of interventions is minimized. If a module for one of the nested
                optimization approaches is used, the global objective function is retrieved from this module.
                The number of SUPPRESS and PROTECT modules is unrestricted and can be combined with the other
                modules, however only one of the modules OPTKNOCK, ROBUSKNOCK and OPTCOUPLE may be used at a time.
                For details, see SDModule.

            solver (optional (str)): (Default: same as defined in model / COBRApy)
                The solver that should be used for preparing and carrying out the strain design computation.
                Allowed values are 'cplex', 'gurobi', 'scip' and 'glpk'.

            M (optional (int)): (Default: None)
                If this value is specified (and non-zero, not None), the computation uses the big-M
                method instead of indicator constraints. Since GLPK does not support indicator constraints it uses
                the big-M method by default (with M=1000). M should be chosen 'sufficiently large' to avoid computational
                artifacts and 'sufficiently small' to avoid numerical issues.

            compress (optional (bool)): (Default: True)
                If 'True', the iterative network compressor is used.

            gene_kos (optional (bool)): (Default: False)
                If 'True', strain designs are computed based on gene-knockouts instead of reaction knockouts. This
                parameter needs not be defined if any of ki_cost, ko_cost, gki_cost, gko_cost and reg_cost is used.
                By default, reactions are considered as knockout targets.

            ko_cost (optional (dict)): (Default: None)
                A dictionary of reaction identifiers and their associated knockout costs. If not specified, all reactions
                are treated as knockout candidates, equivalent to ko_cost = {'r1':1, 'r2':1, ...}. If a subset of reactions
                is listed in the dict, all others are not considered as knockout candidates.

            ki_cost (optional (dict)): (Default: None)
                A dictionary of reaction identifiers and their associated costs for addition. If not specified, all reactions
                are treated as knockout candidates. Reaction addition candidates must be present in the original model with
                the intended flux boundaries **after** insertion. Additions are treated adversely to knockouts, meaning that
                their exclusion from the network is not associated with any cost while their presence entails intervention costs.

            gko_cost (optional (dict)): (Default: None)
                A dictionary of gene identifiers and their associated knockout costs. To reference genes, gene IDs can be
                used,as well as gene names. If not specified, genes are not treated as knockout candidates. An exception is
                the 'gene_kos' argument. If 'gene_kos' is used, all genes are treated as knockout candidates with intervention
                costs of 1. This is equivalent to gko_cost = {'g1':1, 'g2':1, ...}.

            gki_cost (optional (dict)): (Default: None)
                A dictionary of gene identifiers and their associated addition costs. To reference genes, gene IDs can be
                used, as well as gene names. If not specified, none of the genes are treated as addition candidates.

            reg_cost (optional [dict]): ( Default: None)
                Regulatory interventions candidates can be optionally specified as a list. Thereby, the constraint marking the
                regulatory intervention is put as key and the associated intervention cost is used as the corresponding value.
                E.g., reg_cost = {'1 EX_o2_e = -1': 1, ... <other regulatory interventions>}. Instead of strings, constraints
                can also be passed as lists. reg_cost = {[{'EX_o2_e':1}, '=', -1]: 1, ...}

            solution_approach (optional (str)): ( Default: 'best')
                The approach used to find strain designs. Possible values are 'any', 'best' or 'populate'. 'any' is usually the
                fastest option, since optimality is not enforced. Hereby computed MCS are still irreducible intervention sets,
                however, not MCS with the fewest possible number of interventions. 'best' computes globally optimal strain designs,
                that is, MCS with the fewest number of interventions, OptKnock strain designs with the highest possible production
                rate, OptCouple strain designs with the hightest growth coupling potential etc.. 'populate' does the same as 'best',
                but makes use of CPLEX' and Gurobi's populate function to generate multiple strain designs. It is identical to 'best'
                when used with SCIP or GLPK.
                Attention:
                If 'any' used with OptKnock, for instance, the MILP may return the wild type as a possible immediately. Technically,
                the wildtype fulfills the criterion of maximal growth (inner objective) and maximality of the global objective is
                omitted by using 'any', so that carrying no product synthesis is permitted. Additional constraints can be used
                in the OptKnock problem to circumvent this. However, Optknock should generally be used with the 'best' option.
            """

        self.kwargs_milp = None
        self.reac_map = None
        self.compressed_model = None
        self.uncompressed_model = None
        self.orig_model = None
        self.orig_gko_cost = None
        self.orig_gki_cost = None
        logging.info('Preparing strain design computation.')
        self.model = model
        if isinstance(sd_modules, SDModule):
            sd_modules = [sd_modules]
        self.orig_sd_modules = [m.copy() for m in sd_modules]
        self.sd_modules = [m.copy() for m in sd_modules]
        self.solver = select_solver(solver, self.model)
        logging.info(f'  Using {self.solver} for solving LPs during preprocessing.')
        self.big_M = M
        self.compress = compress
        self.gene_kos = gene_kos
        self.ko_cost = ko_cost
        self.ki_cost = ki_cost
        self.gko_cost = gko_cost
        self.gki_cost = gki_cost
        self.reg_cost = reg_cost

        self.solution_approach = solution_approach
        self.advanced = advanced
        self.use_scenario = use_scenario
        self.has_gene_names = False
        # uncompressed costs and model
        self.uncompressed_ko_cost = {}
        self.uncompressed_ki_cost = {}
        self.uncompressed_gko_cost = {}
        self.uncompressed_gki_cost = {}
        self.uncompressed_reg_cost = {}
        # compressed costs and model
        self.compressed_map_reac = None
        self.compressed_ki_cost = None
        self.compressed_ko_cost = None
        self.essential_kis = set()
        self.essential_reactions = set()
        self.orig_ko_cost = {}
        self.orig_ki_cost = {}
        self.orig_reg_cost = {}

        self.orig_gko_cost = {}
        self.orig_gki_cost = {}
        self.steps()
 

    def check_args(self):
        if self.big_M is None:
            if self.solver in [CPLEX, GUROBI]:
                self.big_M = 1000
            else:
                self.big_M = np.inf
        if self.gene_kos or self.gko_cost is not None or self.gki_cost is not None:
            self.check_gene_names()
            if self.gko_cost is None or not self.gko_cost:
                if self.has_gene_names:
                    self.gko_cost = {k: 1.0 for k in self.model.genes.list_attr('name')}
                else:
                    self.gko_cost = {k: 1.0 for k in self.model.genes.list_attr('id')}

        # work out reaction knockouts/ins
        if self.ko_cost is None:
            if not self.gene_kos:
                self.ko_cost = {k: 1.0 for k in self.model.reactions.list_attr('id')}
            else:
                self.ko_cost = {}
        if self.ki_cost is None:
            self.ki_cost = {}
        if self.reg_cost is None:
            self.reg_cost = {}
        self.uncompressed_ko_cost = deepcopy(self.ko_cost)
        self.uncompressed_ki_cost = deepcopy(self.ki_cost)
        self.uncompressed_reg_cost = deepcopy(self.reg_cost)
        self.orig_ko_cost = deepcopy(self.ko_cost)
        self.orig_ki_cost = deepcopy(self.ki_cost)
        self.orig_reg_cost = deepcopy(self.reg_cost)
        if self.solution_approach not in [ANY, BEST, POPULATE]:
            raise Exception("Solution approach must be one of 'any', 'best' or 'populate'.")
            # suppress standard output from copying model
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()), DisableLogger():
            self.orig_model = self.model.copy()
            self.uncompressed_model = self.model.copy()
            self.compressed_model = self.model.copy()

    def steps(self):
        self.check_args()

        if self.gene_kos:
            self.setup_gko()
        self.preprocess_model()
        if self.gene_kos:
            self.modify_for_genes()
            # move gene ko costs to reaction ko costs
            self.uncompressed_ko_cost.update(self.uncompressed_gko_cost)
            self.uncompressed_ki_cost.update(self.uncompressed_gki_cost)
        self.uncompressed_ko_cost.update(extend_model_regulatory(self.compressed_model, self.uncompressed_reg_cost))

        # make copy to modify for compressed
        self.compressed_ko_cost = self.uncompressed_ko_cost
        self.compressed_ki_cost = self.uncompressed_ki_cost

        if self.compress:
            self.compress_model()

        self.setup_milp_args()

    def check_gene_names(self):
        if not hasattr(self.model, 'genes') and self.model.genes:
            self.gene_kos = False
            logging.warning("Gene knockouts were specified, but no genes are defined in the model. ")
            return
        # genes must not begin with number, put a 'g' in front of genes that start with a number
        if any([True for g in self.model.genes if g.id[0].isdigit()]):
            logging.warning("Gene IDs must not start with a digit. Inserting prefix 'g' where necessary.")
            rename_genes(self.model, {g.id: 'g' + g.id for g in self.model.genes if g.id[0].isdigit()})
        used_gene_ids = set()
        if self.gko_cost is not None:
            used_gene_ids.update(self.gko_cost.keys())
        if self.gki_cost is not None:
            used_gene_ids.update(self.gki_cost.keys())
        if self.reg_cost is not None:
            used_gene_ids.update(self.reg_cost.keys())
        # check if genes are in the model
        m_genes = set(g.name for g in self.model.genes)
        overlap = used_gene_ids.intersection(m_genes)
        if np.all([len(g.name) for g in self.model.genes]) and len(overlap):
            self.has_gene_names = True

    def setup_gko(self):
        if self.gki_cost is None:
            self.gki_cost = {}
        if self.gko_cost is None:
            self.gko_cost = {}
        self.uncompressed_gko_cost = deepcopy(self.gko_cost)
        self.uncompressed_gki_cost = deepcopy(self.gki_cost)
        self.orig_gko_cost = deepcopy(self.uncompressed_gko_cost)
        self.orig_gki_cost = deepcopy(self.uncompressed_gki_cost)
        gene_interventions = set(self.uncompressed_gko_cost.keys()).union(self.uncompressed_gki_cost.keys())
        reactions_interventions = set(self.uncompressed_ko_cost.keys()).union(self.uncompressed_ki_cost.keys())
        overlap_1 = False
        for r in reactions_interventions:
            for g in self.uncompressed_model.reactions.get_by_id(r).genes:
                if g in gene_interventions:
                    overlap_1 = True

        overlap_2 = set(self.uncompressed_gko_cost.keys()).intersection(set(self.uncompressed_gki_cost.keys()))
        overlap_3 = set(self.uncompressed_ko_cost.keys()).intersection(set(self.uncompressed_ki_cost.keys()))

        if overlap_1 or len(overlap_2) or len(overlap_3):
            raise Exception('Specified gene and reaction knock-out/-in costs contain overlap. '
                            'Make sure that metabolic interventions are enabled either through reaction or '
                            'through gene interventions and are defined either as knock-ins or as knock-outs.')

    def preprocess_model(self):
        # remove external metabolites
        remove_ext_mets(self.compressed_model)
        # replace model bounds with +/- inf if above a certain threshold
        remove_dummy_bounds(self.model)
        # FVAs to identify blocked, irreversible and essential reactions, as well as non-bounding bounds
        logging.info('  FVA to identify blocked reactions and irreversibility.')
        bound_blocked_or_irrevers_fva(self.model, solver=self.solver)
        logging.info('  FVA(s) to identify essential reactions.')

        for m in self.sd_modules:
            if m[MODULE_TYPE] != SUPPRESS:  # Essential reactions can only be determined from desired
                # or opt-/robustknock modules
                flux_limits = fva(self.compressed_model, solver=self.solver, constraints=m[CONSTRAINTS])
                for (reac_id, limits) in flux_limits.iterrows():
                    if np.min(abs(limits)) > 1e-10 and np.prod(np.sign(limits)) > 0:  # find essential
                        self.essential_reactions.add(reac_id)
        logging.info(f'Essential reactions: {len(self.essential_reactions)}')
        # remove ko-costs (and thus knock ability) of essential reactions
        for er in self.essential_reactions:
            if er in self.uncompressed_ko_cost:
                self.uncompressed_ko_cost.pop(er)


    def modify_for_genes(self):
        # If computation of gene-based intervention strategies
        if self.compress:
            num_genes = len(self.compressed_model.genes)
            num_gpr = len([True for r in self.model.reactions if r.gene_reaction_rule])
            logging.info(f'Preprocessing GPR rules {num_genes} genes, {num_gpr} gpr rules).')
            # removing irrelevant genes will also remove essential reactions from the list of knockable genes
            self.uncompressed_gko_cost = remove_irrelevant_genes(
                self.compressed_model,
                list(self.essential_reactions),
                self.uncompressed_gki_cost,
                self.uncompressed_gko_cost
            )
            c_num_genes = len(self.compressed_model.genes)
            c_num_gpr = len([True for r in self.compressed_model.reactions if r.gene_reaction_rule])
            if  c_num_genes < num_genes or c_num_gpr < num_gpr:
                num_genes = len(self.compressed_model.genes)
                num_gpr = len([True for r in self.compressed_model.reactions if r.gene_reaction_rule])
                logging.info(f'  Simplified to {num_genes} genes and {num_gpr} gpr rules.')
        logging.info('  Extending metabolic network with gpr associations.')
        reac_map = extend_model_gpr(self.compressed_model, self.has_gene_names)
        for i, m in enumerate(self.sd_modules):
            for p in [CONSTRAINTS, INNER_OBJECTIVE, OUTER_OBJECTIVE, PROD_ID]:
                if p in m and m[p] is not None:
                    if p == CONSTRAINTS:
                        for c in m[p]:
                            for k in list(c[0].keys()):
                                v = c[0].pop(k)
                                for n, w in reac_map[k].items():
                                    c[0][n] = v * w
                    if p in [INNER_OBJECTIVE, OUTER_OBJECTIVE, PROD_ID]:
                        for k in list(m[p].keys()):
                            v = m[p].pop(k)
                            for n, w in reac_map[k].items():
                                m[p][n] = v * w



    def compress_model(self):
        logging.info(f'Compressing Network ({len(self.uncompressed_model.reactions)} reactions).')
        # compress network by lumping sequential and parallel reactions alternatively.
        # Exclude reactions named in strain design modules from parallel compression
        no_par_compress_reacs = set()
        for m in self.sd_modules:
            for p in [CONSTRAINTS, INNER_OBJECTIVE, OUTER_OBJECTIVE, PROD_ID]:
                if p in m and m[p] is not None:
                    param = m[p]
                    if p == CONSTRAINTS:
                        for c in param:
                            for k in c[0].keys():
                                no_par_compress_reacs.add(k)
                    if p in [INNER_OBJECTIVE, OUTER_OBJECTIVE, PROD_ID]:
                        for k in param.keys():
                            no_par_compress_reacs.add(k)
        cmp_map_reac = compress_model(self.compressed_model, no_par_compress_reacs)
        # compress information in strain design modules
        self.sd_modules = compress_modules(self.sd_modules, cmp_map_reac)
        # compress ko_cost and ki_cost
        self.compressed_ko_cost, self.compressed_ki_cost, self.compressed_map_reac = compress_ki_ko_cost(
            self.compressed_ko_cost, self.compressed_ki_cost, cmp_map_reac
        )

        # An FVA to identify essentials before building and launching MILP (not sure if this has an effect)
        logging.info('  FVA(s) in compressed model to identify essential reactions.')
        essential_reacs = set()
        for m in self.sd_modules:
            if m[MODULE_TYPE] != SUPPRESS:  # Essential reactions can only be determined from desired
                # or opt-/robustknock modules
                flux_limits = fva(self.compressed_model, solver=self.solver, constraints=m[CONSTRAINTS])
                for (reac_id, limits) in flux_limits.iterrows():
                    if np.min(abs(limits)) > 1e-10 and np.prod(np.sign(limits)) > 0:  # find essential
                        essential_reacs.add(reac_id)

        # remove ko-costs (and thus knockability) of essential reactions
        for er in essential_reacs:
            if er in self.compressed_ko_cost:
                self.compressed_ko_cost.pop(er)

        self.essential_kis = set(self.compressed_ki_cost[er] for er in essential_reacs if er in self.compressed_ki_cost)


    def setup_milp_args(self):
        # Build MILP
        kwargs_milp = {
            SOLVER: self.solver,
            'M' : self.big_M,
            KOCOST:self.compressed_ko_cost,
            KICOST:self.compressed_ki_cost,
            'essential_kis':self.essential_kis,
            }
        self.kwargs_milp = kwargs_milp
        logging.info("Finished preprocessing:")
        logging.info(f"  Model size: {len(self.compressed_model.reactions)} rxns, "
                     f"{len(self.compressed_model.metabolites)} metabolites")
        print(len(self.compressed_ko_cost), len(self.compressed_ki_cost), len(self.essential_kis))
        logging.info(f"  {len(self.compressed_ko_cost) + len(self.compressed_ki_cost) - len(self.essential_kis)} "
                     f"targetable reactions")

    def run(self, max_solutions=5, time_limit=60, max_cost=np.inf, solution_approach=ANY):
        """
        max_cost (optional (int)): (Default: inf):
            The maximum cost threshold for interventions. Every possible intervention is associated with a
            cost value (1, by default). Strain designs cannot exceed the max_cost threshold. Individual
            intervention cost factors may be defined through ki_cost, ko_cost, gki_cost, gko_cost and reg_cost.

        max_solutions (optional (int)): (Default: inf)
            The maximum number of MILP solutions that are generated for a strain design problem. The number of returned
            strain designs is usually larger than the number of max_solutions, since a MILP solution is decompressed
            to multiple strain designs. When the compress-flag is set to 'False' the number of returned solutions is
            equal to max_solutions.

        time_limit (optional (int)): (Default: inf)
            The time limit in seconds for the MILP-solver.


        """
        self.kwargs_milp['max_cost'] = max_cost
        sd_milp = SDMILP(self.compressed_model, self.sd_modules, **self.kwargs_milp)

        kwargs_computation = {
            'show_no_ki': True,
            'max_solutions': max_solutions,
            'time_limit': time_limit,
        }

        # solve MILP
        if solution_approach == ANY:
            cmp_sd_solution = sd_milp.compute(**kwargs_computation)
        elif solution_approach == BEST:
            cmp_sd_solution = sd_milp.compute_optimal(**kwargs_computation)
        elif solution_approach == POPULATE:
            cmp_sd_solution = sd_milp.enumerate(**kwargs_computation)
        else:
            raise ERROR

        logging.info('  Decompressing.')
        if cmp_sd_solution.status in [OPTIMAL, TIME_LIMIT_W_SOL]:
            sd = expand_sd(cmp_sd_solution.get_reaction_sd_mark_no_ki(), self.compressed_map_reac)
            sd = filter_sd_maxcost(sd, max_cost, self.uncompressed_ko_cost, self.uncompressed_ki_cost)
            sd = postprocess_reg_sd(self.uncompressed_reg_cost, sd)
        else:
            sd = []

        setup = deepcopy(cmp_sd_solution.sd_setup)
        setup.update({MODULES: self.orig_sd_modules, KOCOST: self.orig_ko_cost,
                      KICOST: self.orig_ki_cost, REGCOST: self.orig_reg_cost})
        if self.gene_kos:
            setup.update({GKOCOST: self.orig_gko_cost, GKICOST: self.orig_gki_cost})
        sd_solutions = SDSolutions(self.orig_model, sd, cmp_sd_solution.status, setup)
        logging.info(str(len(sd)) + ' solutions found.')

        return sd_solutions





def postprocess_reg_sd(reg_cost, sd):
    """Postprocess regulatory interventions

    Mark regulatory interventions with true or false"""
    for s in sd:
        for k, v in reg_cost.items():
            if k in s:
                s.pop(k)
                s.update({v['str']: True})
            else:
                s.update({v['str']: False})
    return sd