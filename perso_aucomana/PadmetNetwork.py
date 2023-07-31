from padmet.classes.padmetSpec import PadmetSpec
from typing import Set, Dict


def try_key_assignment(dictionary, key):
    try:
        value = dictionary[key]
    except KeyError:
        value = None
    return value


class PadmetNetwork:

    def __init__(self, padmet_network):
        padmet_spec = PadmetSpec(padmet_network)
        self.reactions: Dict[str, Reaction] = dict()
        self.compounds: Dict[str, Compound] = dict()
        self.genes: Dict[str, Gene] = dict()
        self.proteins: Dict[str, Protein] = dict()
        self.classes: Dict[str, Class] = dict()
        self.pathways: Dict[str, Pathway] = dict()
        self.__init_attributes(padmet_spec)

    def __init_attributes(self, p_spec):
        for node in p_spec.dicOfNode.values():
            # print(node.type)
            if node.type == 'reaction':
                self.reactions[node.id] = Reaction(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
            elif node.type == 'compound':
                self.compounds[node.id] = Compound(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
            elif node.type == 'gene':
                self.genes[node.id] = Gene(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
            elif node.type == 'protein':
                self.proteins[node.id] = Protein(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
            elif node.type == 'class':
                self.classes[node.id] = Class(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
            elif node.type == 'pathway':
                self.pathways[node.id] = Pathway(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec)
        for pw_id in self.pathways:
            self.pathways[pw_id].calculate_completion_rate(self.reactions.keys())


class Reaction:
    def __init__(self, rxn_node, rxn_rlt_in, p_spec):
        self.id = rxn_node.id
        self.common_name = try_key_assignment(rxn_node.misc, 'COMMON-NAME')
        self.direction = try_key_assignment(rxn_node.misc, 'DIRECTION')
        self.ec_number = try_key_assignment(rxn_node.misc, 'EC-NUMBER')
        self.name = set()
        self.is_class = set()
        self.xref = None
        self.supp_data = dict()
        self.reconstruction_data = dict()
        self.in_pathways = set()
        self.reactants = set()
        self.products = set()
        self.linked_genes = set()
        self.__init_rlt_in(rxn_rlt_in, p_spec)

    def __init_rlt_in(self, rxn_rlt_in, p_spec):
        if rxn_rlt_in is not None:
            for rlt in rxn_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'consumes':
                    self.reactants.add(rlt.id_out)
                elif rlt.type == 'produces':
                    self.products.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_reconstructionData':
                    self.reconstruction_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_suppData':
                    self.supp_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'is_in_pathway':
                    self.in_pathways.add(rlt.id_out)
                elif rlt.type == 'is_linked_to':
                    self.linked_genes.add(rlt.id_out)
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))


class Compound:

    def __init__(self, cpd_node, cpd_rlt_in, p_spec):
        self.id = cpd_node.id
        self.smiles = try_key_assignment(cpd_node.misc, 'SMILES')
        self.common_name = try_key_assignment(cpd_node.misc, 'COMMON-NAME')
        self.molecular_weight = try_key_assignment(cpd_node.misc, 'MOLECULAR-WEIGHT')
        self.inchi = try_key_assignment(cpd_node.misc, 'INCHI')
        self.inchi_key = try_key_assignment(cpd_node.misc, 'INCHI-KEY')
        self.name = set()
        self.is_class = set()
        self.xref = None
        self.supp_data = dict()
        self.__init_rlt_in(cpd_rlt_in, p_spec)

    def __init_rlt_in(self, cpd_rlt_in, p_spec):
        if cpd_rlt_in is not None:
            for rlt in cpd_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_suppData':
                    self.supp_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))


class Gene:

    def __init__(self, g_node, g_rlt_in, p_spec):
        self.id = g_node.id
        self.transcription_direction = try_key_assignment(g_node.misc, 'TRANSCRIPTION-DIRECTION')
        self.centisome_position = try_key_assignment(g_node.misc, 'CENTISOME-POSITION')
        self.left_end_position = try_key_assignment(g_node.misc, 'LEFT-END-POSITION')
        self.right_end_position = try_key_assignment(g_node.misc, 'RIGHT-END-POSITION')
        self.is_class = set()
        self.name = set()
        self.xref = None
        self.supp_data = dict()
        self.proteins_coding = set()
        self.__init_rlt_in(g_rlt_in, p_spec)

    def __init_rlt_in(self, g_rlt_in, p_spec):
        if g_rlt_in is not None:
            for rlt in g_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_suppData':
                    self.supp_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))
                elif rlt.type == 'codes_for':
                    self.proteins_coding.add(rlt.id_out)


class Protein:

    def __init__(self, prot_node, prot_rlt_in, p_spec):
        self.id = prot_node.id
        self.is_class = set()
        self.name = set()
        self.xref = None
        self.supp_data = dict()
        self.rxn_catalysing = None
        self.linked_species = None
        self.__init_rlt_in(prot_rlt_in, p_spec)

    def __init_rlt_in(self, prot_rlt_in, p_spec):
        if prot_rlt_in is not None:
            for rlt in prot_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_suppData':
                    self.supp_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))


class Pathway:

    def __init__(self, pw_node, pw_rlt_in, p_spec):
        self.id = pw_node.id
        self.common_name = try_key_assignment(pw_node.misc, 'COMMON-NAME')
        self.taxonomic_range = try_key_assignment(pw_node.misc, 'TAXONOMIC-RANGE')
        self.input_cpd = try_key_assignment(pw_node.misc, 'INPUT-COMPOUNDS')
        self.output_cpd = try_key_assignment(pw_node.misc, 'OUTPUT-COMPOUNDS')
        self.rxn_order = try_key_assignment(pw_node.misc, 'REACTIONS-ORDER')
        self.nb_rxn = 0
        if self.rxn_order is not None:
            self.rxn_order = self.rxn_order[0].split(',')
            self.nb_rxn = len(set(self.rxn_order))
        self.is_class = set()
        self.name = set()
        self.xref = None
        self.in_pathways = None
        self.rxn_present = set()
        self.completion_rate = None
        self.completion_str = None
        self.__init_rlt_in(pw_rlt_in, p_spec)

    def __init_rlt_in(self, pw_rlt_in, p_spec):
        if pw_rlt_in is not None:
            for rlt in pw_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))

    def calculate_completion_rate(self, rxn_set):
        if self.rxn_order is not None:
            for r in self.rxn_order:
                if r in rxn_set:
                    self.rxn_present.add(r)
        if self.nb_rxn != 0:
            self.completion_rate = len(self.rxn_present) / self.nb_rxn
        else:
            self.completion_rate = float(1)
        self.completion_str = f'{len(self.rxn_present)}/{self.nb_rxn}'


class Class:

    def __init__(self, cls_node, cls_rlt_in, p_spec):
        self.id = cls_node.id
        self.common_name = try_key_assignment(cls_node.misc, 'COMMON-NAME')
        self.is_class = set()
        self.name = set()
        self.xref = None
        self.supp_data = dict()
        self.__init_rlt_in(cls_rlt_in, p_spec)

    def __init_rlt_in(self, cls_rlt_in, p_spec):
        if cls_rlt_in is not None:
            for rlt in cls_rlt_in:
                if rlt.type == 'is_a_class':
                    self.is_class.add(rlt.id_out)
                elif rlt.type == 'has_xref':
                    self.xref = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_suppData':
                    self.supp_data[rlt.id_out] = p_spec.dicOfNode[rlt.id_out].misc
                elif rlt.type == 'has_name':
                    self.name = self.name.union(set(p_spec.dicOfNode[rlt.id_out].misc['LABEL']))
