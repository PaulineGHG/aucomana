from padmet.classes.padmetSpec import PadmetSpec
import time
start = time.time()


def try_key_assignment(dictionary, key):
    try:
        value = dictionary[key]
    except KeyError:
        value = None
    return value


class PadmetNetwork:

    def __init__(self, padmet_spec):
        self.reactions = set()
        self.compounds = set()
        self.genes = None
        self.proteins = None
        self.classes = None
        self.pathways = None
        self.__init_attributes(padmet_spec)

    def __init_attributes(self, p_spec):
        for node in p_spec.dicOfNode.values():
            if node.type == 'reaction':
                self.reactions.add(Reaction(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec))
            elif node.type == 'compound':
                self.compounds.add(Compound(node, try_key_assignment(p_spec.dicOfRelationIn, node.id), p_spec))


class Reaction:
    def __init__(self, rxn_node, rxn_rlt_in, p_spec):
        self.id = rxn_node.id
        self.name = try_key_assignment(rxn_node.misc, 'COMMON-NAME')
        self.direction = try_key_assignment(rxn_node.misc, 'DIRECTION')
        self.ec_number = try_key_assignment(rxn_node.misc, 'EC-NUMBER')
        self.is_class = set()
        self.xref = None
        self.supp_data = set()
        self.reconstruction_data = set()
        self.in_pathways = set()
        self.reactants = set()
        self.products = set()
        self.linked_genes = set()
        self.__init_rlt_in(rxn_rlt_in, p_spec)

    def __init_rlt_in(self, rxn_rlt_in, p_spec):
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
                self.reconstruction_data.add(rlt.id_out)
            elif rlt.type == 'has_suppData':
                self.supp_data.add(rlt.id_out)
            elif rlt.type == 'is_in_pathway':
                self.in_pathways.add(rlt.id_out)
            elif rlt.type == 'is_linked_to':
                self.linked_genes.add(rlt.id_out)


class Compound:

    def __init__(self, cpd_node, rxn_rlt_in, p_spec):
        self.id = cpd_node.id
        self.smiles = try_key_assignment(cpd_node.misc, 'SMILES')
        self.name = try_key_assignment(cpd_node.misc, 'COMMON-NAME')
        self.molecular_weight = try_key_assignment(cpd_node.misc, 'MOLECULAR-WEIGHT')
        self.inchi = try_key_assignment(cpd_node.misc, 'INCHI')
        self.inchi_key = try_key_assignment(cpd_node.misc, 'INCHI-KEY')
        self.is_class = set()
        self.xref = None
        self.supp_data = set()
        self.__init_rlt_in(rxn_rlt_in, p_spec)

    def __init_rlt_in(self, rxn_rlt_in, p_spec):
        for rlt in rxn_rlt_in:
            if rlt.type == 'is_a_class':
                self.is_class.add(rlt.id_out)
            elif rlt.type == 'has_xref':
                self.xref = p_spec.dicOfNode[rlt.id_out].misc
            elif rlt.type == 'has_suppData':
                self.supp_data.add(p_spec.dicOfNode[rlt.id_out].misc)


class Protein:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None
        self.rxn_catalysing = None
        self.linked_species = None


class Gene:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None
        self.proteins_coding = None


class Pathway:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.in_pathways = None
        self.completion_rate = None


class Class:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None


PADMET_NW_FILE = '/home/phamongi/Documents/Runs/run62/networks/PADMETs/Ascophyllum-nodosum_MALE.padmet'
P = PadmetSpec(PADMET_NW_FILE)
N = PadmetNetwork(P)

print([c.supp_data for c in N.compounds])

end = time.time()
print(end-start)
