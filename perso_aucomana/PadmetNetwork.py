from padmet.classes.padmetSpec import PadmetSpec


def try_key_assignment(dictionary, key):
    try:
        value = dictionary[key]
    except KeyError:
        value = None
    return value


class PadmetNetwork:

    def __init__(self, padmet_spec):
        self.reactions = self.__init_reactions(padmet_spec)
        self.compounds = None
        self.genes = None
        self.proteins = None
        self.classes = None
        self.pathways = None

    @staticmethod
    def __init_reactions(p_spec):
        reaction_set = set()
        for node in p_spec.dicOfNode.values():
            if node.type == 'reaction':
                reaction_set.add(Reaction(node))

        return reaction_set


class Class:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None


class Reaction:
    def __init__(self, rxn_node):
        self.id = rxn_node.id
        self.is_class = None
        self.name = try_key_assignment(rxn_node.misc, 'COMMON-NAME')
        self.xref = None
        self.supp_data = None
        self.reconstruction_data = None
        self.in_pathways = None
        self.reactants = None
        self.products = None
        self.linked_genes = None
        self.direction = try_key_assignment(rxn_node.misc, 'DIRECTION')
        self.ec_number = try_key_assignment(rxn_node.misc, 'EC-NUMBER')


class Compound:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None


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


PADMET_NW_FILE = '/home/phamongi/Documents/Runs/run62/networks/PADMETs/Ascophyllum-nodosum_MALE.padmet'
P = PadmetSpec(PADMET_NW_FILE)
N = PadmetNetwork(P)

print([r.id for r in N.reactions])
