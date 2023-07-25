
class PadmetNetwork:

    def __init__(self, padmet_network_file):
        self.reactions = None
        self.compounds = None
        self.genes = None
        self.proteins = None
        self.classes = None
        self.pathways = None


class Class:

    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None

class Reaction:
    def __init__(self):
        self.id = None
        self.is_class = None
        self.name = None
        self.xref = None
        self.supp_data = None
        self.reconstruction_data = None
        self.in_pathways = None
        self.reactants = None
        self.products = None
        self.linked_genes = None

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