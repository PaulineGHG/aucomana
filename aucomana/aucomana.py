"""
AuCoMAna class
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as m_colors
from typing import Iterable, List
from aucomana.utils.utils import get_grp_set
from aucomana.utils.pathways import Pathways
from aucomana.utils.reactions import Reactions
from supervenn import supervenn
from Bio import AlignIO
from Bio import SeqIO
from io import StringIO
from Bio.Align.Applications import MuscleCommandline


class GroupsAnalysis:
    """
    Attributes
    ----------
    group_template: str
    reactions_tsv: str
    pathways_tsv: str
    genes_tsv: str
    metabolites_tsv: str
    """

    STAT = (" mean", " med", " sd", " min", " max")
    TO_CALCULATE = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")

    def __init__(self, compare_path: str, group_template: str):
        """ Init the AuCoMAna class.

        Parameters
        ----------
        compare_path: str
            Path of the directory with compare padmets results
        group_template: str
            Groups template tsv file
        """
        self.group_template = group_template
        self.reactions_tsv = os.path.join(compare_path, 'reactions.tsv')
        self.pathways_tsv = os.path.join(compare_path, 'pathways.tsv')
        self.genes_tsv = os.path.join(compare_path, 'genes.tsv')
        self.metabolites_tsv = os.path.join(compare_path, 'metabolites.tsv')

    def group_reactions_comparison(self, groups_comp: Iterable[str], group_analysis: str = 'all',
                                   output='groups_comp_rxn.png'):
        species_analysis = get_grp_set(self.group_template, group_analysis)

        reaction = Reactions(file_reactions_tsv=self.reactions_tsv,
                             species_list=list(species_analysis))
        dic_groups = dict()
        for group in list(groups_comp) + [group_analysis]:
            dic_groups[group] = get_grp_set(self.group_template, (group, group_analysis))

        absent = reaction.get_rxn_absent(species=species_analysis, unique=True)
        abs_groups = dict()

        for group in list(groups_comp) + [group_analysis]:
            abs_groups[group] = list()

            for sp in dic_groups[group]:
                abs_groups[group].append(absent[sp][0])

        ind = ["absent"]
        col = list(abs_groups.keys())
        df = pd.DataFrame(columns=col, index=ind)
        for group, numbers_list in abs_groups.items():
            df.loc[ind[0], group] = round(float(np.mean(numbers_list)), 2)

        ax = df.plot.bar()
        for container in ax.containers:
            ax.bar_label(container)
        ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        ax.set_ylabel("Number of pathways")
        ax.set_title(f'Unique absent reactions for {"/".join(groups_comp)} groups')
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(output)

    def group_pathway_completion_comparison(self, groups_comp: Iterable[str], group_analysis: str = 'all',
                                            output='groups_comp_pw.png'):
        species_analysis = get_grp_set(self.group_template, group_analysis)

        pathway = Pathways(file_pathways_tsv=self.pathways_tsv,
                           species_list=list(species_analysis))
        dic_groups = dict()
        for group in list(groups_comp) + [group_analysis]:
            dic_groups[group] = get_grp_set(self.group_template, (group, group_analysis))

        absent = pathway.get_pw_absent(species=species_analysis, unique=True)
        mini = pathway.get_pw_min(species=species_analysis, unique=True)
        incomplete = pathway.get_pw_incomplete(species=species_analysis, unique=True)

        abs_groups = dict()
        min_groups = dict()
        inc_groups = dict()

        for group in list(groups_comp) + [group_analysis]:
            abs_groups[group] = list()
            min_groups[group] = list()
            inc_groups[group] = list()

            for sp in dic_groups[group]:
                abs_groups[group].append(absent[sp][0])
                min_groups[group].append(mini[sp][0])
                inc_groups[group].append(incomplete[sp][0])

        ind = ["absent", "minimal", "incomplete"]
        col = list(abs_groups.keys())
        df = pd.DataFrame(columns=col, index=ind)
        for group, numbers_list in abs_groups.items():
            df.loc[ind[0], group] = round(float(np.mean(numbers_list)), 2)
        for group, numbers_list in min_groups.items():
            df.loc[ind[1], group] = round(float(np.mean(numbers_list)), 2)
        for group, numbers_list in inc_groups.items():
            df.loc[ind[2], group] = round(float(np.mean(numbers_list)), 2)

        ax = df.plot.bar()
        for container in ax.containers:
            ax.bar_label(container)
        ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        ax.set_ylabel("Number of pathways")
        ax.set_title(f'Unique absent/minimal/incomplete pathway for {"/".join(groups_comp)} groups')
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(output)

    def group_supervenn_rxn(self, groups_comp: List[str], fig_size=(16, 8), colors=None, output='supervenn.png'):
        reactions = Reactions(self.reactions_tsv)

        order_sp = []
        dic_sp_groups = dict()
        for group in groups_comp:
            species = get_grp_set(self.group_template, group)
            for sp in species:
                order_sp.append(sp)
                dic_sp_groups[sp] = group

        if colors is None:
            colors = list(m_colors.TABLEAU_COLORS.keys())
        colors_cycle = (len(groups_comp) // len(colors) + 1) * colors

        dic_groups_color = dict()
        for i in range(len(groups_comp)):
            dic_groups_color[groups_comp[i]] = colors_cycle[i]

        rxn_sets = reactions.get_rxn_present()

        plt.figure(figsize=fig_size)
        supervenn([rxn_sets[sp][1] for sp in order_sp], order_sp, side_plots=False,
                  color_cycle=[dic_groups_color[dic_sp_groups[sp]] for sp in order_sp])
        plt.savefig(output)


class SequencesAnalysis:
    def __init__(self, species_sequences_dir: str, compare_path: str):
        self.sp_genome = dict()
        self.sp_proteome = dict()
        self.__init_dicts(species_sequences_dir)
        self.reactions_tsv = os.path.join(compare_path, 'reactions.tsv')
        self.pathways_tsv = os.path.join(compare_path, 'pathways.tsv')
        self.genes_tsv = os.path.join(compare_path, 'genes.tsv')
        self.metabolites_tsv = os.path.join(compare_path, 'metabolites.tsv')

    def __init_dicts(self, species_sequences_dir):
        for species in os.listdir(species_sequences_dir):
            for file in os.listdir(os.path.join(species_sequences_dir, species)):
                if file.endswith('.faa'):
                    self.sp_proteome[species] = os.path.join(species_sequences_dir, species, file)
                elif file.endswith('.fna'):
                    self.sp_genome[species] = os.path.join(species_sequences_dir, species, file)

    def multiple_alignments(self, rxn, source_species: str, output: str, species_list: List[str]=None):
        """ Need "muscle" installed.

        Parameters
        ----------
        rxn
        source_species
        output
        species_list

        Returns
        -------

        """
        # Get reaction instance
        r = Reactions(file_reactions_tsv=self.reactions_tsv)
        # Obtain species list to compare ro th source_species
        if species_list is None:
            species_list = r.species_list
        # Obtain gene associated to rxn and get source sequences from the source species
        genes_assoc = r.get_genes_assoc(rxn)[rxn]
        source_genes =  genes_assoc[source_species]
        seq_source = [SeqIO.to_dict(SeqIO.parse(self.sp_proteome[source_species], 'fasta'))[g] for g in source_genes]

        # Create output directories
        seq_output = os.path.join(output, 'sequences')
        align_output = os.path.join(output, 'alignments')
        self.__create_directories((seq_output, align_output))

        # Align source sequences with sequences of each species :
        all_align_file = os.path.join(align_output, 'all_align.txt')
        with open(all_align_file, 'a') as f:
            f.write(f'All alignments with {source_species} genes :\nseq\tid\tlen(seq)\n')
        all_seq_len = {seq.id: len(seq.seq) for seq in seq_source}
        all_seq = list()

        for sp in species_list:
            if sp != source_species:
                sp_genes = self.__filter_genes(genes_assoc[sp])

                print(f'align {source_genes} with {sp_genes}')
                seq_sp = SeqIO.to_dict(SeqIO.parse(self.sp_proteome[sp], 'fasta'))
                print(seq_sp.keys())
                seq_to_align = [seq_sp[sp_g] for sp_g in sp_genes]
                seq_len = {seq.id: len(seq.seq) for seq in seq_to_align}
                all_seq_len.update(seq_len)
                all_seq += seq_to_align

                # Create sequences fasta files in sequences directory
                seq_file = os.path.join(seq_output, f'{sp}__vs__{source_species}__{rxn}.fasta')
                with open(seq_file, "w") as output_handle:
                    SeqIO.write(list(seq_source) + seq_to_align, output_handle, "fasta")

                # Align all the sequences
                align_file = os.path.join(align_output, f'{sp}__vs__{source_species}__{rxn}.aln')
                muscle_cline = MuscleCommandline(input=seq_file, out=align_file)
                os.system(str(muscle_cline))

                # Write alignment in final file with all alignments
                self.__write_alignment(align_file, all_align_file, source_genes, all_seq_len)

        # Global align
        global_align_fasta = os.path.join(align_output, 'global_align.fasta')
        global_align_file = os.path.join(align_output, 'global_align.txt')
        seq_file = os.path.join(seq_output, 'all_seq.fasta')
        with open(seq_file, "w") as output_handle:
            SeqIO.write(list(seq_source) + all_seq, output_handle, "fasta")
        muscle_cline = MuscleCommandline(input=seq_file, out=global_align_fasta)
        os.system(str(muscle_cline))

        # Write alignment in final file with all alignments
        self.__write_alignment(global_align_fasta, global_align_file, source_genes, all_seq_len)

    @staticmethod
    def __write_alignment(fasta_file, txt_file, source_genes, all_seq_len):
        alignment = AlignIO.read(open(fasta_file), "fasta")
        with open(txt_file, 'a') as f:
            for record in alignment:
                if record.id in source_genes:
                    f.write(f'{record.seq}\t{record.id}\t{str(all_seq_len[record.id])}\n')
            for record in alignment:
                if record.id not in source_genes:
                    f.write(f'{record.seq}\t{record.id}\t{str(all_seq_len[record.id])}\n')

    @staticmethod
    def __create_directories(directories: Iterable[str]):
        for directory in directories:
            if not os.path.exists(directory):
                os.mkdir(directory)

    @staticmethod
    def __filter_genes(sp_genes):
        sp_genes_2 = set()
        for sp_g in sp_genes:
            if not sp_g.startswith('predicted_') and sp_g != 'nan':
                sp_genes_2.add(sp_g)
        return sp_genes_2