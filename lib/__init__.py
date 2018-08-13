from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import *
from skbio.alignment import local_pairwise_align_ssw, make_identity_substitution_matrix
from skbio import TabularMSA, DNA, Protein
from Bio.SubsMat.MatrixInfo import blosum62
from .Parser import MasterFile, MasterCollection
from .utils import split_len
import os

__all__ = ['MasterFile', 'StopChecker', 'MasterCollection']


class StopChecker:
    def __init__(self, mcol, outdir, kmer_length=5, threshold=5, match_score=6, mismatch_penalty=-4, submat=None, **kwargs):
        """Take a list of  collection of MasterFile as input and check stop codon"""
        self.mcol = mcol
        self.alignments(mcol, kmer_length, threshold, match_score, mismatch_penalty)
        self._give_meaning()

    def alignments(self, mcol, kmer_length, threshold, match_score, mismatch_penalty):
        list_of_dictionaries, list_of_extended_dictionaries, reference_species, common_genes = mcol.digest()
        list_of_formatted_dictionaries = []
        # finds common genes between all species.
        # for each specie
        for i in range(len(list_of_dictionaries)):
            genes_of_interest = []
            sequences = []
            extended_sequences = []
            # for each gene
            for j, gene in enumerate(list_of_dictionaries[i]['genes_list']):
                for k in range(len(common_genes)):
                    if list_of_dictionaries[i]['genes_list'][j] == common_genes[k]:
                        genes_of_interest.append(
                            list_of_dictionaries[i]['genes_list'][j])
                        sequences.append(
                            list_of_dictionaries[i]['sequences'][gene])
            formatted_dict = {'species_name': list_of_dictionaries[i]['species_name'],
                              'species_gc': list_of_dictionaries[i]['species_gc'], 'genes_list': genes_of_interest, 'sequences': sequences}
            list_of_formatted_dictionaries.append(formatted_dict)
        
        # Does a local alignment between reference specie's genes and others.
        complete_gene_list = []
        data1_temporary = []
        data2_temporary = []
        data3_temporary = []
        data4_temporary = [] # Start end for pairwise extension j-slide 39 elements
        data5_temporary = [] # Scores for pairwise extension j-slide 39 elements
        data6_temporary = [] # Start end for pairwise extension k-slide 39 elements
        data7_temporary = [] # Scores for pairwise extension k-slide 39 elements
        data8_temporary = [] # Best alignments of j-k slides.
        scores_list = []
        # TODO: Why doesn't this work ? submat = _load_matrix("blosum62.txt")
        submat = make_identity_substitution_matrix(5, -1, alphabet='ARNDCQEGHILKMFPSTWYVBZX*')
        # List includes a dictionary for each gene
        for i in range(len(common_genes)):
            reference_specie_name = ""
            reference_sequence = ""
            reference_sequence_protein = ""
            non_reference_species = []
            genes = []
            original_species_names = []
            original_sequences = []
            original_sequences_protein = []
            non_reference_species_gc = []
            extensions = []
            for j in range(len(list_of_formatted_dictionaries)):
                # Finds the sequence for the reference specie:
                if list_of_formatted_dictionaries[j]['species_name'] == reference_species[i]:
                    reference_specie_name = list_of_formatted_dictionaries[j]['species_name']
                    for k in range(len(list_of_formatted_dictionaries[j]['genes_list'])):
                        if list_of_formatted_dictionaries[j]['genes_list'][k] == common_genes[i]:
                        # Transforms nucleotide sequences into protein sequences
                            # Removes the '!' chars in sequence for proper transformation to proteins.
                            formatted = list_of_formatted_dictionaries[j]['sequences'][k].replace("!","")
                            # Finds id of the specie we are looking at:
                            table_id = list_of_formatted_dictionaries[j]['species_gc']
                            # Inserts the id into table to get accurate transform table:
                            gene = Seq(formatted,generic_dna)
                            # Uses the transform table to transform nucleotide sequence into protein sequence.
                            protein = gene.translate(table = table_id)
                            reference_sequence = formatted
                            # Stores the protein representation of reference sequence.
                            reference_sequence_protein = protein
                            break
                # Finds non-reference sequences for that gene:
                else:
                    for k in range(len(list_of_formatted_dictionaries[j]['genes_list'])):
                        if list_of_formatted_dictionaries[j]['genes_list'][k] == common_genes[i]:
                            # Transforms nucleotide sequences into protein sequences
                            formatted = list_of_formatted_dictionaries[j]['sequences'][k].replace("!","")
                            # Finds id of the specie we are looking at:
                            table_id = list_of_formatted_dictionaries[j]['species_gc']
                            # Inserts the id into table to get accurate transform table:
                            gene = Seq(formatted,generic_dna)
                            protein = gene.translate(table = table_id)
                            non_reference_species.append(list_of_formatted_dictionaries[j]['species_name'])
                            non_reference_species_gc.append(table_id)
                            genes.append(list_of_formatted_dictionaries[j]['genes_list'][k])
                            original_sequences.append(formatted)
                            original_sequences_protein.append(protein)
                            break

            # pairs_name: Stores the pairs in a string format-> Bracteacoccus_minor / Bracteacoccus_aerius i.e
            pairs_name = []
            # pairs_result: is a list of lists. for each pairs_name, it stores score / start_end / alignments.
            pairs_results = []
            # start_end: stores start_end results
            start_end = []
            # Traverses non-reference_sequences and do local alignment with reference, stores their scores.
            for j in range(len(non_reference_species)):
                # results: a list that stores scores, start_end and alignments.
                # results = []
                # ! Put pairs in correct order with start-end order.
                pairs_name.append("%s/%s" % (non_reference_species[j],reference_specie_name))
                original_string = Protein(str(original_sequences_protein[j]))
                reference_string = Protein(str(reference_sequence_protein))
                alignment, score, start_end_positions = local_pairwise_align_ssw(
                original_string,
                reference_string,
                substitution_matrix = submat,
                )
                #results.append(start_end_positions)

                pairs_results.append(start_end_positions)
                start_end.append(start_end_positions)

            # Now do pairwise between all original sequences:
            for j in range(len(non_reference_species)-1):
                for k in range(j+1,len(non_reference_species)):
                    # results = []
                    pairs_name.append("%s/%s" % (non_reference_species[j],non_reference_species[k]))
                    j_sequence = Protein(str(original_sequences_protein[j]))
                    k_sequence = Protein(str(original_sequences_protein[k]))
                    alignment, score, start_end_positions = local_pairwise_align_ssw(
                    j_sequence,
                    k_sequence,
                    substitution_matrix = submat,
                    suppress_sequences = True
                    )
                    # results.append(start_end_positions)
                    pairs_results.append(start_end_positions)
                    start_end.append(start_end_positions)

            # Finds extensions
            for j in range(len(non_reference_species)):
                for k in range(len(list_of_extended_dictionaries)):
                    count = 0
                    if list_of_extended_dictionaries[k]['species_name'] == non_reference_species[j]:
                        for l in range(len(list_of_extended_dictionaries[k]['genes_list'])):
                            if genes[j] == list_of_extended_dictionaries[k]['genes_list'][l]:
                                if count == 1:
                                    break
                                extensions.append(list_of_extended_dictionaries[k]['sequences'][l])
                                count = count + 1

            # Depending on parameter extension_length, trims extensions if necessary
            trimmed_extensions = []
            for j in range(len(extensions)):
                length_of_reference = len(reference_sequence)
                length_of_extension = len(extensions[j])
                length_of_original_sequence = len(original_sequences[j])
                length_of_possible_extension = length_of_reference - length_of_original_sequence
                if  length_of_extension > length_of_possible_extension:
                    trimmed_extension = ""
                    # Maximum extension can get at reference sequence's length.
                    trimmed_extension = extensions[j][:length_of_possible_extension]
                    trimmed_extensions.append(trimmed_extension.replace("!",""))
                else:
                    trimmed_extensions.append(extensions[j].replace("!",""))

            allowed_extensions_protein = []
            for j in range(len(trimmed_extensions)):
                table_id = non_reference_species_gc[j]
                gene = Seq(trimmed_extensions[j],generic_dna)
                protein = gene.translate(table = table_id)

                allowed_extensions_protein.append(protein)
            # TODO: TODO: Order is weird you can fix it
            data1_temporary.append(pairs_name)
            data2_temporary.append(pairs_results)

            # Scores_list will be used for doing necessary calculations with extensions.
            score_dictionary = {'gene_name': common_genes[i], 'original_species': non_reference_species, 'original_sequences': original_sequences, 'original_sequences_protein': original_sequences_protein, 'reference_specie': reference_species[i],
            'reference_sequences': reference_sequence, 'reference_sequence_protein': reference_sequence_protein, 'extensions': trimmed_extensions,
             'extensions_protein': allowed_extensions_protein, 'start_end': start_end}
            scores_list.append(score_dictionary)

        # Here we will start a new, final chapter everyone!
        # First compare all sequences with reference sequence.
        for i in range(len(scores_list)):
            pairs_results_extensions = []
            # Pairwise reference with each non-reference specie.
            for j in range(len(scores_list[i]['original_species'])):
                length_of_extension = len(scores_list[i]['extensions_protein'][j])
                # Divide extensions to k-mers.
                division_amount = 0
                if length_of_extension < kmer_length:
                    division_amount = 1
                else:
                    division_amount = int(length_of_extension/kmer_length)
                # Find the place left in reference sequence.
                end_position_reference = scores_list[i]['start_end'][j][1][1]
                # Find end position of original sequence.
                end_position_original = scores_list[i]['start_end'][j][0][1]
                # Take whats left of reference from initial local alignment.
                reference_sequence_right = scores_list[i]['reference_sequence_protein'][end_position_reference:]
                # Concatenate the last part of original sequence and the extensions for them.
                concatenated_original_seq = scores_list[i]['original_sequences_protein'][j][end_position_original:]
                concatenated = concatenated_original_seq + scores_list[i]['extensions_protein'][j]
                # Length of each k-mer
                amount = int(length_of_extension/division_amount)
                # List for storing k-mers
                divided_extensions = []
                # Store scores to find maximum scored one, preferably closest to the end of
                # initial local alignment.
                scores_alignments = []
                start_end_position_alignments = []
                if not length_of_extension == 0:
                    divided_extensions = split_len(concatenated, amount,0)
                    for k in range(len(divided_extensions)):
                        newStringProtein = Protein(str(divided_extensions[k]))
                        try:
                            refseq = Protein(str(reference_sequence_right))
                            try:
                                alignment, score, start_end_positions = local_pairwise_align_ssw(
                                newStringProtein,
                                refseq,
                                substitution_matrix = submat,
                                )
                                # Check if score is significant enough.
                                if score > threshold:
                                    scores_alignments.append(score)
                                    # Calculate start end position with respect to beginning.
                                    # start_location = start of the extension + other extensions before * their length + end of local alignment
                                    start_location = start_end_positions[0][0] + len(divided_extensions[0])*k + scores_list[i]['start_end'][j][0][1]
                                    # start_location_r = end location of initial alignment + start position
                                    start_location_r = start_end_positions[1][0] + scores_list[i]['start_end'][j][1][1]
                                    # end_location = end of the extension + other extensions before * their length + end of local alignment
                                    end_location = start_end_positions[0][1] + len(divided_extensions[0])*k + scores_list[i]['start_end'][j][0][1]
                                    end_location_r = start_end_positions[1][1] + scores_list[i]['start_end'][j][1][1]
                                    start_end_position_alignments.append("(%s,%s),(%s,%s)" % (start_location,end_location,start_location_r,end_location_r))
                                else:
                                    scores_alignments.append(0)
                                    start_end_position_alignments.append("None")
                            except IndexError:
                                scores_alignments.append(0)
                                start_end_position_alignments.append("None")
                        except ValueError:
                            scores_alignments.append(0)
                            start_end_position_alignments.append("None")
                else:
                    scores_alignments.append(0)
                    start_end_position_alignments.append("None")
                best_start_end = "None"
                best_score = 0

                for j in range(len(scores_alignments)):
                    if scores_alignments[j]/((j+1)**(1/2)) > best_score:
                        best_score = scores_alignments[j]
                        best_start_end = start_end_position_alignments[j]
                pairs_results_extensions.append(best_start_end)
            data3_temporary.append(pairs_results_extensions)

        # pairwise.
        # TODO: Problem with extensions is probably here.
        for i in range(len(scores_list)):
            pairs_results_extensions_j = []
            pairs_results_extensions_j_score = []
            for j in range(len(scores_list[i]['original_species'])-1):
                for k in range(j+1,len(scores_list[i]['original_species'])):
                    length_of_extension_j = len(scores_list[i]['extensions_protein'][j])
                    length_of_extension_k = len(scores_list[i]['extensions_protein'][k])
                    # Divide extensions to k-mers.
                    division_amount_j = 1
                    division_amount_k = 1
                    if not length_of_extension_j < kmer_length:
                        division_amount_j = int(length_of_extension_j/kmer_length)
                    if not length_of_extension_k < kmer_length:
                        division_amount_k = int(length_of_extension_k/kmer_length)

                    # Find the place left in reference sequence
                    orig_spec = len(scores_list[i]['original_species'])
                    j_loc = j + orig_spec
                    # TODO: Check this calculation, possibly wrong.
                    end_position_j = scores_list[i]['start_end'][j_loc][0][1]
                    end_position_k = scores_list[i]['start_end'][j_loc][1][1]
                    # Take whats left of reference from initial local alignment.
                    j_sequence_right = scores_list[i]['original_sequences_protein'][j][end_position_j:]
                    k_sequence_right = scores_list[i]['original_sequences_protein'][k][end_position_k:]
                    j_sequence_total = j_sequence_right + scores_list[i]['extensions_protein'][j]
                    k_sequence_total = k_sequence_right + scores_list[i]['extensions_protein'][k]
                    # Length of each k-mer
                    amount_j = int(length_of_extension_j/division_amount_j)
                    amount_k = int(length_of_extension_k/division_amount_k)

                    # List for storing k-mers
                    divided_extensions_j = []
                    divided_extensions_k = []
                    # Store scores to find maximum scored one, preferably closest to the end of
                    # initial local alignment.
                    scores_alignments = []
                    start_end_position_alignments = []
                    if not length_of_extension_k == 0 and not length_of_extension_j == 0:
                        divided_extensions_j = split_len(j_sequence_total, amount_j,0)
                        divided_extensions_k = split_len(k_sequence_total, amount_k,0)
                        # First slide j
                        for l in range(len(divided_extensions_j)):
                            # Now we have to extend them both in the double for loop and collect scores.
                            protein_j = Protein(str(divided_extensions_j[l]))
                            # protein_k = Protein(str(divided_extensions_k))
                            try:
                                j_sequence = Protein(str(protein_j))
                                k_sequence = Protein(str(k_sequence_total))
                                try:
                                    alignment, score, start_end_positions = local_pairwise_align_ssw(
                                    j_sequence,
                                    k_sequence,
                                    substitution_matrix = submat,
                                    )
                                    if score > threshold:
                                        scores_alignments.append(score)
                                        # Calculate start end position with respect to beginning.
                                        # start_location = start of the extension + other extensions before * their length + end of local alignment
                                        start_location = start_end_positions[0][0] + len(divided_extensions_j[0])*l + scores_list[i]['start_end'][j][0][1]
                                        # start_location_r = end location of initial alignment + start position
                                        start_location_r = start_end_positions[1][0] + scores_list[i]['start_end'][j][1][1]
                                        # end_location = end of the extension + other extensions before * their length + end of local alignment
                                        end_location = start_end_positions[0][1] + len(divided_extensions_j[0])*l + scores_list[i]['start_end'][j][0][1]
                                        end_location_r = start_end_positions[1][1] + scores_list[i]['start_end'][j][1][1]
                                        start_end_position_alignments.append("(%s,%s),(%s,%s)" % (start_location,end_location,start_location_r,end_location_r))
                                    else:
                                        scores_alignments.append(0)
                                        start_end_position_alignments.append("None")
                                except IndexError:
                                    scores_alignments.append(0)
                                    start_end_position_alignments.append("None")
                            except ValueError:
                                scores_alignments.append(0)
                                start_end_position_alignments.append("None")
                    else:
                        scores_alignments.append(0)
                        start_end_position_alignments.append("None")
                    best_start_end = "None"
                    best_score = 0
                    for k in range(len(scores_alignments)):
                        if scores_alignments[k]/((k+1)**(1/3)) > best_score:
                            best_score = scores_alignments[k]
                            best_start_end = start_end_position_alignments[k]
                    pairs_results_extensions_j_score.append(best_score)
                    pairs_results_extensions_j.append(best_start_end)
            data4_temporary.append(pairs_results_extensions_j)
            data6_temporary.append(pairs_results_extensions_j_score)

        for i in range(len(scores_list)):
            pairs_results_extensions_k = []
            pairs_results_extensions_k_score = []
            for j in range(len(scores_list[i]['original_species'])-1):
                for k in range(j+1,len(scores_list[i]['original_species'])):
                    length_of_extension_j = len(scores_list[i]['extensions_protein'][j])
                    length_of_extension_k = len(scores_list[i]['extensions_protein'][k])
                    # Divide extensions to k-mers.
                    division_amount_j = 1
                    division_amount_k = 1
                    if not length_of_extension_j < kmer_length:
                        division_amount_j = int(length_of_extension_j/kmer_length)
                    if not length_of_extension_k < kmer_length:
                        division_amount_k = int(length_of_extension_k/kmer_length)

                    # Find the place left in reference sequence
                    orig_spec = len(scores_list[i]['original_species'])
                    j_loc = j + orig_spec
                    # TODO: Check this calculation, possibly wrong.
                    end_position_j = scores_list[i]['start_end'][j_loc][0][1]
                    end_position_k = scores_list[i]['start_end'][j_loc][1][1]
                    # Take whats left of reference from initial local alignment.
                    j_sequence_right = scores_list[i]['original_sequences_protein'][j][end_position_j:]
                    k_sequence_right = scores_list[i]['original_sequences_protein'][k][end_position_k:]
                    j_sequence_total = j_sequence_right + scores_list[i]['extensions_protein'][j]
                    k_sequence_total = k_sequence_right + scores_list[i]['extensions_protein'][k]
                    # Length of each k-mer
                    amount_j = int(length_of_extension_j/division_amount_j)
                    amount_k = int(length_of_extension_k/division_amount_k)

                    # List for storing k-mers
                    divided_extensions_j = []
                    divided_extensions_k = []
                    # Store scores to find maximum scored one, preferably closest to the end of
                    # initial local alignment.
                    scores_alignments = []
                    start_end_position_alignments = []
                    if not length_of_extension_k == 0 and not length_of_extension_j == 0:
                        divided_extensions_j = split_len(j_sequence_total, amount_j,0)
                        divided_extensions_k = split_len(k_sequence_total, amount_k,0)
                        # Second slide k
                        for l in range(len(divided_extensions_k)):
                            # Now we have to extend them both in the double for loop and collect scores.
                            protein_k = Protein(str(divided_extensions_k[l]))
                            # protein_k = Protein(str(divided_extensions_k))
                            try:
                                j_sequence = Protein(str(j_sequence_total))
                                k_sequence = Protein(str(protein_k))
                                try:
                                    alignment, score, start_end_positions = local_pairwise_align_ssw(
                                    j_sequence,
                                    k_sequence,
                                    substitution_matrix = submat,
                                    )
                                    if score > threshold:
                                        scores_alignments.append(score)
                                        # Calculate start end position with respect to beginning.
                                        # start_location = start of the extension + other extensions before * their length + end of local alignment
                                        start_location = start_end_positions[0][0] + len(divided_extensions_j[0])*l + scores_list[i]['start_end'][j][0][1]
                                        # start_location_r = end location of initial alignment + start position
                                        start_location_r = start_end_positions[1][0] + scores_list[i]['start_end'][j][1][1]
                                        # end_location = end of the extension + other extensions before * their length + end of local alignment
                                        end_location = start_end_positions[0][1] + len(divided_extensions_j[0])*l + scores_list[i]['start_end'][j][0][1]
                                        end_location_r = start_end_positions[1][1] + scores_list[i]['start_end'][j][1][1]
                                        start_end_position_alignments.append("(%s,%s),(%s,%s)" % (start_location,end_location,start_location_r,end_location_r))
                                    else:
                                        scores_alignments.append(0)
                                        start_end_position_alignments.append("None")
                                except IndexError:
                                    scores_alignments.append(0)
                                    start_end_position_alignments.append("None")
                            except ValueError:
                                scores_alignments.append(0)
                                start_end_position_alignments.append("None")
                    else:
                        scores_alignments.append(0)
                        start_end_position_alignments.append("None")
                    best_start_end = "None"
                    best_score = 0
                    for k in range(len(scores_alignments)):
                        if scores_alignments[k]/((k+1)**(1/3)) > best_score:
                            best_score = scores_alignments[k]
                            best_start_end = start_end_position_alignments[k]
                    pairs_results_extensions_k_score.append(best_score)
                    pairs_results_extensions_k.append(best_start_end)
            data5_temporary.append(pairs_results_extensions_k)
            data7_temporary.append(pairs_results_extensions_k_score)

        for i in range(len(data6_temporary)):
            totalPoints6 = 0
            totalPoints7 = 0
            for j in range(len(data6_temporary[i])):
                totalPoints6 += data6_temporary[i][j]
                totalPoints7 += data7_temporary[i][j]
            if totalPoints6 >= totalPoints7:
                for k in range(len(data4_temporary[i])):
                    data3_temporary[i].append(data4_temporary[i][k])
            else:
                for k in range(len(data5_temporary[i])):
                    data3_temporary[i].append(data5_temporary[i][k])

        self.information_dictionary = {'pairs': data1_temporary, 'initials': data2_temporary, 'extensions': data3_temporary}
        return self.information_dictionary

    def _give_meaning(self):

        for i in range(len(self.information_dictionary['pairs'])):
            print("Gene: %s" % self.mcol.common_genes[i]) # For each gene.
            print("")
            for j in range(len(self.information_dictionary['pairs'][i])): # For each pair for that gene
                print("Pair: %s" % self.information_dictionary['pairs'][i][j])
                print("Initial: %s" % self.information_dictionary['initials'][i][j])
                print("Extensions: %s" % self.information_dictionary['extensions'][i][j])
                print("")