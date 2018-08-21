from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import *
from skbio.alignment import local_pairwise_align_ssw, make_identity_substitution_matrix
from skbio import TabularMSA, DNA, Protein
from .Parser import MasterFile, MasterCollection
from .utils import split_len
import os
# import time # TODO: Remove after tests

__all__ = ['MasterFile', 'StopChecker', 'MasterCollection']

# TODO: skbio.alignment.make_identity_substitution_matrix version 0.5.2 has open issue on substitution_matrix, it will be fixed in
# 0.5.3 which has not came up yet. It can be fixed after they fix the issue. (Maybe feed BLOSUM62)
# http://scikit-bio.org/docs/0.5.1/generated/skbio.alignment.make_identity_substitution_matrix.html#skbio.alignment.make_identity_substitution_matrix

class StopChecker:
    def __init__(self, mcol, outdir, kmer_length=15, threshold=10, match_score=10, mismatch_penalty=-4, submat=None, **kwargs):
        """Take a list of  collection of MasterFile as input and check stop codon"""
        # start_time = time.time() # TODO: Remove after tests
        self.mcol = mcol
        self.alignments(mcol, kmer_length, threshold, match_score, mismatch_penalty)
        self._give_meaning()
        # print("--- %s seconds --- For StopChecker" % (time.time() - start_time)) # TODO: Remove after tests

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
        submat = make_identity_substitution_matrix(5, -2, alphabet='ARNDCQEGHILKMFPSTWYVBZX*')
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
            data1_temporary.append(pairs_name)
            data2_temporary.append(pairs_results)

            # Scores_list will be used for doing necessary calculations with extensions.
            score_dictionary = {'gene_name': common_genes[i], 'original_species': non_reference_species, 'original_sequences': original_sequences, 'original_sequences_protein': original_sequences_protein, 'reference_specie': reference_species[i],
            'reference_sequences': reference_sequence, 'reference_sequence_protein': reference_sequence_protein, 'extensions': trimmed_extensions,
             'extensions_protein': allowed_extensions_protein, 'start_end': start_end}
            scores_list.append(score_dictionary)

        # TODO:TODO: Storing variable informations for modelling DT/RF.
        length_of_extensions = []
        length_of_genes = []
        frequency_of_stop_codon = []
        stop_codons_evolutionary = []
        for i in range(len(scores_list)): # For each gene // Genes should be taken apart seperately while collecting data?
            # The length of the extension
            # For the reference specie:
            length_of_extensions.append(0)
            # The length of the gene compared to other homologs:
            length_of_genes.append(len(scores_list[i]['reference_sequences']))
            for j in range(len(scores_list[i]['original_species'])): # For each non-reference-specie
                length_of_extensions.append(len(scores_list[i]['extensions'][j]))
                length_of_genes.append(len(scores_list[i]['original_sequences'][j]))
            # Frequency of the putative stop codon in all coding region
                # Divide into 3-mers, count how many stop codons occur in the all coding region (ignoring frameshift)
            substrings = split_len(scores_list[i]['reference_sequences'][::-1],3,0) # Reverse, divide 3-mers
            for j in range(len(substrings)):
                substrings[j] = substrings[j][::-1] # Reverse again each element. First element is the stop codon.
            countStopCodons = 0
            # Stop codon is first element, add to this specie's stop codon:
            stop_codons_evolutionary.append(substrings[0])
            for j in range(len(substrings)):
                if substrings[j] == substrings[0]:
                    countStopCodons += 1
            frequency_of_stop_codon.append(countStopCodons/len(substrings))
            # Obviously all of them have 1 because its a fucking stop codon? Emmanuel pls help.
            for j in range(len(scores_list[i]['original_species'])):
                substrings = split_len(scores_list[i]['original_sequences'][j][::-1],3,0) # Reverse, divide 3-mers
                for k in range(len(substrings)):
                    substrings[k] = substrings[k][::-1] # Reverse again each element. First element is the stop codon.
                countStopCodons = 0
                # Add into list of stop codons:
                stop_codons_evolutionary.append(substrings[0])
                for k in range(len(substrings)):
                    if substrings[k] == substrings[0]:
                        countStopCodons += 1
                frequency_of_stop_codon.append(countStopCodons/len(substrings))
        # Frequency of the stop codon in all evolutionary close stop codons
        frequency_evolutionary = []
        mean_similarity_extension = []
        for i in range(len(scores_list)):
            rev = scores_list[i]['reference_sequences'][::-1]
            stop_rev = rev[:3]
            stop = stop_rev[::-1]
            stopCodonCount = 0
            for j in range(len(stop_codons_evolutionary)):
                if stop_codons_evolutionary[j] == stop:
                    stopCodonCount += 1
            frequency_evolutionary.append(stopCodonCount/len(stop_codons_evolutionary))
            for j in range(len(scores_list[i]['original_species'])):
                stopCodonCount = 0
                rev = scores_list[i]['original_sequences'][j][::-1]
                stop_rev = rev[:3]
                stop = stop_rev[::-1]
                for k in range(len(stop_codons_evolutionary)):
                    if stop_codons_evolutionary[k] == stop:
                        stopCodonCount += 1
                frequency_evolutionary.append(stopCodonCount/len(stop_codons_evolutionary))
        # Mean similarity of the extension to other genomes
        # How to determine mean similarity of extensions with other genomes?
            # TODO: Store all scores of a non-original specie with others in local alignment. After that sum all of those scores and
            # divide by total number of species to get mean similarity.
            mean_similarity_extension.append(0) # Reference species do not have any extensions.

            for j in range(len(scores_list[i]['original_species'])):
                total_score_for_specie = 0
                for k in range(len(scores_list[i]['original_species'])):
                    if not j == k:
                        j_sequence = Protein(str(scores_list[i]['extensions_protein'][j]))
                        k_sequence = Protein(str(scores_list[i]['extensions_protein'][k]))
                        try:
                            try:
                                alignment, score, start_end_positions = local_pairwise_align_ssw(
                                j_sequence,
                                k_sequence,
                                substitution_matrix = submat,
                                )
                                total_score_for_specie += score
                            except IndexError:
                                total_score_for_specie += 0
                        except TypeError:
                            total_score_for_specie += 0
                mean_similarity_extension.append(total_score_for_specie / len(scores_list[i]['original_species']))

        # 1- length_of_extensions -> Stores length of extensions
        # 2- length_of_genes -> Stores length of genes
        # 3- frequency_of_stop_codon -> Stores frequency of stop codon in all coding region
        # 4- mean_similarity_extension -> Stores mean similarity of extension among all other genes.
        # TODO: 5- frequency_stop_codon_cdr -> Stores frequency of stop codon for each specie among all genes.
        # TODO: 6- mean_evolutionary_closeness -> Stores how close to other species depending on an evolutionary tree. TODO: Help me Emmanuel
        # Mean evolutionary closeness to other species depending on evolutionary tree fed?
        # For mean evolutionary closeness we need an evolutionary tree as an input.
        # Try to take a evolutionary tree as an input? Smh.
        # https://en.wikipedia.org/wiki/Newick_format -> Maybe we can use Newick tree formats?
        # TODO: Try to find more variables? TODO: Help me Emmanuel
        # mean_similarity_extension === extensions in information dictionary.
        # 7- mean_similarity_initials = initials in information dictionary.

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
                        refseq = Protein(str(reference_sequence_right))
                        try:
                            alignment, score, start_end_positions = local_pairwise_align_ssw(
                            newStringProtein,
                            refseq,
                            substitution_matrix = submat,
                            )
                            # Check if score is significant enough.
                            if not score < threshold:
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
                                start_end_position_alignments.append("No extensions")
                        except IndexError:
                            scores_alignments.append(0)
                            start_end_position_alignments.append("No extensions")
                else:
                    scores_alignments.append(0)
                    start_end_position_alignments.append("No extensions")
                best_start_end = "No extensions"
                best_score = 0

                for j in range(len(scores_alignments)):
                    if scores_alignments[j]/((j+1)**(1/2)) > best_score:
                        best_score = scores_alignments[j]
                        best_start_end = start_end_position_alignments[j]
                pairs_results_extensions.append(best_start_end)
            data3_temporary.append(pairs_results_extensions)

        # pairwise.
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
                            j_sequence = Protein(str(protein_j))
                            k_sequence = Protein(str(k_sequence_total))
                            try:
                                alignment, score, start_end_positions = local_pairwise_align_ssw(
                                j_sequence,
                                k_sequence,
                                substitution_matrix = submat,
                                )
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
                            except IndexError:
                                scores_alignments.append(0)
                                start_end_position_alignments.append("No extensions")
                    else:
                        scores_alignments.append(0)
                        start_end_position_alignments.append("No extensions")
                    best_start_end = "No extensions"
                    best_score = 0
                    for k in range(len(scores_alignments)):
                        if scores_alignments[k]/((k+1)**(1/2)) > best_score:
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
                            j_sequence = Protein(str(j_sequence_total))
                            k_sequence = Protein(str(protein_k))
                            try:
                                alignment, score, start_end_positions = local_pairwise_align_ssw(
                                j_sequence,
                                k_sequence,
                                substitution_matrix = submat,
                                )
                                # if not score < threshold:
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
                                # else:
                                #     scores_alignments.append(0)
                                #     start_end_position_alignments.append("None")
                            except IndexError:
                                scores_alignments.append(0)
                                start_end_position_alignments.append("No extensions")
                    else:
                        scores_alignments.append(0)
                        start_end_position_alignments.append("No extensions")
                    best_start_end = "No extensions"
                    best_score = 0
                    for k in range(len(scores_alignments)):
                        if scores_alignments[k]/((k+1)**(1/2)) > best_score:
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
            if totalPoints6 >= totalPoints7 and totalPoints6 >= len(data4_temporary[0])*threshold:
                for k in range(len(data4_temporary[i])):
                    data3_temporary[i].append(data4_temporary[i][k])
            elif totalPoints7 >= len(data4_temporary[0])*threshold:
                for k in range(len(data5_temporary[i])):
                    data3_temporary[i].append(data5_temporary[i][k])
            else:
                for k in range(len(data4_temporary[i])):
                    data3_temporary[i].append('No extensions')

        self.information_dictionary = {'pairs': data1_temporary, 'initials': data2_temporary, 'extensions': data3_temporary}
        return self.information_dictionary

    def _give_meaning(self):
        # Need to gather all the important data for all variables that are determined.
            # Variables required:
                # The length of the extension +
                # The length of the gene compared to other homologs +
                # The length of each gene's reference +

                # Frequency of the putative stop codon in all coding region
                # Frequency of the stop codon in all evolutionary close stop codons

                # Mean similarity of the extension to other genomes
                # Mean evolutionary closeness to other species depending on evolutionary tree fed?

                # CoreTracker uses -> Fisher's p value? Telford score of C coding for X?


        # Calculate their Gini impurities. (?) TODO: Help me Emmanuel
        # Manipulate data in a format such that training/test will be split easily.
        # Gather information with research to find some True Positive dataset. TODO: Help me Emmanuel
        # Create a Decision Tree Model for visualization.
        # Create a RF Model for meaning.
        pass
