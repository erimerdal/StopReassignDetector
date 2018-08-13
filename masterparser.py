import re
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from skbio.alignment import local_pairwise_align_ssw,make_identity_substitution_matrix
from skbio import TabularMSA, DNA, Protein
from pandas import *
import sys
import os

class MasterFile:

    # Master File format for sequence alignment.
    # This class helps parsing Master File Formats.
    # It will be used in the cases that we will also be looking for
    # stop - sense codon reassignments.
    #
    # The MasterFile format is the following:
    # ;; Masterfile modified automatically by mfannot version 1.33
    # ;; on Tue Dec  5 01:45:11 2017 by user enoutahie on host cardamom
    # ;;    - Gene Totals: 63
    # ;;    - List of genes added:
    # ;;      atp6                 atp9                 cob (2 introns)
    # ...
    # ...
    # >Bracteacoccus_minor gc=22
    # ;     G-orf615 <== end
    #   1234  TGACACATAACAATAATAAACAGAAACAAGATTTGT
    # ;     G-nad5 ==> start
    #   1270  ATGTATTTCTGTGCAGTAATGGCGCCTTTCCTAGGCAGCGCACTGGCTGGATTGGCCGGA
    #   ...
    # ;     G-orf615 <== start
    #   3082  GACTACGCTCTAGTACTCAAAGTATGCGTGGTCTTCGGTCTGATGGCACTCGCTTTCCCA
    #   ...
    # ;     G-nad5 ==> end

    def __init__(self, infile, kmer_length = 15, threshold = 0, gap_open_penalty = -2, match_score = 6, mismatch_penalty = -4):
        # Parameters:
            # Infile parameter should be a list of Master Files obtained from MFannot results.
                # TODO: Maybe we can include in our code MFannot command line tools so a person can input as fasta.
            # kmer_length: An integer which decides how many nucleotides/extension block should be.
                # Increasing kmer_length will make it harder to find a better extension.
            # threshold: An integer that decides the minimum score that local alignments of extensions can return.
            # gap_open_penalty, match_score, mismatch_penalty: Not used currently because of Blosum62 Matrix scores.

        self.infile = infile
        list_of_dictionaries = []
        species_names = []
        reference_species = []
        common_genes = []
        # Constructor collects all sequences necessary from the input files and stores them.
        for x in range(len(infile)):
            self.sequences = self._internal_masterparser(infile[x])
            list_of_dictionaries.append(self.sequences)
            species_names.append(self.sequences['species_name'])
        # Constructor collects all the extensions to those sequences and stores them.
        # Also collects the species to take reference, and genes that are in common.
        common_genes, reference_species, list_of_extended_dictionaries = self._internal_extendedparser(infile,list_of_dictionaries)
        # Constructor sends all collected information to align and give a meaning to the results.
        information_dictionary = self._alignments(list_of_dictionaries,list_of_extended_dictionaries, kmer_length, threshold, gap_open_penalty, match_score, mismatch_penalty, reference_species, common_genes)
        self._give_meaning(information_dictionary,common_genes)

    def _alignments(self, list_of_dictionaries,list_of_extended_dictionaries, kmer_length, threshold, gap_open_penalty, match_score, mismatch_penalty, reference_species, common_genes):
        list_of_formatted_dictionaries = []
        # finds common genes between all species.
        # for each specie
        for i in range(len(list_of_dictionaries)):
            genes_of_interest = []
            sequences = []
            extended_sequences = []
            # for each gene
            for j in range(len(list_of_dictionaries[i]['genes_list'])):
                for k in range(len(common_genes)):
                    if list_of_dictionaries[i]['genes_list'][j] == common_genes[k]:
                        genes_of_interest.append(list_of_dictionaries[i]['genes_list'][j])
                        sequences.append(list_of_dictionaries[i]['sequences'][j])
            formatted_dict = {'species_name': list_of_dictionaries[i]['species_name'], 'species_id': list_of_dictionaries[i]['species_id'], 'genes_list': genes_of_interest, 'sequences': sequences}
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
            non_reference_species_id = []
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
                            table_id = list_of_formatted_dictionaries[j]['species_id']
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
                            table_id = list_of_formatted_dictionaries[j]['species_id']
                            # Inserts the id into table to get accurate transform table:
                            gene = Seq(formatted,generic_dna)
                            protein = gene.translate(table = table_id)
                            non_reference_species.append(list_of_formatted_dictionaries[j]['species_name'])
                            non_reference_species_id.append(table_id)
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
                suppress_sequences = True
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
                table_id = non_reference_species_id[j]
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
            if totalPoints6 >= totalPoints7:
                for k in range(len(data4_temporary[i])):
                    data3_temporary[i].append(data4_temporary[i][k])
            else:
                for k in range(len(data5_temporary[i])):
                    data3_temporary[i].append(data5_temporary[i][k])

        information_dictionary = {'pairs': data1_temporary, 'initials': data2_temporary, 'extensions': data3_temporary}
        return information_dictionary

    def _give_meaning(self,information_dictionary,common_genes):
        pass
        # print(len(information_dictionary['pairs'])) = 13
        # print(len(information_dictionary['initials'])) = 13
        # print(len(information_dictionary['extensions']))
        # for i in range(len(information_dictionary['pairs'])):
        #     print("Gene: %s" % common_genes[i]) # For each gene.
        #     print("")
        #     for j in range(len(information_dictionary['pairs'][i])): # For each pair for that gene
        #         print("Pair: %s" % information_dictionary['pairs'][i][j])
        #         print("Initial Local Alignment: %s" % information_dictionary['initials'][i][j])
        #         print("")
        #
        # for i in range(len(information_dictionary['extensions'])):
        #     print(information_dictionary['extensions'][i])
        #         #Using only initial local alignment.


    def _internal_extendedparser(self, infile, list_of_dictionaries):
         # For each gene that is common, we need to find if the gene is extendable or not.
         #     1- We need to find which genes are common in all species, and store them.
         #     2- Filter the interesting genes
         # List of interesting genes: nad / atp / cob / cox / sdh / tat can be extended. Other genes are omitted.
        common_genes_3x = []
        for x in range(len(list_of_dictionaries[0]['genes_list'])):
            for y in range(1,len(list_of_dictionaries)):
                for z in range(len(list_of_dictionaries[y]['genes_list'])):
                    if list_of_dictionaries[0]['genes_list'][x] == list_of_dictionaries[y]['genes_list'][z]:
                        if "G-nad" in list_of_dictionaries[0]['genes_list'][x] or "G-atp" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(list_of_dictionaries[0]['genes_list'][x])
                        if "G-cob" in list_of_dictionaries[0]['genes_list'][x] or "G-cox" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(list_of_dictionaries[0]['genes_list'][x])
                        if "G-sdh" in list_of_dictionaries[0]['genes_list'][x] or "G-tat" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(list_of_dictionaries[0]['genes_list'][x])
        common_genes = unique(common_genes_3x)
        # 3- Now that we have genes that are common and interesting, among all files we need to find
        #      the longest sequenced ones for that specific gene, to use it as a reference gene to compare
        #      with others later on.
        reference_species = []
        sequence_lengths = []
        for x in range(len(common_genes)):
            for y in range(len(infile)):
                for z in range(len(list_of_dictionaries[y]['genes_list'])):
                    if list_of_dictionaries[y]['genes_list'][z] == common_genes[x]:
                        sequence_lengths.append(len(list_of_dictionaries[y]['sequences'][z]))

        # 4- Now find the longest of these for each gene and store which specie it is.
        #      4.a- Divide in chunks of size species.
        chunks = [sequence_lengths[x*len(infile):x*len(infile)+len(infile)] for x in range(0, int(len(sequence_lengths)/len(infile)))]
        #      4.b- Find longest for each and store which specie has longest.
        # TODO: Change this to dynamic finding later. Sorry for this..
        # TODO: Reference_specie is automatically Bracteacoccus_minor now for every gene since its the only gc = 22.
        for x in range(len(chunks)):
            maximum = max(chunks[x])
            for y in range(len(chunks[x])):
                if maximum == chunks[x][y]:
                    reference_species.append(list_of_dictionaries[y]['species_name'])
                    break
        # for y in range(len(common_genes)):
        #     for x in range(len(list_of_dictionaries)):
        #         if "minor" in list_of_dictionaries[x]['species_name']:
        #             reference_species.append(list_of_dictionaries[x]['species_name'])
        #             break

        # 5- For each gene except the longest stored one, find extensions and store them
        list_of_extended_dictionaries = []
        reverse_gene_order = []
        line_numbers = []
        # First parse whole file to find which gene has <== start before our gene.
        # This for loop stores the reverse genes.
        for x in range(len(infile)):
            line_number = 0
            initial_run_reverse = []
            line_number_array = []
            currentWD = os.getcwd()
            handle = open(currentWD + "/" + infile[x], 'r')
            line = handle.readline()
            line_number = line_number + 1;
            while True:
                if not line:
                    break
                if line.startswith(';') and not line.startswith(';;'):
                    # Found a gene
                    geneName = line[1:].rstrip()
                    fullgenename = " ".join(geneName.split())
                    iterate = fullgenename.split(" ")
                    # Necessary informations are parsed
                    genename = iterate[0]
                    updown = iterate[1]
                    startend = iterate[2]
                    # We only collect reverse starts and ends.
                    if updown == "<==":
                        if startend == "start":
                            initial_run_reverse.append(genename)
                            line_number_array.append(line_number)
                            line = handle.readline()
                            line_number = line_number + 1;
                        else:
                            line = handle.readline()
                            line_number = line_number + 1;
                    else:
                        line = handle.readline()
                        line_number = line_number + 1;
                else:
                    line = handle.readline()
                    line_number = line_number + 1;
            reverse_gene_order.append(initial_run_reverse)
            line_numbers.append(line_number_array)

        # This is a big loop for reading extension sequences.
        # For all species (infile):
        for x in range(len(infile)):
            extended_sequences = []
            genes = []
            # For all common genes that these species have:
            for y in range(len(common_genes)):
                if not reference_species[y] == list_of_dictionaries[x]['species_name']:
                    currentWD = os.getcwd()
                    handle = open(currentWD + "/" + infile[x], 'r')
                    line = handle.readline()
                    start_end_count = 0
                    extension = []
                    reverse_extension = []
                    while True:
                        if not line:
                            break
                        if line.startswith(';') and not line.startswith(';;'):
                            geneName = line[1:].rstrip()
                            fullgenename = " ".join(geneName.split())
                            iterate = fullgenename.split(" ")
                            # Necessary informations are parsed
                            genename = iterate[0]
                            updown = iterate[1]
                            startend = iterate[2]
                            # We have found our gene
                            if genename == common_genes[y]:
                                # For the reverse genes we have a different approach than non-reverse genes.
                                if updown == "<==":
                                    if startend == "end":
                                        genes.append(genename)
                                        currentWD = os.getcwd()
                                        reverse_handle = open(currentWD + "/" + infile[x], 'r')
                                        from_line = 0
                                        to_line = 0
                                        for i in range(len(reverse_gene_order[x])):
                                            if reverse_gene_order[x][i] == genename:
                                                to_line = line_numbers[x][i]
                                                from_line = line_numbers[x][i-1]
                                                break

                                        for i in range(from_line):
                                            reverse_line = reverse_handle.readline()

                                        for i in range(from_line,to_line):
                                            string = genename + " "
                                            if string in reverse_line:
                                                break
                                            reverse_extension.append(reverse_line.rstrip())
                                            if reverse_line.startswith(';;') or reverse_line.startswith(';'):
                                                reverse_extension.remove(reverse_line.rstrip())
                                            reverse_line = reverse_handle.readline()
                                        sequence = "".join(reverse_extension).replace(" ", "").replace("\r", "")
                                        sequence = re.sub("\d","",sequence)
                                        # Reversing extensions seem like a better option
                                        sequence_reversed = reverse_sequence(sequence)
                                        extended_sequences.append(sequence_reversed[::-1])

                                start_end_count = start_end_count + 1;

                                # This is for the non-reverse gene extensions
                                if start_end_count == 2:
                                    if startend == "end":
                                        # Go until another gene starts
                                        line = handle.readline()
                                        genes.append(genename)
                                        while True:
                                            if line.startswith(';'):
                                                geneName = line[1:].rstrip()
                                                fullgenename = " ".join(geneName.split())
                                                iterate = fullgenename.split(" ")
                                                # Necessary informations are parsed
                                                genename = iterate[0]
                                                updown = iterate[1]
                                                startend = iterate[2]
                                                if startend == "start":
                                                    if updown == "==>":
                                                        break;
                                                    else:
                                                        line = handle.readline()
                                                else:
                                                    line = handle.readline()

                                            if not line:
                                                break
                                            else:
                                                extension.append(line.rstrip())
                                                # Removes other genes starting or end.
                                                if line.startswith(';;') or line.startswith(';'):
                                                    extension.remove(line.rstrip())
                                                line = handle.readline()
                                        # Format the sequence in a better, readable format.
                                        sequence = "".join(extension).replace(" ", "").replace("\r", "")
                                        sequence = re.sub("\d","",sequence)

                                        extended_sequences.append(sequence)

                                    else:
                                        line = handle.readline()
                                else:
                                    line = handle.readline()
                            else:
                                line = handle.readline()
                        else:
                            line = handle.readline()

            # Add all the necessary information to a dictionary to return.
            extended_dict = {'species_name': list_of_dictionaries[x]['species_name'] , 'genes_list': genes, 'sequences': extended_sequences}

            # Add each dictionary / per specie to a list of dictionaries.
            list_of_extended_dictionaries.append(extended_dict)
        return common_genes, reference_species, list_of_extended_dictionaries


    def _internal_masterparser(self, infile):
         currentWD = os.getcwd()
         handle = open(currentWD + "/" + infile, 'r')
         line = handle.readline()
         while True:
             if not line:
                 break
             while True:
                 line = handle.readline()
                 if not line:
                     break
                 if line.startswith('>'):
                     break
             speciesname = line[1:].rstrip()
             fullspeciesname = " ".join(speciesname.split())
             iterate = fullspeciesname.split(" ")
             specie_id_unformatted = iterate[1]
             # specie_id gives 22 if gc='22', 1 if gc='1' etc.
             specie_id = specie_id_unformatted[3:]
             reverse = []
             genes = []
             sequences = []
             line = handle.readline()
             # Storing genes and if they should be reversed.
             while True:
                 if not line:
                     break
                 if line.startswith(';') and not line.startswith(';;'):
                     geneName = line[1:].rstrip()
                     fullgenename = " ".join(geneName.split())
                     iterate = fullgenename.split(" ")
                     # Necessary informations are parsed
                     genename = iterate[0]
                     updown = iterate[1]
                     startend = iterate[2]
                     # We should store the gene in genes with these conditions:
                     #     1- If gene name has ==> and start in it
                     #     2- If gene name has <== and end in it, then reverse it.
                     # We will be removing introns and exons from gene names.
                     if not "-I" in genename:
                         if not "-E" in genename:
                             if updown == "==>" and startend == "start":
                                genes.append(genename)
                                reverse.append(False)
                             if updown == "<==" and startend == "end":
                                genes.append(genename)
                                reverse.append(True)
                     line = handle.readline()
                 else:
                     line = handle.readline()

             # Traverse file for every gene to store their sequences
             for x in range(len(genes)):
                 currentWD = os.getcwd()
                 handle = open(currentWD + "/" + infile, 'r')
                 line = handle.readline()
                 while True:
                     if not line:
                         break
                     if line.startswith(';') and not line.startswith(';;'):
                         geneName = line[1:].rstrip()
                         fullgenename = " ".join(geneName.split())
                         iterate = fullgenename.split(" ")
                         # """ Necessary informations are parsed """
                         genename = iterate[0]
                         if genename == genes[x]:
                             # """ Now we can store the sequence for it """
                             lines = []
                             line = handle.readline()
                             search = genename + " "
                             while True:
                                 if not line:
                                     break
                                 if search in line:
                                     break
                                 else:
                                     flag = False
                                     # If they are lowercase, this means they belong
                                     # to an intron which should not be taken into the sequence.
                                     for c in line:
                                         if c.islower():
                                             flag = True
                                             break
                                     if not flag:
                                         lines.append(line.rstrip())
                                         # Removes other genes starting or end.
                                         if line.startswith(';'):
                                             lines.remove(line.rstrip())
                                         line = handle.readline()
                                     else:
                                         line = handle.readline()

                             sequence = "".join(lines).replace(" ", "").replace("\r", "")
                             sequence = re.sub("\d","",sequence)
                             sequences.append(sequence)
                             line = handle.readline()

                         else:
                             line = handle.readline()
                     else:
                         line = handle.readline()

             # """ Now we should reverse 5'->3' strands to 3'->5' strands. """
             reversedsequences = [None]*(len(genes))
             for x in range(len(genes)):
               if reverse[x]:
                   reversedsequences[x] = reverse_sequence(sequences[x][::-1])
               else:
                   reversedsequences[x] = sequences[x]

             master_dict = {'species_name': speciesname , 'species_id': specie_id, 'genes_list': genes, 'sequences': reversedsequences}
             return master_dict


# reverse_sequence(): Reverses the "sequence" to its complementary strand. (5'->3'::3'->5')
def reverse_sequence(sequence):
    # """ Adenosine -> 1, then 1 -> Timine """
    tempseq = sequence.replace('A','1')
    # """ Timine -> 2, then 2 -> Adenosine """
    tempseq1 = tempseq.replace("T","2")
    # """ Cytosine -> 3, then 3 -> Guanine """
    tempseq2 = tempseq1.replace("C","3")
    # """ Guanine -> 4, then 4 -> Cytosine """
    tempseq3 = tempseq2.replace("G","4")
    # """ Now reverse back """
    tempseq4 = tempseq3.replace("1","T")
    tempseq5 = tempseq4.replace("2","A")
    tempseq6 = tempseq5.replace("3","G")
    tempseq7 = tempseq6.replace("4","C")
    reversed_sequence = tempseq7
    return reversed_sequence



# split_len(): Splits the sequence "seq" into "length" amount of seperate sequences.
def split_len(seq, length,position):
    return [seq[i:i+length] for i in range(position, len(seq), length)]



# getSubStrings(): Divides the given sequence "RNA" into 3-mers i.e [0-3,1-4,2-5 etc.]. Position should be 0,1 or 2.
def getSubStrings(RNA, position):
    return [RNA[i:i+3] for i in range(position, len(RNA), 3)]



# returnBlosumMatrix(): Returns a 2-dimensional dictionary of the Blosum62 Matrix.
def returnBlosumMatrix():
    Blosum62Matrix = {'A' : {'A':4, 'R':-1, 'N':-2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1,
     'L': -1, 'K': -1, 'M': -1,'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': 2,
     'Z': -1, 'X': 0, '*': -4},

      'R' : {'A':-1, 'R':5, 'N':0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3,
       'L': -2, 'K': 2, 'M': -1,'F': -3, 'P': -2, 'S': 1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1,
       'Z': 0, 'X': -1, '*': -4},

       'N' : {'A':-2, 'R':0, 'N':6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3,
        'L': -3, 'K': 0, 'M': -2,'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3,
        'Z': 0, 'X': -1, '*': -4},

        'D' : {'A':-2, 'R':-2, 'N':1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3,
         'L': -4, 'K': -1, 'M': -3,'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4,
         'Z': 1, 'X': -1, '*': -4},

         'C' : {'A':0, 'R':-3, 'N':-3, 'D': -3, 'C': 9, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': -1,
          'L': -1, 'K': -3, 'M': -1,'F': -2, 'P': -3, 'S': 1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3,
          'Z': -3, 'X': -2, '*': -4},

          'Q' : {'A':-1, 'R':1, 'N':0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3,
           'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0,
           'Z': 3, 'X': -1, '*': -4},

           'E' : {'A':1, 'R':0, 'N':0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3,
            'L': -3, 'K': 1, 'M': -2,'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1,
            'Z': 4, 'X': -1, '*': -4},

            'G' : {'A':0, 'R':-2, 'N':-0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4,
             'L': -4, 'K': -2, 'M': -3,'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1,
             'Z': -2, 'X': -1, '*': -4},

             'H' : {'A':-2, 'R':0, 'N':1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3,
              'L': -3, 'K': -1, 'M': -2,'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0,
              'Z': 0, 'X': -1, '*': -4},

              'I' : {'A':-1, 'R':-3, 'N':-3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4,
                'L': 2, 'K': -3, 'M': 1,'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3,
                'Z': -3, 'X': -1, '*': -4},

               'L' : {'A':-1, 'R':-2, 'N':-3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2,
                'L': 4, 'K': -2, 'M': 2,'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4,
                'Z': -3, 'X': -1, '*': -4},

                'K' : {'A':-1, 'R':2, 'N':0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3,
                 'L': -2, 'K': 5, 'M': -1,'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0,
                 'Z': 1, 'X': -1, '*': -4},

                 'M' : {'A':-1, 'R':-1, 'N':-2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1,
                  'L': 2, 'K': -1, 'M': 5,'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3,
                  'Z': -1, 'X': -1, '*': -4},

                  'F' : {'A':-2, 'R':-3, 'N':-3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0,
                   'L': 0, 'K': -3, 'M': 0,'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3,
                   'Z': -3, 'X': -1, '*': -4},

                   'P' : {'A':-1, 'R':-2, 'N':-2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3,
                    'L': -3, 'K': -1, 'M': -2,'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2,
                    'Z': -1, 'X': -2, '*': -4},

                    'S' : {'A':1, 'R':-1, 'N':1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2,
                     'L': -2, 'K': 0, 'M': -1,'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0,
                     'Z': 0, 'X': 0, '*': -4},

                     'T' : {'A':0, 'R':-1, 'N':0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1,
                      'L': -1, 'K': -1, 'M': -1,'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1,
                      'Z': -1, 'X': 0, '*': -4},

                      'W' : {'A':-3, 'R':-3, 'N':-4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3,
                       'L': -2, 'K': -3, 'M': -1,'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': -11, 'Y': 2, 'V': -3, 'B': -4,
                       'Z': -3, 'X': -2, '*': -4},

                       'Y' : {'A':-2, 'R':-2, 'N':-2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1,
                        'L': -1, 'K': -2, 'M': -1,'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3,
                        'Z': -2, 'X': -1, '*': -4},

                        'V' : {'A':0, 'R':-3, 'N':-3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': -3,
                         'L': 1, 'K': -2, 'M': 1,'F': -1, 'P': -2, 'S':-2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3,
                         'Z': -2, 'X': -1, '*': -4},

                         'B' : {'A':-2, 'R':-1, 'N':3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3,
                          'L': -4, 'K': 0, 'M': -3,'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4,
                          'Z': 1, 'X': -1, '*': -4},

                          'Z' : {'A':-1, 'R':0, 'N':0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3,
                           'L': -3, 'K': 1, 'M': -1,'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1,
                           'Z': 4, 'X': -1, '*': -4},

                           'X' : {'A':0, 'R':-1, 'N':-1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1,
                            'L': -1, 'K': -1, 'M': -1,'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1,
                            'Z': -1, 'X': -1, '*': -4},

                            '*' : {'A':-4, 'R':-4, 'N':-4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4,
                             'L': -4, 'K': -4, 'M': -4,'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4,
                             'Z': -4, 'X': -4, '*': 1}
      }
    return Blosum62Matrix

# function to get unique values
def unique(list1):

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
currentWD = os.getcwd()
try:
    os.mkdir(os.path.expanduser(currentWD + "/Output"))
except FileExistsError:
    pass

files_list = os.listdir(currentWD + "/MasterFiles")
for element in files_list:
    if "DS" in element:
        files_list.remove(element)

for i in range(len(files_list)):
    files_list[i] = "/MasterFiles/" + files_list[i]

instance = MasterFile(files_list)
