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

    def __init__(self, infile, kmer_length = 20, extension_length = 1, threshold = 0, gap_open_penalty = -2, match_score = 6, mismatch_penalty = -4, allow = 1):
        # Parameters:
            # Infile parameter should be a list of Master Files obtained from MFannot results.
                # TODO: Maybe we can include in our code MFannot command line tools so a person can input as fasta.
            # kmer_length: An integer which decides how many nucleotides/extension block should be.
                # Increasing kmer_length will make it harder to find a better extension.
            # extension_length: An integer decides that length(original sequence + extension) / length(reference_sequence)
            # threshold: An integer that decides the minimum score that local alignments of extensions can return.
            # gap_open_penalty, match_score, mismatch_penalty: Not used currently because of Blosum62 Matrix scores.
            # allow: An integer boolean that decides if we allow original sequence + extensions to be longer than reference sequence.

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
        self._alignments(list_of_dictionaries,list_of_extended_dictionaries, kmer_length, extension_length, threshold, gap_open_penalty, match_score, mismatch_penalty, allow, reference_species, common_genes)





    def _alignments(self, list_of_dictionaries,list_of_extended_dictionaries, kmer_length, extension_length, threshold, gap_open_penalty, match_score, mismatch_penalty, allow, reference_species, common_genes):
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
        scores_dictionary = [] # A dictionary that for each gene holds scores.
        for i in range(len(common_genes)):
            reference_sequence = ""
            reference_sequence_protein = ""
            scores = []
            start_end = []
            non_reference_species = []
            genes = []
            alignments = []
            original_sequences = []
            original_sequences_protein = []
            non_reference_species_id = []
            extensions = []
            for j in range(len(list_of_formatted_dictionaries)):
                # Finds the sequence for the reference specie:
                if list_of_formatted_dictionaries[j]['species_name'] == reference_species[i]:
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


            Blosum62Matrix = returnBlosumMatrix()
            # Traverses non-reference_sequences and do local alignment with reference, stores their scores.
            for j in range(len(original_sequences_protein)):
                submat = make_identity_substitution_matrix(4, -6, alphabet='ARNDCQEGHILKMFPSTWYVBZX*')
                    # TODO: Try to add scores for match/mismatch etc.
                original_string = Protein(str(original_sequences_protein[j]))
                reference_string = Protein(str(reference_sequence_protein))
                alignment, score, start_end_positions = local_pairwise_align_ssw(
                original_string,
                reference_string,
                substitution_matrix = submat
                )
                scores.append(score)
                start_end.append(start_end_positions)
                alignments.append(alignment)

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
                length_of_extension = len(extensions[j])
                length_of_original_sequence = len(original_sequences[j])
                length_of_possible_extension = int(length_of_original_sequence * extension_length)
                if length_of_possible_extension > length_of_extension:
                    trimmed_extensions.append(extensions[j])
                else:
                    trimmed_ex = ""
                    for k in range(length_of_possible_extension):
                        trimmed_ex = trimmed_ex + extensions[j][k]
                    trimmed_extensions.append(trimmed_ex)
            # Free extensions
            extensions.clear()
            # Depending on parameter allow - trims the extensions if necessary:
            allowed_extensions = []
            for j in range(len(trimmed_extensions)):
                if allow == 0:
                    length_of_reference_sequence = len(reference_sequence)
                    length_of_original_sequence = len(original_sequences[j])
                    length_of_extended = len(trimmed_extensions[j])
                    if length_of_extended + length_of_original_sequence > length_of_reference_sequence:
                        allowed_extension = ""
                        for k in range(length_of_reference_sequence - length_of_original_sequence):
                            allowed_extension = allowed_extension + trimmed_extensions[j][k]
                        replaced = allowed_extension.replace("!","")
                        allowed_extensions.append(replaced)
                    else:
                        replaced = trimmed_extensions[j].replace("!","")
                        allowed_extensions.append(replaced)
                else:
                    replaced = trimmed_extensions[j].replace("!","")
                    allowed_extensions.append(replaced)

            # Free trimmed extensions, not necessary anymore.
            trimmed_extensions.clear()
            # Transform extensions into proteins too.
            allowed_extensions_protein = []
            for j in range(len(allowed_extensions)):
                table_id = non_reference_species_id[j]
                gene = Seq(allowed_extensions[j],generic_dna)
                protein = gene.translate(table = table_id)
                allowed_extensions_protein.append(protein)

            score = {'gene_name': common_genes[i], 'original_species': non_reference_species, 'original_sequences': original_sequences, 'original_sequences_protein': original_sequences_protein, 'reference_specie': reference_species[i], 'reference_sequences': reference_sequence, 'reference_sequence_protein': reference_sequence_protein, 'scores': scores, 'start_end': start_end, 'alignments': alignments, 'extensions': allowed_extensions,
             'extensions_protein': allowed_extensions_protein}
            scores_dictionary.append(score)

        calculations = []
        # for each gene
        # (gene1, [specie1, specie2, ..], [max_score1,max_score2, ..], [extension1, extension2, ..])
        genes_already_printed = []
        for i in range(len(scores_dictionary)):
            if not scores_dictionary[i]['gene_name'] in genes_already_printed:
                # Create a gene file and move it to Output File.
                currentWD = os.getcwd()
                filename = scores_dictionary[i]['gene_name'] + ".txt"
                f= open(filename,"w+")
                os.rename(currentWD + "/" + filename, currentWD + "/Output/" + filename)
                f.write("Reference Specie: %s\n\n" % scores_dictionary[i]['reference_specie'])
                f.write("----- INITIAL LOCAL ALIGNMENT -----\n\n")

                species = []
                max_scores = []
                extended_alignments = []
                genes_already_printed.append(scores_dictionary[i]['gene_name'])
                # for each specie that is not reference
                for j in range(len(scores_dictionary[i]['original_species'])):
                    necessary_total = scores_dictionary[i]['original_sequences'][j] + scores_dictionary[i]['extensions'][j]
                    f.write("Specie: %s\n" % scores_dictionary[i]['original_species'][j])
                    species.append(scores_dictionary[i]['original_species'][j])
                    length_of_extension = len(scores_dictionary[i]['extensions'][j])
                    f.write("Protein Length Original: %s\n" % len(scores_dictionary[i]['original_sequences_protein'][j]))
                    f.write("Protein Length (Original + Extension): %s\n" % (length_of_extension + len(scores_dictionary[i]['original_sequences_protein'][j])))
                    f.write("Protein Length Reference: %s\n" % len(scores_dictionary[i]['reference_sequence_protein']))
                    f.write('Start/End L.A= %s\n' % scores_dictionary[i]['start_end'][j])
                    #f.write('Original Detected Stop Codon: %s' % detectedStop)
                    f.write('---\n\n')
                    # their extended sequence should be divided into k-mers
                    # If it has extensions, divide them and concatenate each of them with
                    # original sequence and then do local alignment, then store their scores.
                    scores = []
                    alignments = []
                    end_start = []

                    division_amount = 0
                    # Determine division
                    if length_of_extension < kmer_length:
                        division_amount = 1
                    else:
                        division_amount = int(length_of_extension/kmer_length)

                    # Find the place left in reference sequence
                    end_position_reference = scores_dictionary[i]['start_end'][0][1][1]
                    # Concatenate the last part of original sequence and the extensions for them
                    end_position_original = scores_dictionary[i]['start_end'][j][0][1]
                    concatenate_original = scores_dictionary[i]['original_sequences_protein'][j][end_position_original:]

                    reference_sequence_right = scores_dictionary[i]['reference_sequence_protein'][end_position_reference:]
                    divisible_extensions = concatenate_original + scores_dictionary[i]['extensions_protein'][j]
                    extensions = []
                    amount = int(length_of_extension/division_amount)
                    newString = []
                    if not length_of_extension == 0:
                        extensions = split_len(divisible_extensions, amount,0)
                        for k in range(len(extensions)):
                            newString.append(str(extensions[k]))
                            newStringProtein = Protein(str(''.join(newString)))
                            try:
                                refseq = Protein(str(reference_sequence_right))
                                try:
                                    alignment, score, start_end_positions = local_pairwise_align_ssw(
                                    newStringProtein,
                                    refseq,
                                    substitution_matrix = submat,
                                    )
                                    if score > threshold:
                                        # First remove the last extension which passed the threshold:
                                        newString.remove(str(extensions[k]))
                                        for i in range(len(str(extensions[k]))):
                                            # Now add proteins one by one and try to find the scores that are passing threshold:
                                            newString.append(str(extensions[k][i:i+1]))
                                            newStringProtein = Protein(str(''.join(newString)))
                                            alignment, score, start_end_positions = local_pairwise_align_ssw(
                                            newStringProtein,
                                            refseq,
                                            substitution_matrix = submat,
                                            )
                                            if score > threshold:
                                                # f.write("newStringProtein = %s\n" % str(newStringProtein))
                                                # f.write("LastAA = %s\n" % str(extensions[k][i:i+1]))
                                                # f.write("RefSeq = %s\n" % refseq)
                                                # f.write("Alignment = %s\n" % alignment)
                                                # f.write("Score = %s\n" % score)
                                except IndexError:
                                    pass
                            except ValueError:
                                pass
                    else:
                        pass

                    # divide_nmer = 5
                    # starts_lists = []
                    # before_maximum_lists = []
                    # if not maximum_score_index == -1:
                    #     # Divide the maximum score in divide_nmers starting from 0,1,2,.. to include frameshifts.
                    #     for i in range(divide_nmer):
                    #         before_maximum = ""
                    #         print(extensions[maximum_score_index])
                    #         starts_i = split_len(str(extensions[maximum_score_index]),divide_nmer, i)
                    #         print(starts_i)
                    #         starts_lists.append(starts_i)
                    #         for j in range(maximum_score_index):
                    #             before_maximum = before_maximum + extensions[j]
                    #         before_maximum = before_maximum + extensions[maximum_score_index][:i]
                    #         print(before_maximum)
                    #         before_maximum_lists.append(before_maximum)
                    #         newString = []
                    #         newString.append(str(before_maximum))
                    #     # Now for each of these divisions we will gather scores until we find the maximum one.
                    #     # For starts_0:
                    #     # Concatenate maximum extension with before extensions:
                    #     newString = []
                    #     newString.append(str(before_maximum))
                    #     scores_0, alignments_0, end_start_0 = [],[],[]
                    #     for k in range(len(starts_0)-1):
                    #         newString.append(str(starts_0[k]))
                    #         newStringProtein = Protein(str(''.join(newString)))
                    #         try:
                    #             refseq = Protein(str(reference_sequence_right))
                    #             try:
                    #                 alignment, score, start_end_positions = local_pairwise_align_ssw(
                    #                 newStringProtein,
                    #                 refseq,
                    #                 substitution_matrix = submat,
                    #                 )
                    #                 scores_0.append(score)
                    #                 alignments_0.append(alignment)
                    #                 end_start_0.append(start_end_positions)
                    #             except IndexError:
                    #                 scores_0.append(0)
                    #                 alignments_0.append("Nothing")
                    #                 end_start_0.append(0)
                    #         except ValueError:
                    #             scores_0.append(0)
                    #             alignments_0.append("Nothing")
                    #             end_start_0.append(0)
                    #     # Find best for scores_0
                    #     bestScore_0 = 0
                    #     bestIndex_0 = 0
                    #     for k in range(len(scores_0)):
                    #         if scores_0[k] > bestScore_0:
                    #             bestScore_0 = scores_0[k]
                    #             bestIndex_0 = k
                    #     # For starts_1:
                    #     newString = []
                    #     newString.append(str(before_maximum))
                    #     scores_1, alignments_1, end_start_1 = [],[],[]
                    #     for k in range(len(starts_1)-1):
                    #         newString.append(str(starts_1[k]))
                    #         newStringProtein = Protein(str(''.join(newString)))
                    #         try:
                    #             refseq = Protein(str(reference_sequence_right))
                    #             try:
                    #                 alignment, score, start_end_positions = local_pairwise_align_ssw(
                    #                 newStringProtein,
                    #                 refseq,
                    #                 substitution_matrix = submat,
                    #                 )
                    #                 scores_1.append(score)
                    #                 alignments_1.append(alignment)
                    #                 end_start_1.append(start_end_positions)
                    #             except IndexError:
                    #                 scores_1.append(0)
                    #                 alignments_1.append("Nothing")
                    #                 end_start_1.append(0)
                    #         except ValueError:
                    #             scores_1.append(0)
                    #             alignments_1.append("Nothing")
                    #             end_start_1.append(0)
                    #     # Find best for scores_0
                    #     bestScore_1 = 0
                    #     bestIndex_1 = 0
                    #     for k in range(len(scores_1)):
                    #         if scores_1[k] > bestScore_1:
                    #             bestScore_1 = scores_1[k]
                    #             bestIndex_1 = k
                    #
                    #     # For starts_2:
                    #     newString = []
                    #     newString.append(str(before_maximum))
                    #     scores_2, alignments_2, end_start_2 = [],[],[]
                    #     for k in range(len(starts_2)-1):
                    #         newString.append(str(starts_2[k]))
                    #         newStringProtein = Protein(str(''.join(newString)))
                    #         try:
                    #             refseq = Protein(str(reference_sequence_right))
                    #             try:
                    #                 alignment, score, start_end_positions = local_pairwise_align_ssw(
                    #                 newStringProtein,
                    #                 refseq,
                    #                 substitution_matrix = submat,
                    #                 )
                    #                 scores_2.append(score)
                    #                 alignments_2.append(alignment)
                    #                 end_start_2.append(start_end_positions)
                    #             except IndexError:
                    #                 scores_2.append(0)
                    #                 alignments_2.append("Nothing")
                    #                 end_start_2.append(0)
                    #         except ValueError:
                    #             scores_2.append(0)
                    #             alignments_2.append("Nothing")
                    #             end_start_2.append(0)
                    #     # Find best for scores_0
                    #     bestScore_2 = 0
                    #     bestIndex_2 = 0
                    #     for k in range(len(scores_2)):
                    #         if scores_2[k] > bestScore_2:
                    #             bestScore_2 = scores_2[k]
                    #             bestIndex_2 = k
                    #     # Find best frameshift and best index of the alignment.
                    #     best_total_shift = -1
                    #     best_total_index = -1
                    #     best_total_score = -1
                    #     if bestScore_0 >= bestScore_1:
                    #         if bestScore_0 >= bestScore_2:
                    #             best_total_shift = 0
                    #             best_total_index = bestIndex_0
                    #             best_total_score = bestScore_0
                    #         else:
                    #             best_total_shift = 2
                    #             best_total_index = bestIndex_2
                    #             best_total_score = bestScore_2
                    #     else:
                    #         if bestScore_1 >= bestScore_2:
                    #             best_total_shift = 1
                    #             best_total_index = bestIndex_1
                    #             best_total_score = bestScore_1
                    #         else:
                    #             best_total_shift = 2
                    #             best_total_index = bestIndex_2
                    #             best_total_score = bestScore_2
                    #     # Determine possible stop codon. Possible stop codon is the best_total_index of the best_total_shift
                    #     possibleStopCodon = ""
                    #
                    #     left_extension = end_start[maximum_score_index][0][0]
                    #     right_extension = end_start[maximum_score_index][0][1]
                    #     left_reference = end_start[maximum_score_index][1][0]
                    #     right_reference = end_start[maximum_score_index][1][1]
                    #     if best_total_shift == 0:
                    #         possibleStopCodon = starts_0[best_total_index]
                    #     elif best_total_shift == 1:
                    #         possibleStopCodon = starts_1[best_total_index]
                    #     else:
                    #         possibleStopCodon = starts_2[best_total_index]
                    #
                    #     try1_original_left = (end_position_original + left_extension) * 3
                    #     try1_original_right = (try1_original_left + 9)
                    #     try2 = necessary_total[try1_original_left-12:try1_original_right+12]
                    #     gene = Seq(try2,generic_dna)
                    #     protein = gene.translate(table = 1)
                    #     # Now search possible Stop Codon in protein and determine start-end locations correctly.
                    #     # Maybe include 1 more aminoacid to the right while taking best 3-mer.
                    #     # => 4 or 3 aminoacids for each gene, print them as possible stop codons.
                    #     threeMer,threeMer1,threeMer2 = "a","b","c"
                    #     for i in range(len(str(protein))-2):
                    #         if str(protein)[i:i+3] == possibleStopCodon:
                    #             startLoc = i
                    #             threeMer = try2[startLoc+3*i-3:startLoc+3*i+6]
                    #             threeMer1 = try2[startLoc+3*i-2:startLoc+3*i+7]
                    #             threeMer2 = try2[startLoc+3*i-4:startLoc+3*i+5]
                    #             break
                    #     # Considering Frameshifts.
                    #     gene = Seq(threeMer,generic_dna)
                    #     protein = gene.translate(table = 1)
                    #     gene1 = Seq(threeMer1,generic_dna)
                    #     protein1 = gene1.translate(table = 1)
                    #     gene2 = Seq(threeMer2,generic_dna)
                    #     protein2 = gene2.translate(table = 1)
                    #
                    #     if str(protein) == possibleStopCodon:
                    #         f.write("--List--\n")
                    #         f.write("Protein: %s\n" % protein)
                    #         possibleStopCodonList = [threeMer[i:i+3] for i in range(0, len(threeMer), 3)]
                    #         for i in range(len(possibleStopCodonList)):
                    #             f.write("%s\n" % possibleStopCodonList[i])
                    #         f.write("\n")
                    #     elif str(protein1) == possibleStopCodon:
                    #         f.write("--List--\n")
                    #         f.write("Protein: %s\n" % protein1)
                    #         possibleStopCodonList = [threeMer1[i:i+3] for i in range(0, len(threeMer1), 3)]
                    #         for i in range(len(possibleStopCodonList)):
                    #             f.write("%s\n" % possibleStopCodonList[i])
                    #         f.write("\n")
                    #     else:
                    #         f.write("--List--\n")
                    #         f.write("Protein: %s\n" % protein2)
                    #         possibleStopCodonList = [threeMer2[i:i+3] for i in range(0, len(threeMer2), 3)]
                    #         for i in range(len(possibleStopCodonList)):
                    #             f.write("%s\n" % possibleStopCodonList[i])
                    #         f.write("\n")
                        # else:
                        #     f.write("None found. Probably the last codon in sequence is correct stop codon.\n")
                        #     f.write("\n")



                                # lastCodon = try2[-3:] # Now we have the last codon.
                                # # Complete sequence to look for
                                # checkBefore = try1[:3*(right_extension + end_position_original + x*amount + 1)]
                                # # Complete sequence to look for in a reversed format
                                # checkBeforeReverse = checkBefore[::-1]
                                # # Divide to 3-mers from beginning.
                                # mers3 = getSubStrings(checkBeforeReverse, 0)
                                # # Search in the sequence if there is a same codon. Frameshifts are ignored.
                                # found = False
                                # for search in range(1,len(mers3)):
                                #     if mers3[search] == lastCodon[::-1]:
                                #         found = True
                                #         break
                                #
                                # if not found:
                                #     f.write("\n")
                                #     f.write("LOCAL ALIGNMENTS WITH EXTENSIONS:\n")
                                #     f.write("-------------------------------\n")
                                #     f.write("Specie: %s\n" % scores_dictionary[i]['original_species'][j])
                                #     f.write("Reference Distance From End: %s\n" % distanceReference)
                                #     f.write("Original Distance From End: %s\n" % distanceOriginal)
                                #     f.write("Score: %.2f\n" % newScore)
                                #     f.write("Start/End E.: [(%s, %s), (%s, %s)]\n" % (left_extension + end_position_original + x*amount, right_extension + end_position_original + x*amount + 1, left_reference + end_position_reference, right_reference + end_position_reference))
                                #     tempOriginal = scores_dictionary[i]['original_sequences_protein'][j] + scores_dictionary[i]['extensions_protein'][j]
                                #     #Original S.P: %s" % tempOriginal[end_position_original + amount*x + left_extension: end_position_original + right_extension + amount*x + 1])
                                #     #Reference S.P: %s" % scores_dictionary[i]['reference_sequence_protein'][left_reference + end_position_reference:right_reference + end_position_reference + 1])
                                #     f.write("Alignment: %s\n" % alignments[x])
                                #     f.write("POSSIBLE REASSIGNED STOP CODON: %s\n" % lastCodon)
                                #     f.write("\n")
                                #
                                # found = False

            # # For each possible gene and their reference sequences
            # #print("REFERENCE SEQUENCE REASSIGNMENT RESULT:")
            # #print("-------------------------------")
            # for i in range(len(scores_dictionary)):
            #     #print("Gene name: %s" % scores_dictionary[i]['gene_name'])
            #     #print("Reference specie: %s" % scores_dictionary[i]['reference_specie'])
            #     reference_sequence = scores_dictionary[i]['reference_sequences']
            #     reference_sequence_protein = scores_dictionary[i]['reference_sequence_protein']
            #     # Divide reference sequence protein into extensions.
            #     divisions_protein = split_len(reference_sequence_protein, amount)
            #     # For each non-reference species:
            #     for j in range(len(scores_dictionary[i]['original_sequences_protein'])):
            #         # Store score of the initial local alignment.
            #         initial_score = scores_dictionary[i]['scores'][j]
            #         initial_start_end = scores_dictionary[i]['start_end'][j]
            #         # Now go backwards in divisions_protein and check if score is better than original.
            #         count = 0
            #         for k in range(len(divisions_protein),-1,-1):
            #             new_reference_protein = new_reference_protein + divisions_protein[count]
            #


    def _internal_extendedparser(self, infile, list_of_dictionaries):
         # For each gene that is common, we need to find if the gene is extendable or not.
         #     1- We need to find which genes are common in all species, and store them.
         #     2- Filter the interesting genes
         # List of interesting genes: nad / atp / cob / cox / sdh / tat can be extended. Other genes are omitted.
        common_genes = []
        for x in range(len(list_of_dictionaries[0]['genes_list'])):
            for y in range(1,len(list_of_dictionaries)):
                for z in range(len(list_of_dictionaries[y]['genes_list'])):
                    if list_of_dictionaries[0]['genes_list'][x] == list_of_dictionaries[y]['genes_list'][z]:
                        if "G-nad" in list_of_dictionaries[0]['genes_list'][x] or "G-atp" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes.append(list_of_dictionaries[0]['genes_list'][x])
                        if "G-cob" in list_of_dictionaries[0]['genes_list'][x] or "G-cox" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes.append(list_of_dictionaries[0]['genes_list'][x])
                        if "G-sdh" in list_of_dictionaries[0]['genes_list'][x] or "G-tat" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes.append(list_of_dictionaries[0]['genes_list'][x])

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
        # for x in range(len(chunks)):
        #     maximum = max(chunks[x])
        #     for y in range(len(chunks[x])):
        #         if maximum == chunks[x][y]:
        #             reference_species.append(list_of_dictionaries[y]['species_name'])
        #             break
        for y in range(len(common_genes)):
            for x in range(len(list_of_dictionaries)):
                if "minor" in list_of_dictionaries[x]['species_name']:
                    reference_species.append(list_of_dictionaries[x]['species_name'])
                    break

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
