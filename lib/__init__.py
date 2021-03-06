from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import *
from skbio.alignment import local_pairwise_align_ssw, make_identity_substitution_matrix, global_pairwise_align_protein
from skbio import TabularMSA, DNA, Protein
from .Parser import MasterFile, MasterCollection
from .utils import split_len,_one_hot_encode,_create_one_hot_array,_stop_mapper
import os
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
from numpy import array
from numpy import argmax
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt

__all__ = ['MasterFile', 'StopChecker', 'MasterCollection']

class StopChecker:
    def __init__(self, mcol, outdir, kmer_length=15, threshold=10, match_score=10, mismatch_penalty=-4, submat=None, **kwargs):
        """Take a list of  collection of MasterFile as input and check stop codon"""
        self.mcol = mcol
        self.alignments(mcol, kmer_length, threshold, match_score, mismatch_penalty)
        self._give_meaning()
        self._learn()

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
        e_values_total_reference = []
        e_values_total_original = []
        pident_total_reference = []
        pident_total_original = []
        self.scores_list = []
        # List includes a dictionary for each gene
        for i in range(len(common_genes)):
            reference_specie_name = ""
            reference_sequence = ""
            reference_sequence_protein = ""
            non_reference_species = []
            original_species_names = []
            genes = []
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
                            # # Debug: -1- Passed.
                            # print("Formatted: %s" % formatted)
                            # print("Protein: %s" % protein)
                            original_sequences_protein.append(protein)
                            break
            # print("Common genes[i] = %s" % common_genes[i])
            # print("Non-reference-species: %s" % non_reference_species)

            # pairs_name: Stores the pairs in a string format-> Bracteacoccus_minor / Bracteacoccus_aerius i.e
            pairs_name = []
            # pairs_result: is a list of lists. for each pairs_name, it stores score / start_end / alignments.
            pairs_results = []
            # start_end: stores start_end results
            start_end = []
            # scores: stores bitscores
            scores = []
            # pident: stores how much identical are two sequences.
            pident = []
            # e-values: stores initial e-values.
            e_values_initials = []
            # Traverses non-reference_sequences and do local alignment with reference, stores their scores.
            for j in range(len(non_reference_species)):
                # results: a list that stores scores, start_end and alignments.
                # results = []
                # ! Put pairs in correct order with start-end order.
                pairs_name.append("%s/%s" % (non_reference_species[j],reference_specie_name))
                original_string = SeqRecord(original_sequences_protein[j],
                           id="seq1")
                reference_string = SeqRecord(reference_sequence_protein,
                           id="seq2")
                SeqIO.write(original_string, "seq1.fasta", "fasta")
                SeqIO.write(reference_string, "seq2.fasta", "fasta")

                # Run BLAST and parse the output
                command = './blastp -outfmt 6 -query seq1.fasta -subject seq2.fasta'
                p = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE)
                result = p.stdout.readlines()[0].decode('utf-8')
                la_data = result.split('\t')
                # # Debug: -2- Passed.
                # print("Local Alignment Data: %s" % la_data)

                # 1.	 qseqid	 query (e.g., gene) sequence id
                # 2.	 sseqid	 subject (e.g., reference genome) sequence id
                # 3.	 pident	 percentage of identical matches
                # 4.	 length	 alignment length
                # 5.	 mismatch	 number of mismatches
                # 6.	 gapopen	 number of gap openings
                # 7.	 qstart	 start of alignment in query
                # 8.	 qend	 end of alignment in query
                # 9.	 sstart	 start of alignment in subject
                # 10.	 send	 end of alignment in subject
                # 11.	 evalue	 expect value
                # 12.	 bitscore	 bit score
                start_end_original = []
                start_end_original.append(int(la_data[6]))
                start_end_original.append(int(la_data[7]))
                start_end_reference = []
                start_end_reference.append(int(la_data[8]))
                start_end_reference.append(int(la_data[9]))
                total = []
                total.append(start_end_original)
                total.append(start_end_reference)
                pairs_results.append(total)
                start_end.append(total)
                # Stored bitscore, percentage identicality and e-values for initials.
                # Between reference and non-reference pairs.
                scores.append(float(la_data[11]))
                pident.append(float(la_data[2]))
                e_values_initials.append(float(la_data[10]))


            # Now do pairwise between all original sequences:
            for j in range(len(non_reference_species)-1):
                for k in range(j+1,len(non_reference_species)):
                    pairs_name.append("%s/%s" % (non_reference_species[j],non_reference_species[k]))
                    j_sequence = SeqRecord(original_sequences_protein[j],
                               id="seq1")
                    k_sequence = SeqRecord(original_sequences_protein[k],
                               id="seq2")
                    SeqIO.write(j_sequence, "seq1.fasta", "fasta")
                    SeqIO.write(k_sequence, "seq2.fasta", "fasta")

                    # Run BLAST and parse the output
                    command = './blastp -outfmt 6 -query seq1.fasta -subject seq2.fasta'
                    p = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE)
                    result = p.stdout.readlines()[0].decode('utf-8')
                    la_data = result.split('\t')
                    # # Debug: -3- Passed.
                    # print("Local Alignment Data: %s" % la_data)
                    start_end_original = []
                    start_end_original.append(int(la_data[6]))
                    start_end_original.append(int(la_data[7]))
                    start_end_reference = []
                    start_end_reference.append(int(la_data[8]))
                    start_end_reference.append(int(la_data[9]))
                    total = []
                    total.append(start_end_original)
                    total.append(start_end_reference)

                    pairs_results.append(total)
                    start_end.append(total)
                    # Stored bitscore, percentage identicality and also e-values for initials.
                    # Between non-reference and non-reference pairs.
                    scores.append(float(la_data[11]))
                    pident.append(float(la_data[2]))
                    e_values_initials.append(float(la_data[10]))

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
             'extensions_protein': allowed_extensions_protein, 'start_end': start_end, 'scores': scores, 'pident': pident, 'evalues':e_values_initials}
            self.scores_list.append(score_dictionary)

        # Here we will start a new, final chapter everyone!
        # First compare all sequences with reference sequence.
        for i in range(len(self.scores_list)):
            e_values_extensions = []
            pairs_results_extensions = []
            pairs_results_extensions_scores = []
            pident_extensions = []
            # Pairwise reference with each non-reference specie.
            for j in range(len(self.scores_list[i]['original_species'])):
                # Find the place left in reference sequence.
                end_position_reference = self.scores_list[i]['start_end'][j][1][1]
                # Find end position of original sequence.
                end_position_original = self.scores_list[i]['start_end'][j][0][1]
                # Take whats left of reference from initial local alignment.
                reference_sequence_right = self.scores_list[i]['reference_sequence_protein'][end_position_reference:]
                # Concatenate the last part of original sequence and the extensions for them.
                concatenated_original_seq = self.scores_list[i]['original_sequences_protein'][j][end_position_original:]
                concatenated = concatenated_original_seq + self.scores_list[i]['extensions_protein'][j]
                # # Debug: -4- Passed.
                # print("End p. R: %d, End p. O: %d" % (end_position_reference,end_position_original))
                # print("Length of Reference Sequence: %d, L.o.Original.s: %d" % (len(self.scores_list[i]['reference_sequence_protein']),len(self.scores_list[i]['original_sequences_protein'][j])))
                # print("Concatenated: %s" % concatenated)
                # print("Ref.: %s" % reference_sequence_right)
                newStringProtein = SeqRecord(concatenated,
                           id="seq1")
                refseq = SeqRecord(reference_sequence_right,
                           id="seq2")
                SeqIO.write(newStringProtein, "seq1.fasta", "fasta")
                SeqIO.write(refseq, "seq2.fasta", "fasta")
                # print("Org: %s" % str(newStringProtein))
                # print("Ref: %s" % str(refseq))
                # Run BLAST and parse the output
                command = './blastp -outfmt 6 -query seq1.fasta -subject seq2.fasta'
                p = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE)
                try:
                    result = p.stdout.readlines()[0].decode('utf-8')
                    la_data = result.split('\t')
                except IndexError:
                    la_data = []
                    for any in range(12):
                        la_data.append(0)
                # # Debug: -5- : Always 0,0,0,0,0..
                # print("Local Alignment Data: %s" % la_data)
                start_end_original = []
                start_end_original.append(int(la_data[6]))
                start_end_original.append(int(la_data[7]))
                start_end_reference = []
                start_end_reference.append(int(la_data[8]))
                start_end_reference.append(int(la_data[9]))
                e_values_extensions.append(float(la_data[10]))
                pident_extensions.append(float(la_data[2]))
                total = []
                total.append(start_end_original)
                total.append(start_end_reference)

                # Calculate start end position with respect to beginning.
                # start_location = start of the extension + other extensions before * their length + end of local alignment

                start_location = total[0][0] + self.scores_list[i]['start_end'][j][0][0]
                # start_location_r = end location of initial alignment + start position
                start_location_r = total[1][0] + self.scores_list[i]['start_end'][j][1][0]
                # end_location = end of the extension + other extensions before * their length + end of local alignment
                end_location = total[0][1] + self.scores_list[i]['start_end'][j][0][1]
                end_location_r = total[1][1] + self.scores_list[i]['start_end'][j][1][1]
                start_end_position_reference = []
                start_end_position_reference.append(start_location_r)
                start_end_position_reference.append(end_location_r)
                start_end_position_original = []
                start_end_position_original.append(start_location)
                start_end_position_original.append(end_location)

                start_end_position = []
                start_end_position.append(start_end_position_reference)
                start_end_position.append(start_end_position_original)
                # # Debug: -6- Passed.
                # print("Start Ref: %d, End Ref: %d" % (start_location_r,end_location_r))
                # print("Start Org: %d, End Org: %d" % (start_location,end_location))
                # print("Bitscore = Score: %d" % la_data[11])
                pairs_results_extensions.append(start_end_position)
                pairs_results_extensions_scores.append(float(la_data[11]))
            data3_temporary.append(pairs_results_extensions)
            data5_temporary.append(pairs_results_extensions_scores)
            e_values_total_reference.append(e_values_extensions)
            pident_total_reference.append(pident_extensions)

        # pairwise.
        for i in range(len(self.scores_list)):
            e_values_extensions = []
            pident_extensions = []
            pairs_results_extensions_j = []
            pairs_results_extensions_j_score = []
            for j in range(len(self.scores_list[i]['original_species'])-1):
                for k in range(j+1,len(self.scores_list[i]['original_species'])):
                    # Find the place left in reference sequence
                    orig_spec = len(self.scores_list[i]['original_species'])
                    j_loc = j + orig_spec

                    end_position_j = self.scores_list[i]['start_end'][j_loc][0][1]
                    end_position_k = self.scores_list[i]['start_end'][j_loc][1][1]
                    # Take whats left of reference from initial local alignment.
                    j_sequence_right = self.scores_list[i]['original_sequences_protein'][j][end_position_j:]
                    k_sequence_right = self.scores_list[i]['original_sequences_protein'][k][end_position_k:]
                    j_sequence_total = j_sequence_right + self.scores_list[i]['extensions_protein'][j]
                    k_sequence_total = k_sequence_right + self.scores_list[i]['extensions_protein'][k]

                    j_sequence = SeqRecord(j_sequence_total,
                               id="seq1")
                    k_sequence = SeqRecord(k_sequence_total,
                               id="seq2")
                    SeqIO.write(j_sequence, "seq1.fasta", "fasta")
                    SeqIO.write(k_sequence, "seq2.fasta", "fasta")

                    # Run BLAST and parse the output
                    command = './blastp -outfmt 6 -query seq1.fasta -subject seq2.fasta'
                    p = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE)
                    try:
                        result = p.stdout.readlines()[0].decode('utf-8')
                        la_data = result.split('\t')
                    except IndexError:
                        la_data = []
                        for any in range(12):
                            la_data.append(0)

                    start_end_original = []
                    start_end_original.append(int(la_data[6]))
                    start_end_original.append(int(la_data[7]))
                    start_end_reference = []
                    start_end_reference.append(int(la_data[8]))
                    start_end_reference.append(int(la_data[9]))
                    total = []
                    total.append(start_end_original)
                    total.append(start_end_reference)

                    # start_location = start of the extension + other extensions before * their length + end of local alignment
                    start_location = total[0][0] + self.scores_list[i]['start_end'][j][0][1]
                    # start_location_r = end location of initial alignment + start position
                    start_location_r = total[1][0] + self.scores_list[i]['start_end'][j][1][1]
                    # end_location = end of the extension + other extensions before * their length + end of local alignment
                    end_location = total[0][1] + self.scores_list[i]['start_end'][j][0][1]
                    end_location_r = total[1][1] + self.scores_list[i]['start_end'][j][1][1]
                    start_end_position_reference = []
                    start_end_position_reference.append(start_location_r)
                    start_end_position_reference.append(end_location_r)
                    start_end_position_original = []
                    start_end_position_original.append(start_location)
                    start_end_position_original.append(end_location)

                    start_end_position = []
                    start_end_position.append(start_end_position_reference)
                    start_end_position.append(start_end_position_original)

                    pairs_results_extensions_j_score.append(float(la_data[11]))
                    pairs_results_extensions_j.append(start_end_position)
                    e_values_extensions.append(float(la_data[10]))
                    pident_extensions.append(float(la_data[2]))
                    # # Debug: -7- Passed.
                    # print("Start/end position reference: %s" % start_end_position_reference)
                    # print("Start/end position original: %s" % start_end_position_original)
                    # print("Bitscore: %s" % la_data[11])

            data4_temporary.append(pairs_results_extensions_j)
            e_values_total_original.append(e_values_extensions)
            pident_total_original.append(pident_extensions)
            data6_temporary.append(pairs_results_extensions_j_score)

        for i in range(len(data6_temporary)):
            for k in range(len(data4_temporary[i])):
                data3_temporary[i].append(data4_temporary[i][k])

        # Storing variable informations for modelling DT/RF.
        mean_length_of_extensions = [] # Stores mean length of extensions
        length_of_genes = [] # Stores mean length of genes.
        frequency_of_stop_codon = [] # Check evolutionary close species' frequencies of that stop codon
        stop_codons_evolutionary = [] # Useless
        frequency_evolutionary = []  # Frequency of the stop codon in all evolutionary close stop codons
        mean_similarity_extension = [] # bitscores obtained from blast for extensions.
        mean_e_values_extensions = [] # e-values obtained from blast for extensions.
        mean_similarity_initials = [] # Similar to mean similarity of extensions but for initials.
        mean_e_values_initials = [] # Similar to mean_e_values_extensions but for initials.
            # 3.	 pident	 percentage of identical matches
            # 5.	 mismatch	 number of mismatches
            # 6.	 gapopen	 number of gap openings
        mean_identical_match_percentage_initials = [] # Stores percentage identical values.
        mean_identical_match_percentage_extensions = [] # Similar to mean_identical_match_percentage_initials but for extensions.

        for i in range(len(self.scores_list)): # For each gene // Genes should be taken apart seperately while collecting data
            # The length of the extension
            # For the reference specie:
            # We should do the same order when we were calculating scores/start.end positions.
            reference_extension_total = 0
            for j in range(len(self.scores_list[i]['original_species'])):
                reference_extension_total += int(data3_temporary[i][j][0][1]) - int(data2_temporary[i][j][0][1])
            # print(self.scores_list[i]['original_species'])
            # print("")
            mean_length_of_extensions.append(reference_extension_total/(len(self.scores_list[i]['original_species'])))
            # We need a better version of data structure, too much calculation involved while trying to find all non-references.
            # The data structure should be storing all information for all genes for each specie seperately.

            for j in range(len(self.scores_list[i]['original_species'])):
                mloe_dataset = []
                mloe_dataset_in = []
                # in the order of the original species
                for k in range(len(data1_temporary[i])):
                    data1_split = data1_temporary[i][k].split('/')
                    if self.scores_list[i]['original_species'][j] in data1_split[0]: # Left side
                        mloe_dataset.append(data3_temporary[i][k][0][1]) # Left end
                        mloe_dataset_in.append(data2_temporary[i][k][0][1])
                    if self.scores_list[i]['original_species'][j] in data1_split[1]: # Right side
                        mloe_dataset.append(data3_temporary[i][k][1][1]) # Right end
                        mloe_dataset_in.append(data2_temporary[i][k][1][1])
                mloe = 0
                # Now that we have the dataset for each gene, we can keep calculating.
                for k in range(len(mloe_dataset)):
                    mloe += float(mloe_dataset[k]) - float(mloe_dataset_in[k])
                mean_length_of_extensions.append(mloe/(len(self.scores_list[i]['original_species'])))

            # The length of the gene compared to other homologs:
            length_of_genes.append(len(self.scores_list[i]['reference_sequences']))
            for j in range(len(self.scores_list[i]['original_species'])): # For each non-reference-specie
                length_of_genes.append(len(self.scores_list[i]['original_sequences'][j]))
            # Frequency of the putative stop codon in all coding region
                # Divide into 3-mers, count how many stop codons occur in the all coding region (ignoring frameshift)
            substrings = split_len(self.scores_list[i]['reference_sequences'][::-1],3,0) # Reverse, divide 3-mers
            for j in range(len(substrings)):
                substrings[j] = substrings[j][::-1] # Reverse again each element. First element is the stop codon.
            countStopCodons = 0
            # Stop codon is first element, add to this specie's stop codon:
            stop_codons_evolutionary.append(substrings[0])
            for j in range(len(substrings)):
                if substrings[j] == substrings[0]:
                    countStopCodons += 1
            frequency_of_stop_codon.append(countStopCodons/len(substrings))
            for j in range(len(self.scores_list[i]['original_species'])):
                substrings = split_len(self.scores_list[i]['original_sequences'][j][::-1],3,0) # Reverse, divide 3-mers
                for k in range(len(substrings)):
                    substrings[k] = substrings[k][::-1] # Reverse again each element. First element is the stop codon.
                countStopCodons = 0
                # Add into list of stop codons:
                stop_codons_evolutionary.append(substrings[0])
                for k in range(len(substrings)):
                    if substrings[k] == substrings[0]:
                        countStopCodons += 1
                frequency_of_stop_codon.append(countStopCodons/len(substrings))

        counter = []
        onlyOnce = True # Store counter only once not genes times.
        for i in range(len(self.scores_list)):
            if onlyOnce:
                for j in range(len(self.scores_list[i]['original_species'])-1):
                    for k in range(j+1,len(self.scores_list[i]['original_species'])):
                        counter.append("(%s,%s)" % (j,k))
                onlyOnce = False

            rev = self.scores_list[i]['reference_sequences'][::-1]
            stop_rev = rev[:3]
            stop = stop_rev[::-1]
            stopCodonCount = 0
            for j in range(len(stop_codons_evolutionary)):
                if stop_codons_evolutionary[j] == stop:
                    stopCodonCount += 1
            frequency_evolutionary.append(stopCodonCount/len(stop_codons_evolutionary))
            for j in range(len(self.scores_list[i]['original_species'])):
                stopCodonCount = 0
                rev = self.scores_list[i]['original_sequences'][j][::-1]
                stop_rev = rev[:3]
                stop = stop_rev[::-1]
                for k in range(len(stop_codons_evolutionary)):
                    if stop_codons_evolutionary[k] == stop:
                        stopCodonCount += 1
                frequency_evolutionary.append(stopCodonCount/len(stop_codons_evolutionary))
        # Mean similarity of the extension to other genomes
        # How to determine mean similarity of extensions with other genomes?
            meansim_ref = 0
            for j in range(len(data5_temporary[i])):
                meansim_ref += data5_temporary[i][j]
            meansim_ref /= len(self.scores_list[i]['original_species'])
            mean_similarity_extension.append(meansim_ref) # Reference species similarities are stored.

            # Test
            # print("Length of data6_temporary: %s" % len(data6_temporary))
            # print("Length of data6_temporary[0]/[1]/[2]: %s/%s/%s" % (len(data6_temporary[0]),len(data6_temporary[1]),len(data6_temporary[2])))
            # print("")
            # print("Length of data5_temporary: %s" % len(data5_temporary))
            # print("Length of data5_temporary[0]/[1]/[2]: %s/%s/%s" % (len(data5_temporary[0]),len(data5_temporary[1]),len(data5_temporary[2])))

            # Implement mean similarity extension for non-reference sequences.
            for j in range(len(self.scores_list[i]['original_species'])):
                meansim_org = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        meansim_org += float(data6_temporary[i][j])
                meansim_org += data5_temporary[i][j]
                meansim_org /= len(self.scores_list[i]['original_species'])
                mean_similarity_extension.append(meansim_org)

            meansim_org = 0
            # For reference specie:
            # 'scores': scores
            for j in range(len(self.scores_list[i]['original_species'])): # 4
                meansim_org += self.scores_list[i]['scores'][j]
            meansim_avg = meansim_org / len(self.scores_list[i]['original_species'])
            mean_similarity_initials.append(meansim_avg)

            # For non-reference species:
            # Should start from 4th element, up to last element. Again counter can help us.
            # First organise scores:
            scores_organised = []
            for j in range(len(self.scores_list[i]['original_species']),len(self.scores_list[i]['scores'])):
                scores_organised.append(self.scores_list[i]['scores'][j])
            # Now counter method:
            for j in range(len(self.scores_list[i]['original_species'])):
                meansim_org = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        meansim_org += scores_organised[j]
                meansim_org += self.scores_list[i]['scores'][j]
                meansim_org /= len(self.scores_list[i]['original_species'])
                mean_similarity_initials.append(meansim_org)

            # e-values for reference:
            e_value_total = 0
            for j in range(len(self.scores_list[i]['original_species'])):
                e_value_total += e_values_total_reference[i][j]
            e_values_average = e_value_total / len(self.scores_list[i]['original_species'])
            mean_e_values_extensions.append(e_values_average)

            # e_values for non-reference species:
            for j in range(len(self.scores_list[i]['original_species'])):
                e_value_total = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        e_value_total += e_values_total_original[i][j]
                e_value_total += e_values_total_reference[i][j]
                e_values_average = e_value_total / len(self.scores_list[i]['original_species'])
                mean_e_values_extensions.append(e_values_average)

            # e-values initials for reference:
            # 'evalues':e_values_initials
            e_value_total = 0
            for j in range(len(self.scores_list[i]['original_species'])):
                e_value_total += self.scores_list[i]['evalues'][j]
            e_values_average = e_value_total / len(self.scores_list[i]['original_species'])
            mean_e_values_initials.append(e_values_average)

            # e_values initials for non-reference species:
            # First again organise:
            evalues_organised = []
            for j in range(len(self.scores_list[i]['original_species']),len(self.scores_list[i]['evalues'])):
                evalues_organised.append(self.scores_list[i]['evalues'][j])
            # Now counter method:
            for j in range(len(self.scores_list[i]['original_species'])):
                e_values_total = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        e_values_total += evalues_organised[j]
                e_values_total += self.scores_list[i]['evalues'][j]
                e_values_total /= len(self.scores_list[i]['original_species'])
                mean_e_values_initials.append(e_values_total)

            # Percentage identicals for initials, reference specie:
            # 'pident': pident
            pident_total = 0
            for j in range(len(self.scores_list[i]['original_species'])):
                pident_total += self.scores_list[i]['pident'][j]
            pident_average = pident_total / len(self.scores_list[i]['original_species'])
            mean_identical_match_percentage_initials.append(pident_average)
            # Organise pident values:
            pident_organised = []
            for j in range(len(self.scores_list[i]['original_species']),len(self.scores_list[i]['pident'])):
                pident_organised.append(self.scores_list[i]['pident'][j])
            # Counter method:
            for j in range(len(self.scores_list[i]['original_species'])):
                pident_total = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        pident_total += pident_organised[j]
                pident_total += self.scores_list[i]['pident'][j]
                pident_total /= len(self.scores_list[i]['original_species'])
                mean_identical_match_percentage_initials.append(pident_total)

            # Percentage identicals for extensions:
            # for reference:
            pident_total = 0
            for j in range(len(pident_total_reference[i])):
                pident_total += pident_total_reference[i][j]
            pident_average = pident_total / len(self.scores_list[i]['original_species'])
            mean_identical_match_percentage_extensions.append(pident_average)
            # for originals:
            # Counter method:
            for j in range(len(self.scores_list[i]['original_species'])):
                pident_total = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        pident_total += pident_total_original[i][j]
                pident_total += pident_total_reference[i][j]
                pident_total /= len(self.scores_list[i]['original_species'])
                mean_identical_match_percentage_extensions.append(pident_total)

        # Codon specific variables:
        # codon_frequency = {'G-nad5':,'G-nad4':,'G-cob':,'G-nad6':,'G-nad1':,'G-nad2':,'G-cox2':,'G-atp6':,'G-cox1':,'G-nad3':,'G-cox3':,'G-nad4L':,'G-atp9':}
        # codon_length = {'G-nad5':,'G-nad4':,'G-cob':,'G-nad6':,'G-nad1':,'G-nad2':,'G-cox2':,'G-atp6':,'G-cox1':,'G-nad3':,'G-cox3':,'G-nad4L':,'G-atp9':}
        # We need a list of all possible 64 codons:
        codons_64 = _create_one_hot_array()
        # Now we can calculate their frequencies:
        frequency_list = []
        distance_list = []
        for i in range(len(self.scores_list)): # For each gene:
            frequency_of_gene = 0
            frequency_species = []
            distance_species = []
            codons = split_len(self.scores_list[i]['reference_sequences'],3,0) # For reference specie
            codon_counts = [] # Will store each count
            codon_distances = [] # Will store distances
            for j in range(64): # should also have 64 elements
                codon_counts.append(0)
                codon_distances.append(0)
            # For all elements in codon_64, we have to traverse codons array and count.
            for j in range(len(codons_64)):
                for k in range(len(codons)):
                    if codons_64[j] == codons[k]:
                        codon_counts[j] += 1
            # Now we get all the codon_counts for a gene in a specie.
            total_codons = len(codons)
            frequency_codons = []
            for j in range(len(codon_counts)):
                frequency_codons.append(codon_counts[j]/total_codons)
            frequency_species.append(frequency_codons)
            # distances for reference:
            codons_reverse = codons[::-1]
            for j in range(len(codons_64)): # For each codon in codons_64
                for k in range(len(codons_reverse)):
                    if codons_64[j] == codons_reverse[k]:
                        index = k
                        position = len(codons_reverse) - index
                        codon_distances[j] = position
                        break
            distance_species.append(codon_distances)
            for j in range(len(self.scores_list[i]['original_species'])): # For each non-reference specie:
                # Divide the DNA sequence of the specie into 3-mers:
                codons = split_len(self.scores_list[i]['original_sequences'][j],3,0)
                codon_counts = [] # Will store each count
                codon_distances = []
                for k in range(64): # should also have 64 elements
                    codon_counts.append(0)
                    codon_distances.append(0)
                # For all elements in codon_64, we have to traverse codons array and count.
                for k in range(len(codons_64)):
                    for l in range(len(codons)):
                        if codons_64[k] == codons[l]:
                            codon_counts[k] += 1
                # Now we get all the codon_counts for a gene in a specie.
                total_codons = len(codons)
                frequency_codons = []
                for k in range(len(codon_counts)):
                    frequency_codons.append(codon_counts[k]/total_codons)
                frequency_species.append(frequency_codons)
                # distances for non-reference:
                codons_reverse = codons[::-1]
                for j in range(len(codons_64)): # For each codon in codons_64
                    for k in range(len(codons_reverse)):
                        if codons_64[j] == codons_reverse[k]:
                            index = k
                            position = len(codons_reverse) - index
                            codon_distances[j] = position
                            break
                distance_species.append(codon_distances)
            frequency_list.append(frequency_species)
            distance_list.append(distance_species)


        # ############# Test for Dataset Debug -8- Passed
        # print(len(mean_length_of_extensions)) # Mean extension lengths
        # print("")
        # print(len(mean_similarity_extension)) # Mean bitscores
        # print("")
        # print(len(frequency_evolutionary)) # Frequencies of stop codons
        # print("")
        # print(len(length_of_genes)) # Lengths of genes
        # print("")
        # print(len(mean_e_values_extensions)) # e-values determined from blast
        # print("")
        # print(len(mean_similarity_initials))
        # print("")
        # print(len(mean_e_values_initials))
        # print("")
        # print(len(mean_identical_match_percentage_initials))
        # print("")
        # print(len(mean_identical_match_percentage_extensions))
        # #############

        self.information_dictionary = {'mloe': mean_length_of_extensions, 'fe': frequency_evolutionary, 'log': length_of_genes,
        'mse': mean_similarity_extension, 'msi': mean_similarity_initials, 'mevi': mean_e_values_initials, 'meve': mean_e_values_extensions,
        'mimpi': mean_identical_match_percentage_initials, 'mimpe': mean_identical_match_percentage_extensions, 'flist': frequency_list,
        'dlist': distance_list}

        return self.information_dictionary

    def _give_meaning(self):
        # First create 3 arrays, all specie names, all genes, and their respective positions in each genes:
            # If Mychonastes_homosphaera / nad6 is reference => 0, otherwise 1,2,3,4,.. depending on the place on the original species.
        all_species = [] # N = specie size
        all_genes = [] # M = size of common genes
        respective_positions = [] # N*M = specie size * size of common genes
        all_species.append(self.scores_list[0]['reference_specie'])
        for i in range(len(self.scores_list[0]['original_species'])):
            all_species.append(self.scores_list[0]['original_species'][i])
        for i in range(len(self.scores_list)):
            all_genes.append(self.scores_list[i]['gene_name'])
        for i in range(len(all_genes)):
            for j in range(len(all_species)):
                if all_species[j] == self.scores_list[i]['reference_specie']:
                    respective_positions.append(-1)
                else:
                    for k in range(len(self.scores_list[i]['original_species'])):
                        if all_species[j] == self.scores_list[i]['original_species'][k]:
                            respective_positions.append(k)
                            break

        # print("All genes: ")
        # print(all_genes)
        # print("All species: ")
        # print(all_species)
        # print("Respective positions: ")
        # print(respective_positions)
        # Len: AG/AS/RP = 6/19/114
        # Respective Positions:
        # [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        #  17, 0, 1, 2, 3, 4, -1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        #   17, 0, 1, 2, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        #    17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, -1,
        #     17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, -1, 13, 14, 15, 16,
        #      17, 0, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

        # TODO: Fix all these collection of datas depending on their respective positions.
        # TODO: Also along the way structure as the final data structure.

        one_hot_array = _create_one_hot_array()
        # Stores 64 possible codons by their encoded vector format.
        one_hot_encoded_array = []
        for i in range(64):
            one_hot_encoded_array.append(_one_hot_encode(one_hot_array,one_hot_array[i]))
        # Also encode the data for each genome and for each specie in a numpy vector:
        total_species = len(self.scores_list[0]['original_species']) + 1 # Non-reference species count + 1 reference specie count
        traverse = int(len(self.information_dictionary['mloe']) / total_species) # gene count

        stop_mapper = _stop_mapper()

        df_total = pandas.DataFrame() # This will be the last dataframe when all others will concatenate in the loop.
        for i in range(traverse): # For each gene
            name_gene = self.scores_list[i]['gene_name'] # Gene's name
            name_gene = name_gene[2:]
            df_species = pandas.DataFrame()
            for j in range(total_species): # For all genomes for that gene
            # There should be 64 codons that have information.
                name_specie = all_species[j]
                # Filter name specie and remove gc part.
                names = name_specie.split(' ')
                name_specie = names[0]
                #print(name_specie)
                names_list = []

                respective_place_for_gene = respective_positions[j] + 1 # Determines the place in self.information_dictionary.
                #'mloe': mean_length_of_extensions, 'fe': frequency_evolutionary, 'log': length_of_genes,
                #'mse': mean_similarity_extension, 'msi': mean_similarity_initials, 'mevi': mean_e_values_initials,
                #'meve': mean_e_values_extensions, 'mimpi': mean_identical_match_percentage_initials,
                #'mimpe': mean_identical_match_percentage_extensions,
                mloe_list = []
                fe_list = []
                log_list = []
                mse_list = []
                msi_list = []
                mevi_list = []
                meve_list = []
                mimpi_list = []
                mimpe_list = []
                mloe_value = self.information_dictionary['mloe'][respective_place_for_gene + i*total_species]
                fe_value = self.information_dictionary['fe'][respective_place_for_gene + i*total_species]
                log_value = self.information_dictionary['log'][respective_place_for_gene + i*total_species]
                mse_value = self.information_dictionary['mse'][respective_place_for_gene + i*total_species]
                msi_value = self.information_dictionary['msi'][respective_place_for_gene + i*total_species]
                mevi_value = self.information_dictionary['mevi'][respective_place_for_gene + i*total_species]
                meve_value = self.information_dictionary['meve'][respective_place_for_gene + i*total_species]
                mimpi_value = self.information_dictionary['mimpi'][respective_place_for_gene + i*total_species]
                mimpe_value = self.information_dictionary['mimpe'][respective_place_for_gene + i*total_species]
                for any in range(64):
                    names_list.append(name_specie)
                    mloe_list.append(mloe_value)
                    fe_list.append(fe_value)
                    log_list.append(log_value)
                    mse_list.append(mse_value)
                    msi_list.append(msi_value)
                    mevi_list.append(mevi_value)
                    meve_list.append(meve_value)
                    mimpi_list.append(mimpi_value)
                    mimpe_list.append(mimpe_value)
                # Create Y:
                # Find which codon that gene for that specie stands for:
                codon = stop_mapper[name_specie][name_gene]
                #print(codon)
                # For the list of codons, find where it corresponds
                Y = [] # Create Y list
                for ycount in range(64):
                    if one_hot_array[ycount] == codon:
                        Y.append(1)
                        #print("Location: %s" % ycount)
                    else:
                        Y.append(0)
                #print(Y)

                df = pandas.DataFrame()
                codons_list = one_hot_array
                # df['Specie Name'] = names_list TODO: Unnecessary in Random Forest
                # df['Codons'] = one_hot_array TODO: Open it back with one-hot-encoding
                df['Frequencies'] = self.information_dictionary['flist'][i][respective_place_for_gene]
                df['Distances'] = self.information_dictionary['dlist'][i][respective_place_for_gene]
                # TODO: NOTE: Problem is that prediction is never 1 with these variables.
                # Prediction is always 0.

                df['Mean Extension Length'] = mloe_list
                df['Evolutionary Frequency'] = fe_list
                df['Gene length'] = log_list
                df['Mean Extension Similarity'] = mse_list
                df['Mean Initial Similarity'] = msi_list
                df['Mean Extension E-values '] = meve_list
                df['Mean Initial E-values '] = mevi_list
                df['Mean Extension Match Percent'] = mimpe_list
                df['Mean Initial Match Percent'] = mimpi_list
                df['Y'] = Y # Add Y to data. TODO: Output

                # Since we know we have to only get rows 48/49/56/52/61, we can drop all the others,
                df = df.drop(df.index[[50,51,53,54,55,57,58,59,60,62,63]])
                df = df[48:]
                df_species = pandas.concat([df, df_species])
            # gene_list = []
            # for what in range(16*total_species):
            #     gene_list.append(name_gene)
            # df_species['Gene Name'] = gene_list TODO: Open this back up with one-hot-encoding
            df_total = pandas.concat([df_total, df_species])
        # print(df_total)
        self.df_total = df_total

        # Label is the feature that we want to predict
        labels = np.array(df_total['Y'])
        # print("Labels: ")
        # print(labels)
        # print("")

        # Remove label from feature
        features = df_total.drop(['Y'],axis=1)

        # Save feature names for later use
        feature_list = list(features.columns)
        # print("Feature List: ")
        # print(feature_list)
        # print("")

        # Convert to numpy array
        features = np.array(features)
        # print("Features: ")
        # print(features)
        # print("")

        # Using Skicit-learn to split data into training and testing sets:
        train_features, test_features, train_labels, test_labels = train_test_split(features,labels,test_size = 0.20, random_state = 42)

        # Test if shapes are okay: !! Expecting training features number of columns to match the testing feature
        # number of columns and the number of rows to match for the respective training and testing features. NOTE: Pass.

        # print('Training Features Shape:', train_features.shape)
        # print('Training Labels Shape:', train_labels.shape)
        # print('Testing Features Shape:', test_features.shape)
        # print('Testing Labels Shape:', test_labels.shape)

        # Instantiate model with 1000 decision trees
        rf = RandomForestClassifier(n_estimators = 1000, random_state = 42)
        # Train the model on training data
        rf.fit(train_features,train_labels)
        # print("Model Trained")
        feature_importance(rf,features_list = feature_list)
        # Use the forest's predict method on the test data
        predictions = rf.predict(test_features)
        # Calculate the absolute errors
        mistakes = 0
        for pred in range(len(predictions)):
            print("TL: %d / Pred: %d " % (test_labels[pred], predictions[pred]))
            if int(test_labels[pred]) != int(predictions[pred]):
                mistakes += 1
                print("Mistake: %d" % mistakes)
        # Determine performance metrics
        print("Total mistakes: %d" % mistakes)
        print("Total predictions: %d" % len(predictions))
        percentage_accuracy = ((len(predictions) - mistakes) / len(predictions)) * 100
        # # Calculate and display accuracy
        print("Accuracy:", round(percentage_accuracy,2), '%.')

    def _learn(self):
        pass

def feature_importance(rf, outfile="importance.png", features_list=[]):
    """Show each feature importance"""
    importances = rf.feature_importances_
    if len(features_list) > 0 and len(features_list) != len(importances):
        raise ValueError("Number of features does not fit!")

    indices = np.argsort(importances)[::-1]
    n_feats = len(features_list)
    np.savetxt(outfile + ".txt", np.array([tree.feature_importances_
                                           for tree in rf.estimators_]), delimiter=',', fmt='%1.3e')
    std = np.std(
        [tree.feature_importances_ for tree in rf.estimators_], axis=0)
    plt.figure()
    plt.title("Feature importances")
    print(importances[indices])
    print(features_list)
    plt.bar(range(n_feats), importances[indices], color="b", yerr=std[indices], align="center")
    if len(features_list) > 0:
        features_list = np.asarray(features_list)[indices]
        plt.xticks(range(n_feats), features_list, rotation='vertical')
    plt.xlim([-1, n_feats])
    plt.margins(0.2)

    plt.subplots_adjust(bottom=0.15)
    plt.savefig(outfile, bbox_inches='tight')
