from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from pandas import *
from skbio.alignment import local_pairwise_align_ssw, make_identity_substitution_matrix, global_pairwise_align_protein
from skbio import TabularMSA, DNA, Protein
from .Parser import MasterFile, MasterCollection
from .utils import split_len
import os
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess

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

        # Here we will start a new, final chapter everyone!
        # First compare all sequences with reference sequence.
        for i in range(len(scores_list)):
            pairs_results_extensions = []
            pairs_results_extensions_scores = []
            # Pairwise reference with each non-reference specie.
            for j in range(len(scores_list[i]['original_species'])):
                # Find the place left in reference sequence.
                end_position_reference = scores_list[i]['start_end'][j][1][1]
                # Find end position of original sequence.
                end_position_original = scores_list[i]['start_end'][j][0][1]
                # Take whats left of reference from initial local alignment.
                reference_sequence_right = scores_list[i]['reference_sequence_protein'][end_position_reference:]
                # Concatenate the last part of original sequence and the extensions for them.
                concatenated_original_seq = scores_list[i]['original_sequences_protein'][j][end_position_original:]
                concatenated = concatenated_original_seq + scores_list[i]['extensions_protein'][j]
                # Protein format of the original extension.
                newStringProtein = Protein(str(concatenated))
                # Protein format of the part left in reference.
                refseq = Protein(str(reference_sequence_right))
                # Search reference sequence in original sequence with LOCAL? alignment.
                # print("Length of ref_seq: %s" % len(refseq))
                # print("Length of org_seq: %s" % len(newStringProtein))
                try:
                    alignment, score, start_end_positions = local_pairwise_align_ssw(
                    refseq,
                    newStringProtein,
                    substitution_matrix = submat,
                    )
                    # Calculate start end position with respect to beginning.
                    # start_location = start of the extension + other extensions before * their length + end of local alignment
                    start_location = start_end_positions[0][0] + scores_list[i]['start_end'][j][0][1]
                    # start_location_r = end location of initial alignment + start position
                    start_location_r = start_end_positions[1][0] + scores_list[i]['start_end'][j][1][1]
                    # end_location = end of the extension + other extensions before * their length + end of local alignment
                    end_location = start_end_positions[0][1] + scores_list[i]['start_end'][j][0][1]
                    end_location_r = start_end_positions[1][1] + scores_list[i]['start_end'][j][1][1]
                    start_end_position_reference = []
                    start_end_position_reference.append(start_location_r)
                    start_end_position_reference.append(end_location_r)
                    start_end_position_original = []
                    start_end_position_original.append(start_location)
                    start_end_position_original.append(end_location)
                    # "(%s,%s),(%s,%s)" % (start_location_r,end_location_r,start_location,end_location)
                    start_end_position = []
                    start_end_position.append(start_end_position_reference)
                    start_end_position.append(start_end_position_original)
                except IndexError:
                    score = 0
                    start_end_position = "None"
                pairs_results_extensions.append(start_end_position)
                pairs_results_extensions_scores.append(score)
            data3_temporary.append(pairs_results_extensions)
            data5_temporary.append(pairs_results_extensions_scores)

        # pairwise.
        for i in range(len(scores_list)):
            pairs_results_extensions_j = []
            pairs_results_extensions_j_score = []
            for j in range(len(scores_list[i]['original_species'])-1):
                for k in range(j+1,len(scores_list[i]['original_species'])):
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

                    # Now we have to extend them both in the double for loop and collect scores.
                    j_sequence = Protein(str(j_sequence_total))
                    k_sequence = Protein(str(k_sequence_total))
                    try:
                        alignment, score, start_end_positions = local_pairwise_align_ssw(
                        j_sequence,
                        k_sequence,
                        substitution_matrix = submat,
                        )
                        # Calculate start end position with respect to beginning.
                        # start_location = start of the extension + other extensions before * their length + end of local alignment
                        start_location = start_end_positions[0][0] + scores_list[i]['start_end'][j][0][1]
                        # start_location_r = end location of initial alignment + start position
                        start_location_r = start_end_positions[1][0] + scores_list[i]['start_end'][j][1][1]
                        # end_location = end of the extension + other extensions before * their length + end of local alignment
                        end_location = start_end_positions[0][1]+ scores_list[i]['start_end'][j][0][1]
                        end_location_r = start_end_positions[1][1] + scores_list[i]['start_end'][j][1][1]
                        start_end_position_reference = []
                        start_end_position_reference.append(start_location_r)
                        start_end_position_reference.append(end_location_r)
                        start_end_position_original = []
                        start_end_position_original.append(start_location)
                        start_end_position_original.append(end_location)
                        # "(%s,%s),(%s,%s)" % (start_location_r,end_location_r,start_location,end_location)
                        start_end_position = []
                        start_end_position.append(start_end_position_reference)
                        start_end_position.append(start_end_position_original)
                    except IndexError:
                        score = 0
                        start_end_position = "None"

                    pairs_results_extensions_j_score.append(score)
                    pairs_results_extensions_j.append(start_end_position)
                    # print("j: %s k: %s" % (j,k)) seems okay
            data4_temporary.append(pairs_results_extensions_j)

            data6_temporary.append(pairs_results_extensions_j_score)

        for i in range(len(data6_temporary)):
            for k in range(len(data4_temporary[i])):
                data3_temporary[i].append(data4_temporary[i][k])

        # TODO:TODO: Storing variable informations for modelling DT/RF.
        mean_length_of_extensions = []
        length_of_genes = []
        frequency_of_stop_codon = []
        stop_codons_evolutionary = []
        for i in range(len(scores_list)): # For each gene // Genes should be taken apart seperately while collecting data
            # The length of the extension
            # For the reference specie:
            # We should do the same order when we were calculating scores/start.end positions.
            reference_extension_total = 0
            for j in range(len(scores_list[i]['original_species'])):
                reference_extension_total += data3_temporary[i][j][0][1] - data2_temporary[i][j][0][1]
            mean_length_of_extensions.append(reference_extension_total/(len(scores_list[i]['original_species'])))
            # We need a better version of data structure, too much calculation involved while trying to find all non-references.
            # The data structure should be storing all information for all genes for each specie seperately.

            for j in range(len(scores_list[i]['original_species'])):
                mloe_dataset = []
                mloe_dataset_in = []
                # in the order of the original species
                for k in range(len(data1_temporary[i])):
                    data1_split = data1_temporary[i][k].split('/')
                    if scores_list[i]['original_species'][j] in data1_split[0]: # Left side
                        mloe_dataset.append(data3_temporary[i][k][0][1]) # Left end
                        mloe_dataset_in.append(data2_temporary[i][k][0][1])
                    if scores_list[i]['original_species'][j] in data1_split[1]: # Right side
                        mloe_dataset.append(data3_temporary[i][k][1][1]) # Right end
                        mloe_dataset_in.append(data2_temporary[i][k][1][1])
                mloe = 0
                # Now that we have the dataset for each gene, we can keep calculating.
                for k in range(len(mloe_dataset)):
                    mloe += mloe_dataset[k] - mloe_dataset_in[k]
                mean_length_of_extensions.append(mloe/(len(scores_list[i]['original_species'])))

            # The length of the gene compared to other homologs:
            length_of_genes.append(len(scores_list[i]['reference_sequences']))
            for j in range(len(scores_list[i]['original_species'])): # For each non-reference-specie
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
        # print(len(mean_length_of_extensions)) = 65, 13 genes for 5 species correct.
        # Frequency of the stop codon in all evolutionary close stop codons
        frequency_evolutionary = []
        mean_similarity_extension = []
        counter = []
        onlyOnce = True # Store counter only once not genes times.
        for i in range(len(scores_list)):
            if onlyOnce:
                for j in range(len(scores_list[i]['original_species'])-1):
                    for k in range(j+1,len(scores_list[i]['original_species'])):
                        counter.append("(%s,%s)" % (j,k))
                onlyOnce = False

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
            meansim_ref = 0
            for j in range(len(data5_temporary[i])):
                meansim_ref += data5_temporary[i][j]
            meansim_ref /= len(scores_list[i]['original_species'])
            mean_similarity_extension.append(meansim_ref) # Reference species similarities are stored.

            # Implement mean similarity extension for non-reference sequences.
            for j in range(len(scores_list[i]['original_species'])):
                meansim_org = 0
                for k in range(len(counter)):
                    if str(j) in counter[k]:
                        meansim_org += data6_temporary[i][j]
                meansim_org /= len(scores_list[i]['original_species'])
                mean_similarity_extension.append(meansim_org)

        ############## Test for BLAST
        seq1 = SeqRecord(Seq("FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF"),
                   id="seq1")
        seq2 = SeqRecord(Seq("FQTWEEFSRAEKLYLADPMKVRVVLRYRHVDGNLCIKVTDDLICLVYRTDQAQDVKKIEKF"),
                   id="seq2")
        SeqIO.write(seq1, "seq1.fasta", "fasta")
        SeqIO.write(seq2, "seq2.fasta", "fasta")

        # Run BLAST and parse the output
        command = './blastp -outfmt 6 -query seq1.fasta -subject seq2.fasta'
        p = subprocess.Popen(command , shell=True, stdout=subprocess.PIPE)
        print(p.stdout.readlines()[0].decode('utf-8'))

        ############## Test for Dataset
        # print(mean_length_of_extensions)
        # print("")
        # print(mean_similarity_extension)
        # print("")
        # print(frequency_evolutionary)
        # print("")
        # print(length_of_genes)
        ##############

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


        # TODO: Do we have to take into account the distance scoring function? Maybe so, how? Distance is not included,
        # extension is not trimmed until the length of the reference, sometimes gives stupidly far results therefore.
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
