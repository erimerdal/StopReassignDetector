from .utils import reverse_sequence
from collections import defaultdict as ddict
import re


class MasterCollection:
    def __init__(self, infiles):
        """Parse a list of masterfiles and store their informations"""
        # Infile parameter should be a list of Master Files obtained from MFannot results.
        # TODO: Maybe we can include in our code MFannot command line tools so a person can input as fasta.

        self.infiles = infiles
        self.species_names = []
        self.reference_species = []
        self.common_genes = []
        self.mfilemap = []
        self.extended_sequences = []
        # Constructor collects all sequences necessary from the input files and parse them.
        for x in self.infiles:
            m = MasterFile(x)
            self.mfilemap.append(m)
            self.species_names.append(m.sequences['species_name'])
        # Constructor collects all the extensions to those sequences and stores them.
        # Also collects the species to take reference, and genes that are in common.
        self.common_genes, self.reference_species, self.extended_sequences = MasterFile.extendedparser(
            self.infiles, self.get_list_of_dict())
        # Constructor sends all collected information to align and give a meaning to the results.

    def get_list_of_dict(self):
        """Return the collection of MasterFile objects as a list of dict"""
        return [x.sequences for x in self.mfilemap]

    def digest(self):
        return self.get_list_of_dict(), self.extended_sequences, self.reference_species, self.common_genes


class MasterFile:

    def __init__(self, infile):
        # Parameters:
            # kmer_length: An integer which decides how many nucleotides/extension block should be.
            # Increasing kmer_length will make it harder to find a better extension.
            # extension_length: An integer decides that length(original sequence + extension) / length(reference_sequence)
            # threshold: An integer that decides the minimum score that local alignments of extensions can return.
            # gap_open_penalty, match_score, mismatch_penalty: Not used currently because of Blosum62 Matrix scores.
            # allow: An integer boolean that decides if we allow original sequence + extensions to be longer than reference sequence.

        self.infile = infile
        self.sequences = self._internal_mfileparser(self.infile)

    def format_genome(self, sformat='fasta'):
        """Return genome under a specific format
        Only fasta is currently supported
        """
        complete_genomes = ""
        if not sformat == 'fasta':
            raise NotImplementedError('Other format are not implemented')

        for g in self.sequences['genes_list']:
            seq = self.sequences['sequences'].get(g, '')
            cur_header = '>{gname} {specname}'.format(
                gname=g, specname=self.sequences['species_name'])
            pos = self.sequences['gpos'].get(g)
            if pos:
                cur_header += "{size} ({start}:{end})".format(
                    size=len(seq), start=pos[0], end=pos[1])
            complete_genomes += cur_header + "\n" + seq + "\n"

        return complete_genomes

    def _internal_mfileparser(self, infile):
        """Parse a single MasterFile"""
        is_reverse = {}
        genes = []
        gene_tracker = []
        comments = []
        sequences = ddict(str)  # map each gene name to a sequence
        gpos = ddict(tuple)
        master_dict = {}
        speciesname = ''
        species_gc = 1

        with open(infile, 'r') as handle:
            line = handle.readline()
            while line and not line.startswith('>'):
                # Try to scan for the list of potential genes
                # if line.startswith(';;'):
                #   line = line.strip()
                # nevermind, not useful
                line = handle.readline()
                # skip to first line with '>'
            # Set the required specname and gc code for the genome
            if line:
                # skip to genomic seq
                speciesname = line[1:].rstrip()
                species_gc = speciesname.split()[-1]  # last item
                if species_gc and species_gc != speciesname:
                    species_gc = species_gc.split('=')[-1].strip()

                line = handle.readline()
                # Storing genes and if they should be reversed.
                while line:
                    line = line.strip()
                    if line.startswith(';;'):
                        pass
                    elif line.startswith(';'):
                        # Necessary informations are parsed

                        line = line.strip('; ')
                        if ';;' in line:
                            comments.append(line.rsplit(';;')[-1])
                        else:
                            comments.append('')
                        line = line.split(';;')[0].strip('; ')
                        try:
                            genename, updown, startend = line.split()[0:3]
                            startend = startend.split()[0]
                            is_starting = False

                            # We should store the gene in genes with these conditions:
                            #     1- If gene name has ==> and start in it
                            #     2- If gene name has <== and end in it, then reverse it.
                            # We will be removing introns and exons from gene names.
                            if not ("-I" in genename or '-E' in genename):
                                genes.append(genename)
                                if updown == "==>" and startend == "start":
                                    is_reverse[genename] = False
                                    is_starting = True
                                if updown == "<==" and startend == "end":
                                    is_reverse[genename] = True
                                    is_starting = True
                                if genename not in gene_tracker and is_starting:
                                    gene_tracker.append(genename)
                                else:
                                    gene_tracker = [
                                        gn for gn in gene_tracker if gn != genename]

                        except ValueError:
                            pass
                            # this is one of the gene like rnl that we don't need anyway

                    else:
                        # If they are lowercase, this means they belong
                        # to an intron which should not be taken into the sequence.
                        pos, seq = line.split()
                        if not seq.islower():  # sequence is exon
                            for g in gene_tracker:  # if the gene is not removed already, it's its sequence
                                sequences[g] += seq
                                cur_pos = gpos.get(g)
                                if not cur_pos:
                                    gpos[g] = (int(pos), int(pos)+len(seq))
                                else:
                                    gpos[g] = (cur_pos[0], cur_pos[1]+len(seq))
                    line = handle.readline()

                # """ Now we should reverse 5'->3' strands to 3'->5' strands. """
                for g, seq in sequences.items():
                    if is_reverse.get(g):
                        sequences[g] = reverse_sequence(seq)

                master_dict = {'species_name': speciesname, 'species_gc': species_gc,
                               'genes_list': genes, 'sequences': sequences, 'comments': comments, 'gpos': gpos}

            return master_dict

    @classmethod
    def extendedparser(cls, infile, list_of_dictionaries):
         # For each gene that is common, we need to find if the gene is extendable or not.
         #     1- We need to find which genes are common in all species, and store them.
         #     2- Filter the interesting genes
         # List of interesting genes: nad / atp / cob / cox / sdh / tat can be extended. Other genes are omitted.
        common_genes_3x = []
        for x in range(len(list_of_dictionaries[0]['genes_list'])):
            for y in range(1, len(list_of_dictionaries)):
                for z in range(len(list_of_dictionaries[y]['genes_list'])):
                    if list_of_dictionaries[0]['genes_list'][x] == list_of_dictionaries[y]['genes_list'][z]:
                        if "G-nad" in list_of_dictionaries[0]['genes_list'][x] or "G-atp" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(
                                list_of_dictionaries[0]['genes_list'][x])
                        if "G-cob" in list_of_dictionaries[0]['genes_list'][x] or "G-cox" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(
                                list_of_dictionaries[0]['genes_list'][x])
                        if "G-sdh" in list_of_dictionaries[0]['genes_list'][x] or "G-tat" in list_of_dictionaries[0]['genes_list'][x]:
                            common_genes_3x.append(
                                list_of_dictionaries[0]['genes_list'][x])
        # no need of additional function
        common_genes = list(set(common_genes_3x))
        # 3- Now that we have genes that are common and interesting, among all files we need to find
        #      the longest sequenced ones for that specific gene, to use it as a reference gene to compare
        #      with others later on.
        reference_species = []
        sequence_lengths = []
        for x in range(len(common_genes)):
            for y in range(len(infile)):
                for z in range(len(list_of_dictionaries[y]['genes_list'])):
                    if list_of_dictionaries[y]['genes_list'][z] == common_genes[x]:
                        sequence_lengths.append(
                            len(list_of_dictionaries[y]['sequences'][z]))

        # 4- Now find the longest of these for each gene and store which specie it is.
        #      4.a- Divide in chunks of size species.
        chunks = [sequence_lengths[x*len(infile):x*len(infile)+len(infile)]
                  for x in range(0, int(len(sequence_lengths)/len(infile)))]
        #      4.b- Find longest for each and store which specie has longest.
        # TODO: Change this to dynamic finding later. Sorry for this..
        # TODO: Reference_specie is automatically Bracteacoccus_minor now for every gene since its the only gc = 22.
        for x in range(len(chunks)):
            maximum = max(chunks[x])
            for y in range(len(chunks[x])):
                if maximum == chunks[x][y]:
                    reference_species.append(
                        list_of_dictionaries[y]['species_name'])
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
            handle = open(infile[x], 'r')
            line = handle.readline()
            line_number = line_number + 1
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
                            line_number = line_number + 1
                        else:
                            line = handle.readline()
                            line_number = line_number + 1
                    else:
                        line = handle.readline()
                        line_number = line_number + 1
                else:
                    line = handle.readline()
                    line_number = line_number + 1
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
                    handle = open(infile[x], 'r')
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
                                        reverse_handle = open(infile[x], 'r')
                                        from_line = 0
                                        to_line = 0
                                        for i in range(len(reverse_gene_order[x])):
                                            if reverse_gene_order[x][i] == genename:
                                                to_line = line_numbers[x][i]
                                                from_line = line_numbers[x][i-1]
                                                break

                                        for i in range(from_line):
                                            reverse_line = reverse_handle.readline()

                                        for i in range(from_line, to_line):
                                            string = genename + " "
                                            if string in reverse_line:
                                                break
                                            reverse_extension.append(
                                                reverse_line.rstrip())
                                            if reverse_line.startswith(';;') or reverse_line.startswith(';'):
                                                reverse_extension.remove(
                                                    reverse_line.rstrip())
                                            reverse_line = reverse_handle.readline()
                                        sequence = "".join(reverse_extension).replace(
                                            " ", "").replace("\r", "")
                                        sequence = re.sub("\d", "", sequence)
                                        # Reversing extensions seem like a better option
                                        sequence_reversed = reverse_sequence(
                                            sequence)
                                        extended_sequences.append(
                                            sequence_reversed[::-1])

                                start_end_count = start_end_count + 1

                                # This is for the non-reverse gene extensions
                                if start_end_count == 2:
                                    if startend == "end":
                                        # Go until another gene starts
                                        line = handle.readline()
                                        genes.append(genename)
                                        while True:
                                            if line.startswith(';'):
                                                geneName = line[1:].rstrip()
                                                fullgenename = " ".join(
                                                    geneName.split())
                                                iterate = fullgenename.split(
                                                    " ")
                                                # Necessary informations are parsed
                                                genename = iterate[0]
                                                updown = iterate[1]
                                                startend = iterate[2]
                                                if startend == "start":
                                                    if updown == "==>":
                                                        break
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
                                                    extension.remove(
                                                        line.rstrip())
                                                line = handle.readline()
                                        # Format the sequence in a better, readable format.
                                        sequence = "".join(extension).replace(
                                            " ", "").replace("\r", "")
                                        sequence = re.sub("\d", "", sequence)

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
            extended_dict = {
                'species_name': list_of_dictionaries[x]['species_name'], 'genes_list': genes, 'sequences': extended_sequences}

            # Add each dictionary / per specie to a list of dictionaries.
            list_of_extended_dictionaries.append(extended_dict)
        return common_genes, reference_species, list_of_extended_dictionaries
