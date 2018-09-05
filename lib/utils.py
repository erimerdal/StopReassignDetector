import os
import re
import numpy as np

REV_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def reverse_sequence(sequence):
    """Return the reversed sequence (5'->3'::3'->5')"""
    return "".join(REV_DICT.get(nuc, nuc) for nuc in sequence[::-1])


# split_len(): Splits the sequence "seq" into "length" amount of seperate sequences.
def split_len(seq, length, position):
    return [seq[i:i+length] for i in range(position, len(seq), length)]


# getSubStrings(): Divides the given sequence "RNA" into 3-mers i.e [0-3,1-4,2-5 etc.]. Position should be 0,1 or 2.
def getSubStrings(RNA, position):
    return [RNA[i:i+3] for i in range(position, len(RNA), 3)]

def list_filter(l1, reflist):
    return [x for y in reflist for x in l1  if y in x]

def _one_hot_encode(one_hot_array, codon):
    one_hot_vector = np.array([])
    for i in range(len(one_hot_array)):
        if one_hot_array[i] == codon:
            one_hot_vector = np.append(one_hot_vector,1)
        else:
            one_hot_vector = np.append(one_hot_vector,0)
    return one_hot_vector

def _create_one_hot_array():
    one_hot_array = []
    one_hot_array.append('AAA')
    one_hot_array.append('AAG')
    one_hot_array.append('AAC')
    one_hot_array.append('AAT')
    one_hot_array.append('AGA')
    one_hot_array.append('AGG')
    one_hot_array.append('AGC')
    one_hot_array.append('AGT')
    one_hot_array.append('ACA')
    one_hot_array.append('ACG')
    one_hot_array.append('ACC')
    one_hot_array.append('ACT')
    one_hot_array.append('ATA')
    one_hot_array.append('ATG')
    one_hot_array.append('ATC')
    one_hot_array.append('ATT')

    one_hot_array.append('GAA')
    one_hot_array.append('GAG')
    one_hot_array.append('GAC')
    one_hot_array.append('GAT')
    one_hot_array.append('GGA')
    one_hot_array.append('GGG')
    one_hot_array.append('GGC')
    one_hot_array.append('GGT')
    one_hot_array.append('GCA')
    one_hot_array.append('GCG')
    one_hot_array.append('GCC')
    one_hot_array.append('GCT')
    one_hot_array.append('GTA')
    one_hot_array.append('GTG')
    one_hot_array.append('GTC')
    one_hot_array.append('GTT')

    one_hot_array.append('CAA')
    one_hot_array.append('CAG')
    one_hot_array.append('CAC')
    one_hot_array.append('CAT')
    one_hot_array.append('CGA')
    one_hot_array.append('CGG')
    one_hot_array.append('CGC')
    one_hot_array.append('CGT')
    one_hot_array.append('CCA')
    one_hot_array.append('CCG')
    one_hot_array.append('CCC')
    one_hot_array.append('CCT')
    one_hot_array.append('CTA')
    one_hot_array.append('CTG')
    one_hot_array.append('CTC')
    one_hot_array.append('CTT')

    one_hot_array.append('TAA')
    one_hot_array.append('TAG')
    one_hot_array.append('TAC')
    one_hot_array.append('TAT')
    one_hot_array.append('TGA')
    one_hot_array.append('TGG')
    one_hot_array.append('TGC')
    one_hot_array.append('TGT')
    one_hot_array.append('TCA')
    one_hot_array.append('TCG')
    one_hot_array.append('TCC')
    one_hot_array.append('TCT')
    one_hot_array.append('TTA')
    one_hot_array.append('TTG')
    one_hot_array.append('TTC')
    one_hot_array.append('TTT')
    return one_hot_array

def _stop_mapper():
    stop_mapper = {"Sorghum_bicolor": {"nad1": "TAA", "atp1": "TGA", "nad9": "TAA", "nad3": "TAA", "rps12": "TGA", "rps2b": "TGA", "cox3": "TGA", "rps13": "TGA", "mttb": "TAG",
     "nad5": "TAA", "mat": "TAG", "rps1": "TAG", "ccmfn": "TAG", "nad6": "TAG", "nad7": "TAG", "rps7": "TAA", "cob": "TAG", "ccmc": "TAG", "atp8": "TAA", "nad2": "TAA", "rps4": "TAA",
      "cox1": "TAA", "nad4l": "TAA", "ccmb": "TGA", "atp4": "TAA", "cox2": "TAA", "atp6": "TAA", "nad4": "TGA", "atp9": "TAG", "rpl16": "TAA", "rps3": "TAG", "ccmfc": "CGA"},

    "Brassica_rapa_subsp_oleifera": {"cox2": "TAA", "orf101f": "TAA", "ccmc": "TGA", "nad2": "TAA", "rps4": "TAA", "nad7": "TAG", "orf120": "TAG", "orf115a": "TGA", "orf114": "TGA",
     "atp8": "TGA", "orf118": "TGA", "orf293": "TGA", "orf113a": "TAG", "orf110": "TAA", "orf100a": "TAA", "orf101a": "TGA", "orf123": "TGA", "orf113b": "TAA", "nad5": "TAA",
      "nad9": "TAA", "ccmfn2": "TGA", "atp6": "TAA", "atp4": "TAA", "nad4l": "TAA", "orf101b": "TAA", "orf265": "TAA", "cox1": "TAA", "orf257": "TAG", "orf101e": "TAG",
       "rps3": "TAG", "rpl16": "TAA", "rpl5": "TAA", "rps14": "TAG", "cob": "TGA", "orf115c": "TGA", "orf108a": "TAA", "ccmfc": "TAA", "orf146": "TAA", "cox3": "TGA",
        "orf147": "TAG", "orf164": "TGA", "orf117a": "TGA", "orf100b": "TAG", "orf108b": "TAA", "orf115d": "TAG", "ccmfn1": "TGA", "atp9": "TGA", "orf305": "TAA", "nad3": "TAA",
         "rps12": "TGA", "ccmb": "TAA", "orf101c": "TGA", "orf132": "TAG", "orf122": "TGA", "nad1": "TAA", "orf135": "TAA", "nad6": "TAA", "orf448": "TAG", "nad4": "TGA",
          "orf116": "TAA", "tatc": "TAG", "rpl2": "TGA", "atp1": "TAG", "orf109": "TAA", "orf101d": "TAG", "orf128": "TGA", "orf106b": "TAG", "orf112": "TAG", "rps7": "TAA",
           "orf119": "TGA", "orf125": "TAA", "orf161": "TGA", "orf159": "TGA", "matr": "TAG", "orf115b": "TAG", "orf104": "TGA", "orf195": "TGA"},

    "Arabidopsis_thaliana": {"nad2": "TAG", "nad5": "TAA", "nad9": "TAA", "rpl16": "TAA", "rps3": "TAG", "ccb206": "TAA", "cox2": "TAA", "ccb452": "TAA", "rpl5": "TAA",
     "cob": "TGA", "nad6": "TAA", "rps4": "TAA", "atp6": "TAA", "orfb": "TGA", "nad7": "TAG", "nad1": "TGA", "matr": "TAG", "rpl2": "TGA", "orfx": "TGA",
      "nad4": "TGA", "orf25": "TAA", "nad4l": "TAA", "cox3": "TGA", "ccb382": "TGA", "ccb256": "TAA", "ccb203": "TGA", "rpsl2": "TGA", "nad3": "TAA", "atp9": "TGA",
       "atp1": "TAG", "rps7": "TAA", "cox1": "TAA"},

    "Chara_vulgaris": {"nad4l": "TAA", "rps4": "TAA", "mttb": "TGA", "atp9": "TAA", "orf576": "TAA", "rps11": "TAA", "rpl6": "TAA", "rps14": "TAG", "rpl5": "TAG",
     "rpl14": "TAA", "rpl16": "TAA", "rps3": "TAA", "rps19": "TAA", "rpl2": "TAA", "rps10": "TGA", "nad7": "TAA", "nad3": "TAA", "orf558": "TAG", "orf760": "TAG",
      "atp1": "TAA", "nad9": "TAG", "nad5": "TAA", "nad4": "TAA", "nad2": "TAA", "rps7": "TAA", "rps12": "TGA", "atp6": "TAA", "nad6": "TAG", "cox2": "TAA",
       "cox3": "TAA", "nad1": "TAA", "cob": "TAG", "orf448": "TAA", "yejr": "TAA", "yeju": "TAA", "yejv": "TAA", "rps2": "TAA", "rps1": "TAG", "atp8": "TGA",
        "ymf39": "TGA", "cox1": "TAA", "orf389": "TGA", "orf550": "TGA", "sdh4": "TAA", "sdh3": "TAA", "orf259": "TGA"},

    "Beta_vulgaris_subsp_vulgaris": {"nad9": "TAA", "orf155a": "TAA", "orf114a": "TGA", "orf146": "TAA", "orf106a": "TAG", "orf315": "TAA", "orf126": "TAG", "orf129a": "TAG", "orf133a": "TAG",
     "orf329": "TAG", "orf152": "TAG", "orf105a": "TAG", "orf114b": "TGA", "ccb577": "TGA", "ccb438": "TGA", "orf189": "TGA", "orf103a": "TAG", "orf165": "TAA", "orf518": "TAG",
      "orf25": "TAG", "nad4l": "TAA", "orf116": "TGA", "orf113": "TGA", "orf764": "TAA", "orf100a": "TAG", "orf104a": "TGA", "atp9": "CGA", "orf107a": "TAA", "nad2": "TAA",
       "orf300": "TAA", "orf107b": "TAA", "orf122": "TAA", "orf104b": "TAA", "orf199": "TAA", "orf111": "TAA", "orf103b": "TAA", "nad1": "TAA", "orf104c": "TGA", "orf117": "TAG",
        "rps13": "TGA", "orf215": "TAA", "orf106b": "TGA", "orf310": "TAA", "orf409": "TAG", "orf148": "TAA", "orf348": "TAA", "orf670": "TGA", "orf119a": "TAA", "orf100b": "TAG",
         "orf134": "TAG", "nad5": "TAA", "orf129b": "TAG", "orf145": "TGA", "orf109": "TAA", "orf120": "TGA", "orf147": "TGA", "orf129c": "TAG", "ccb206": "TGA", "orf162": "TAA",
          "orf119b": "TAG", "orf170": "TAA", "orf110a": "TAG", "orf115a": "TAA", "orf317": "TAG", "orf399": "TGA", "orf214": "TAG", "cob": "TGA", "orf112a": "TAA", "rpl5": "TAA",
           "orf108a": "TAG", "rps4": "TAA", "orf169": "TAA", "nad6": "TAA", "orf100c": "TAG", "orfx": "TAG", "orf115b": "TAG", "orf114c": "TGA", "atp6": "CAA", "orf100d": "TAG",
            "nad7": "TAG", "orf155b": "TGA", "orf138": "TGA", "orf121": "TGA", "rps7": "TAA", "orf102a": "TGA", "cox2": "TAA", "orf270": "TAG", "cox1": "TAA", "mat": "TAG",
             "orf198": "TAA", "orf125": "TAA", "orf123": "TAA", "orf102b": "TGA", "orf112b": "TAA", "orf143a": "TGA", "orf192": "TGA", "orf100e": "TGA", "nad3": "TAA", "rps12": "TGA",
              "orf100f": "TGA", "orf187": "TAA", "rps3": "TAA", "orf246": "TAG", "nad4": "TGA", "orf288": "TAA", "orf190": "TGA", "orf352": "TAG", "orf217": "TAA", "orf166": "TGA",
               "orf202": "TAG", "orf290": "TAA", "cox3": "TGA", "orf136": "TGA", "atp8": "TAA", "orf103c": "TAA", "orf131": "TGA", "orf211": "TGA", "orf245": "TGA", "orf105b": "TAG",
                "orf105c": "TAA", "orf133b": "TAA", "orf108b": "TAG", "orf393": "TGA", "orf204": "TGA", "orf143b": "TAA", "orf115c": "TAG", "orf135": "TAA", "orf110b": "TGA", "atp1": "TAA",
                 "orf124": "TGA", "orf114d": "TAA", "orf119c": "TAA", "orf106c": "TGA"},

    "Vitis_vinifera": {"nad5": "TAA", "nad1": "TAA", "matr": "TAG", "peta": "TAG", "ycf4": "TGA", "rbcl": "TAA", "orf106": "TAA", "psbj": "TAG", "petl": "TGA",
     "petg": "TGA", "ccmfn": "TGA", "nad7": "TAG", "orf185": "TAA", "nad2": "TAA", "ccmfc": "TAG", "rna_pol": "TAA", "orf159": "TAA", "sdh3": "TGA", "rpl5": "TAA",
      "rps14": "TAG", "cob2": "TGA", "cox1": "TAA", "orf104": "TGA", "rps7": "TAA", "atp6": "TGA", "rps15": "TAA", "psba": "TAA", "nad4l": "TAA", "atp4": "TAG",
       "nad9": "TAA", "ccmc": "TGA", "rpl2": "TAA", "rps19": "TAA", "rps3": "TAG", "rpl16": "TAG", "orf333": "TAA", "rps13": "TGA", "orf102": "TAG", "rpl32": "TAA",
        "psac": "TGA", "ndhe": "TAG", "rps10": "TGA", "atp8": "TAA", "cox3": "TGA", "sdh4": "TGA", "rps1": "TAA", "rps12": "TGA", "nad3": "TAA", "ccmb": "TGA",
         "rpl20": "TAA", "rps18": "TAG", "rpl33": "TAG", "psaj": "TAA", "psbt": "TAG", "psbn": "TAG", "atpa": "TAA", "atp1": "TGA", "orf86": "TGA", "nad6": "TAA",
          "rps4": "TAA", "rpl14": "TAA", "infa": "TAG", "rpl36": "TAA", "rps11": "TAG", "cox2": "TAA", "atp9": "TGA", "psbm": "TAA", "nad4": "TGA", "orf145": "TAA"},

    "Glycine_max": {"nad2": "TAA", "ccmc": "TGA", "cox3": "TGA", "rp15": "TAA", "rps14": "TAG", "cob": "TGA", "nad3": "TAA", "rps12": "TGA", "atp4": "TAA",
     "rps10": "TGA", "cox1": "TAA", "nad1": "TAA", "ccmfc": "CGA", "atp6": "TAG", "nad4": "TGA", "nad4l": "TAA", "nad6": "TGA", "ccmb": "TGA", "nad7": "TAG",
      "atp1": "TGA", "ccmfn": "TGA", "rps4": "TAA", "rps1": "TAA", "nad5": "TAA", "nad9": "TAA", "rpl16": "TAA", "rps3": "TAG", "atp9": "TAA", "matr": "TGA",
       "cox2": "TAA", "atp8": "TAA", "mttb": "TGA"},

    "Malus_domestica": {"ccmfc": "TGA", "cox1": "TAA", "nad4": "TGA", "nad5": "TAA", "ccmc": "TGA", "ccmb": "TGA", "atp9": "TAA", "ccmfn": "TAG", "nad1": "TAA",
     "rpl5": "TAA", "rps14": "TAG", "cob": "TGA", "rps12": "TGA", "nad3": "TAA", "nad2": "TAA", "nad6": "TAA", "atp1": "TGA", "nad7": "TAG", "rps1": "TAA",
      "matr": "TAG", "rps13": "TGA", "nad9": "TAA", "atp6": "TAG", "sdh4": "TAA", "cox3": "TGA", "atp8": "TAA", "rps3": "TGA", "rnaseh": "TAG", "cox2": "TAA",
       "nad4l": "TAA", "atp4": "TAG"},

    "Oryza_sativa_Indica_Group": {"nad1": "TAA", "cox3": "TGA", "nd5": "AAA", "rps7": "TAA", "nad6": "TAA", "ccmc": "TAG", "nad7": "TAG", "rps3": "TAG", "rpl16": "TAA",
     "nad3": "TAA", "rps12": "TGA", "rps2": "TAG", "nad4": "TGA", "cox2": "TAA", "atp6": "TAA", "rps13": "TGA", "rps4": "TAA", "atp9": "CGA", "rpl2": "TAG",
      "rps19": "TAA", "nad4l": "TAA", "cob": "TAG", "mat": "TAG", "rps1": "TGA", "ccmfn": "TAG", "ccmfc": "AAT", "cox1": "TAG", "rpl5": "TAA", "ccmb": "TGA",
       "nad2": "TAA", "atp1": "TAG", "nad9": "TAA"},

    "Marchantia_polymorpha_subsp_ruderalis": {"nad5": "TAA", "nad4": "TAA", "nad2": "TGA", "rps12": "TGA", "rps7": "TAA", "atp6": "TAA", "nad6": "TGA", "sdh3": "TAA", "nad3": "TAA",
     "rps10": "TAG", "rpl2": "TAG", "rps19": "TAG", "rps3": "TAA", "rpl16": "TAA", "rpl5": "TAG", "rps14": "TAG", "rps8": "TAA", "rpl6": "TAA", "rps13": "TGA",
      "rps11": "TAA", "rps1": "TAA", "atp8": "TGA", "sdh4": "TAA", "nad4l": "TAA", "tatc": "TAA", "cox2": "TAG", "cox3": "TAG", "nad1": "TAA", "cob": "TGA",
       "mat": "TGA", "nad9": "TAA", "atp1": "TAA", "cox1": "TAG", "atp4": "TGA", "rps4": "TAA", "atp9": "TAA", "rps2": "TGA", "ccmb": "TGA", "ccmc": "TGA",
        "ccmfn": "TAA", "ccmfc": "TAA", "rpl10": "TGA"},

    "Zea_mays_subsp_mays": {"rps3": "TAG", "rps4": "TAA", "nad4l": "TAA", "nad1": "TAA", "nad2": "TAA", "nad4": "TGA", "ccmfc": "CGA", "nad9": "TAA", "atp4": "TAG",
     "rpl16": "TAA", "ccmc": "TAG", "nad5": "TAA", "atp8": "TAA", "nad7": "TAG", "rps12": "TGA", "nad3": "TAA", "mat": "TAG", "rps1": "TAA", "ccmfn": "TGA",
      "atp9": "TAA", "cob": "TAG", "rps13": "TGA", "rps7": "TAA", "ccmb": "TGA", "cox1": "TAG", "rps2b": "TGA", "cox3": "TGA", "atp1": "TGA", "rps2a": "TAA",
       "cox2": "TAA", "atp6": "TAG", "nad6": "TAA"},

    "Lotus_japonicus": {"cox1": "TAA", "nad1": "TAA", "mttb": "TAG", "nad7": "TAG", "matr": "TAG", "rpl16": "TAA", "rps3": "TAG", "nad2": "TAA", "nad5": "TAA",
     "rps14": "TAG", "rpl5": "TAA", "ccmfn": "TGA", "atp4": "TAA", "nad4l": "TAA", "atp6": "CAA", "nad4": "TGA", "rps4": "TAA", "atp1": "TAA", "ccmfc": "CGA",
      "atp9": "TAA", "nad6": "TAA", "nad9": "TAA", "cob": "TGA", "ccmb": "TGA", "cox3": "TGA", "atp8": "TAA", "nad3": "TAA", "rps12": "TGA", "ccmc": "TAG",
       "cox2": "TGA", "rps10": "CGA"},

    "Phoenix_dactylifera": {"nad2": "TAA", "atp6": "TGA", "rpl5": "TAA", "rps14": "TAG", "cob": "TGA", "rna_pol": "TAA", "rps7": "TAA", "cox1": "TAA", "nad7": "TAG",
     "orf186": "TAA", "rps13": "TGA", "nad1": "TAA", "nad6": "TAG", "rpl2": "TAG", "rps19": "TAA", "rps3": "TAG", "rpl16": "TAA", "atp4": "TAA", "nad4l": "TAA",
      "ccmfc": "TAA", "atp1": "TAA", "atp9": "TAA", "rps1": "TAA", "nad4": "TGA", "cox3": "TGA", "rps2": "TAG", "rps4": "TAA", "rps11": "TAA", "ccmfn": "TAG",
       "mttb": "TAG", "cox2": "TAA", "nad5": "TAA", "ccmb": "TGA", "orf100": "TGA", "matr": "TAG", "orf192": "TAA", "atp8": "TAA", "orf142": "TGA", "nad9": "TAA",
        "ccmc": "TAG", "nad3": "TAA", "rps12": "TGA", "orf118": "TGA"},

    "Chlorokybus_atmophyticus": {"orf170": "TAG", "orf274": "TGA", "nad3": "TAA", "atp4": "TAA", "atp8": "TAA", "rps1": "TAG", "nad4l": "TAA", "sdh4": "TAA", "sdh3": "TAA",
     "atp1": "TAA", "rps11": "TAA", "rps12": "TAA", "rps7": "TAA", "cob": "TAA", "orf260": "TAA", "atp9": "TAG", "orf845": "TAA", "nad10": "TAA", "nad4": "TAA",
      "atp6": "TAA", "nad6": "TAA", "cox2": "TAG", "cox3": "TAA", "rps3": "TAA", "rpl16": "TAA", "rpl14": "TAA", "rpl2": "TAA", "rps19": "TAA", "rps10": "TAA",
       "nad7": "TAA", "rpl5": "TAA", "mttb": "TAA", "rps4": "TAA", "rps2": "TAA", "nad5": "TAG", "rpl6": "TAA", "rps13": "TAA", "rps14": "TAG", "rps8": "TAG",
        "nad1": "TAA", "nad2": "TAA", "nad9": "TAA", "cox1": "TAA", "orf301": "TGA", "orf755": "TAG", "orf296": "TGA"},

    "Marchantia_paleacea": {"nad5": "TAA", "nad4": "TAA", "nad2": "TGA", "rps12": "TGA", "rps7": "TAA", "atp6": "TAA", "nad6": "TGA", "sdh3": "TAA", "nad3": "TAA",
     "rps10": "TAG", "rpl2": "TAG", "rps19": "TAG", "rps3": "TAA", "rpl16": "TAA", "rpl5": "TAG", "rps14": "TAG", "rps8": "TAA", "rpl6": "TAA", "rps13": "TGA",
      "rps11": "TAA", "rps1": "TAA", "atp8": "TGA", "sdh4": "TAA", "nad4l": "TAA", "tatc": "TAA", "cox2": "TAG", "cox3": "TAG", "nad1": "TAA", "cob": "TGA",
       "rtl": "TAG", "nad9": "TAA", "atp1": "TAA", "cox1": "TAG", "atp4": "TGA", "rps4": "TAA", "atp9": "TAA", "rps2": "TGA", "ccmb": "TGA", "ccmc": "TGA",
        "ccmf": "TAA"},

    "Carica_papaya": {"nad2": "TAA", "nad1": "TAA", "cox2": "TAA", "ccmfc": "CGA", "matr": "TAG", "rps19": "TAA", "ccmc": "CTC", "atp1": "TGA", "nad5": "TAA",
     "nad6": "TAA", "rps4": "TAA", "sdh3": "TGA", "atp6": "CAA", "rps7": "TAA", "rps13": "TGA", "nad7": "TAG", "ccmb": "TGA", "cox1": "TAA", "rps10": "TGA",
      "ccmfn": "TGA", "nad3": "TAA", "rps12": "TGA", "rpl2": "TAA", "rps3": "TAA", "rpl16": "TAA", "atp4": "TAA", "nad4l": "TAA", "atp8": "TGA", "cox3": "TGA",
       "sdh4": "CGA", "cob": "TGA", "rps14": "TAG", "rpl5": "TAA", "nad4": "TGA", "nad9": "TAA", "mttb": "TAG", "atp9": "CGA", "rps1": "CAA"},

    "Nephroselmis_olivacea": {"atp9": "TAA", "atp8": "TAG", "cox2": "TAA", "cox3": "TAA", "atp1": "TAA", "cox1": "TAA", "atp6": "TAA", "atp4": "TAA", "nad4L": "TAA",
     "nad9": "TAA", "cob": "TAA", "nad3": "TAA", "nad2": "TAA", "nad1": "TAA", "nad10": "TAA", "nad7": "TAA", "nad6": "TAG", "nad5": "TAA", "nad4": "TAA"},

    "Physcomitrella_patens": {"cox1": "TAG", "orf622": "TGA", "atp4": "TGA", "atp8": "TGA", "rps1": "TGA", "rps2": "TAG", "ccmb": "TAG", "ccmc": "TAA", "ccmfn": "TAA",
     "ccmf": "TGA", "rps4": "TAA", "tatc": "TAA", "nad4l": "TAA", "sdh4": "TAA", "sdh3": "TAA", "nad5": "TAA", "orf533": "TGA", "nad4": "TAA", "nad2": "TAA",
      "rps12": "TGA", "rps7": "TAA", "atp6": "TAA", "nad6": "TGA", "cox2": "TAA", "cox3": "TAA", "nad1": "TAA", "cob": "TGA", "nad9": "TAA", "atp1": "TAA",
       "nad3": "TAA", "nad7": "TAG", "rpl2": "TAG", "rps19": "TAG", "rps3": "TAA", "rpl16": "TAA", "rpl5": "TAG", "rps14": "TAG", "rpl6": "TAA",
        "rps13": "TGA", "rps11": "TAA", "atp9": "TAA"},

    "Chaetosphaeridium_globosum": {"orf170": "TAA", "orf104": "TGA", "orf202": "TAA", "sdh3": "TAA", "sdh4": "TAA", "nad4l": "TAA", "ymf16": "TAA", "rps4": "TGA", "rps11": "TAA",
     "rps13": "TAA", "rpl6": "TAA", "rps14": "TAA", "rpl5": "TGA", "rpl16": "TGA", "rps3": "TAA", "rps19": "TAA", "rpl2": "TAA", "rps10": "TAA", "nad7": "TAA",
      "atp1": "TAG", "cox3": "TAG", "cox2": "TAA", "orf269": "TAA", "nad6": "TAA", "atp6": "TAA", "rps12": "TAA", "rps7": "TAG", "nad2": "TAA", "nad4": "TGA",
       "nad5": "TAA", "orf126": "TGA", "nad9": "TAA", "nad1": "TAA", "cob": "TGA", "orf207": "TAG", "rps2": "TAA", "rps1": "TAA", "atp8": "TAA", "ymf39": "TAA",
        "cox1": "TAG", "orf303": "TAG", "orf267": "TAG", "atp9": "TAA", "nad3": "TAA"},

    "Mesostigma_viride": {"orf133": "TAA", "orf162": "TGA", "atp9": "TAA", "atp1": "TAA", "atp6": "TAA", "rpl14": "TAA", "atp8": "TAA", "nad7": "TAA", "nad3": "TAA",
     "rpl5": "TAA", "nad6": "TAA", "nad2": "TAA", "ymf39": "TAA", "ymf16": "TAA", "rps1": "TAG", "nad4l": "TAA", "nad1": "TAA", "rps13": "TAA", "rpl6": "TAA",
      "cob": "TAA", "sdh4": "TAA", "rps19": "TAA", "cox2": "TAA", "cox1": "TAA", "cox3": "TAA", "rps4": "TAA", "rps3": "TAA", "sdh3": "TAA", "rps2": "TAA",
       "rps12": "TAA", "rps7": "TAA", "rpl16": "TAA", "nad5": "TAA", "nad4": "TAA", "rps10": "TAA", "rps14": "TAA", "nad9": "TAA", "rps11": "TAA"},

    "Bracteacoccus_aerius": {"cox2a": "TAA", "atp9": "TAA", "cox2": "TAA", "cox3": "TCA", "cox1": "TAA", "atp6": "TAA", "nad3": "TCA", "cob": "TAA", "nad4L": "TCA",
    "nad2": "TCA", "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Pedinomonas_minor gc=4": {"atp8": "TAA", "cox1": "TAA", "atp6": "TAA", "nad4L": "TAA", "cob": "TAA", "nad3": "TAA", "nad2": "TAA", "nad1": "TAA", "nad6": "TAA",
    "nad5": "TAA", "nad4": "TAA"},

    "Kirchneriella_aperta": {"cox2a": "TCA", "atp9": "TCA", "cox2": "TCA", "cox3": "TCA", "cox1": "TCA", "atp6": "TCA", "nad3": "TCA", "cob": "TCA", "nad4L": "TCA",
     "nad2": "TCA", "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Pycnococcus_provasolii": {"rps12": "TTA", "atp9": "TTG", "atp8": "TTG", "cob": "TAA", "cox2": "TTA", "cox3": "TAA", "cox1": "TAG", "atp6": "TAA", "atp4": "TTG",
     "nad3": "TAA", "rps4": "TTA", "rps3": "TTG", "nad4L": "TTA", "nad2": "TTG", "nad1": "TTG", "nad6": "TAA", "nad5": "TAA", "nad4": "TAA"},

    "Chlorotetraedron_incus": {"cox2a": "TCA", "atp9": "TAA", "cox2": "TCA", "cox3": "TAA", "cox1": "TCA", "atp6": "TCA", "nad3": "TAA", "cob": "TCA", "nad4L": "TAA",
     "nad2": "TAA", "nad1": "TAA", "nad6": "TAA", "nad5": "TCA", "nad4": "TAA"},

    "Mychonastes_homosphaera": {"cox2a": "TCA", "atp9": "TCA", "cox2": "TCA", "cox3": "TAA", "cox1": "TCA", "atp6": "TCA", "nad4L": "TCA", "cob": "TCA", "nad3": "TCA",
    "nad2": "TAA", "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Ourococcus_multisporus": {"cox2a": "TCA", "atp9": "TAA", "cox2": "TCA", "cox3": "TCA", "cox1": "TCA", "atp6": "TAA", "nad3": "TCA", "cob": "TAA", "nad4L": "TCA",
     "nad2": "TCA", "nad1": "TCA", "nad6": "TAA", "nad5": "TAA", "nad4": "TCA"},

    "Tetradesmus_obliquus": {"atp9": "TCA", "cox2": "TCA", "cox3": "TCA", "cox1": "TCA", "atp6": "TCA", "nad4L": "TCA", "cob": "TCA", "nad3": "TCA", "nad2": "TCA",
     "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Pseudomuriella_schumacherensis": {"cox2a": "TCA", "atp9": "TCA", "cox2": "TCA", "cox3": "TCA", "cox1": "TCG", "atp6": "TCA", "nad3": "TCG", "cob": "TCA", "nad4L": "TCG",
     "nad2": "TCG", "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Neochloris_aquatica": {"cox2a": "TAA", "atp9": "TCA", "cox2": "TAA", "cox3": "TAA", "cox1": "TAA", "atp6": "TCA", "nad4L": "TAA", "cob": "TCA", "nad3": "TCA",
     "nad2": "TCA", "nad1": "TCA", "nad6": "TAA", "nad5": "TCA", "nad4": "TAA"},

    "Prototheca_wickerhamii": {"atp9": "TAA", "atp8": "TAA", "cox2": "TAA", "cox3": "TAA", "atp1": "TAA", "cox1": "TAA", "atp6": "TAA", "atp4": "TAA", "nad4L": "TAA",
     "nad9": "TAA", "cob": "TAA", "nad3": "TAA", "nad2": "TAA", "nad1": "TAA", "nad7": "TAA", "nad6": "TAA", "nad5": "TAA", "nad4": "TAA"},

    "Bracteacoccus_minor": {"cox2a": "TCA", "atp9": "TAA", "cox2": "TCA", "cox3": "TCA", "cox1": "TCA", "atp6": "TCA", "nad4L": "TCA", "cob": "TCA", "nad3": "TCA",
     "nad2": "TCA", "nad1": "TCA", "nad6": "TCA", "nad5": "TCA", "nad4": "TCA"},

    "Pseudendoclonium_akinetum": {"atp9": "TAA", "atp8": "TAA", "cox2": "TAA", "cox3": "TAA", "atp1": "TAA", "cox1": "TAA", "atp6": "TAA", "atp4": "TAA", "nad3": "TAG",
    "cob": "TAA", "nad4L": "TAG", "nad2": "TAG", "nad1": "TAA", "nad7": "TAA", "nad6": "TAA", "nad5": "TAG", "nad4": "TAA"},

    "Chromochloris_zofingiensis": {"cox2a": "TCA", "atp9": "TAA", "cox2": "TCA", "cox3": "TCA", "cox1": "TAA", "atp6": "TAA", "nad4L": "TAA", "cob": "TAA", "nad3": "TAA",
     "nad2": "TCA", "nad1": "TCA", "nad6": "TCA", "nad5": "TAA", "nad4": "TCA"}}
    return stop_mapper
