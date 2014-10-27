import argparse
from Bio import SearchIO

# Пример как запускать:
#python3 gene-check.py -sff //home/anonym/PycharmProjects/untitled/pseudomolecules.fa -saf //home/anonym/PycharmProjects/untitled/ex_out -pn exonerate -m 1 


# DEBUG_ выставляйте в True при дебаге
DEBUG_ = 0

codon_start_dna = 'ATG'
codon_stop_dna  = ['TAA','TAG','TGA']

def isStopCodon(codon):
    assert(len(codon) == 3)

    if codon in codon_stop_dna:
        return True

    return False

def isStartCodon(codon):
    assert(len(codon) == 3)

    if codon == codon_start_dna:
        return True

    return False

def haveStopCodonInSeq(seq):
    assert ((len(seq) % 3) == 0)

    for i in range(0, len(seq), 3):
        if isStopCodon(seq[i : i + 3]):
            return True

    return False

"""
Protein_Entries is class that contains name of protein and his hsps.
HSP is class that representing high-scoring region(s) between query and hit.
full documentation of class:
http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp.HSP-class.html
and
http://biopython.org/DIST/docs/api/Bio.SearchIO.ExonerateIO-module.html
"""

class Protein_Entries:

    def __init__(self, name):
        self.name = name
        self.hsps = []

    def __hash__(self):
        return self.name.__hash__()

    def addHSP(self, hsp):
        self.hsps.append(hsp)

def GetWrappedEntries(exn_file, format):
    """ Wrap BioPython exonerate SearchIO because not enough clear
     struncture of response.
     First argument is file path to exonerate file and second is format of
     output.
     Return is list of Protein_Entries"""

    qresult = SearchIO.parse(exn_file, format)
    all_results = list(qresult)

    proteins_entries = []

    for query_result in all_results:
        for hit in query_result.hits:
            for hsp in hit.hsps:
                prot_name = hsp.query_id
                protein_entry = [entry for entry in proteins_entries
                if entry.name == prot_name]
                if protein_entry == []:
                    proteins_entries.append(Protein_Entries(prot_name))
                    protein_entry = [entry for entry in proteins_entries
                    if entry.name == prot_name]
                protein_entry = protein_entry[0]
                protein_entry.addHSP(hsp)
    """Example of how proteins_entries looks like:
    Protein  p03  hsps:  [HSP(hit_id='chr1', query_id='p03', 5 fragments),
    HSP(hit_id='chr1', query_id='p03', 2 fragments),
    HSP(hit_id='chr2', query_id='p03', 5 fragments),
    HSP(hit_id='chr2', query_id='p03', 1 fragments)]


    Protein  p01  hsps:  [HSP(hit_id='chr1', query_id='p01', 2 fragments),
    HSP(hit_id='chr1', query_id='p01', 2 fragments)]

    'fragments' mean, that introns between them"""
    return proteins_entries

# Для выода диагностических сообщейний при дебаге
def printdbg(*args):
    if DEBUG_:
        print(*args)

# Для нашей работы нам наужно сам геном и файл выравнивания с белками
def TotalCheck(_genome, _protein_alignment, criterion = 0):

    def CheckORFAndEIS(genome, protein_alignment):
        print("Step 3.1 - CheckOpenReadingFrame And ExoneIntroneStructure")
        protein_alignment_result = []

        for protein in protein_alignment:
            print("--------------------------------------")
            print("Сейчас проверяется белок", protein.name)
            # Проверяем кратность 3м для найденного выравнивания
            # А также расположение старт и стоп кодонов
            for hsp_ in protein.hsps:
                assert (hsp_.hit_end - hsp_.hit_start > 0)
                if isStartCodon(genome[hsp_.hit_id][hsp_.hit_start: hsp_.hit_start + 3]):
                    printdbg("-Старт кодон на месте", [protein.name, hsp_.hit_id, hsp_.hit_start, hsp_.hit_end,
                                                      genome[hsp_.hit_id][hsp_.hit_start: hsp_.hit_start + 3]])

                    if isStopCodon(genome[hsp_.hit_id][hsp_.hit_end - 3: hsp_.hit_end]):
                        printdbg("--Стоп кодон на месте", [protein.name, hsp_.hit_id, hsp_.hit_start, hsp_.hit_end,
                                                      genome[hsp_.hit_id][hsp_.hit_end - 3: hsp_.hit_end]])
                        # Склеиваем экзоны
                        r = []
                        r.append(hsp_.hit_start)
                        for rng in hsp_.hit_inter_ranges:
                            r.append(rng[0])
                            r.append(rng[1])
                        r.append(hsp_.hit_end - 3)

                        buf = ""
                        for i in range(len(r) - 1):
                            if (i + 1) % 2:
                                buf += genome[hsp_.hit_id][r[i] : r[i + 1]]

                        if len(buf) % 3 == 0:
                            printdbg("---Последовательность между старт и стоп кодоном кратна трем")

                            if not haveStopCodonInSeq(buf):
                                printdbg("----Последовательность не имеет стоп кодонов посреди экзонов")
                            else:
                                break

                            protein_alignment_result.append(Protein_Entries(protein.name))
                            protein_alignment_result[len(protein_alignment_result) - 1].addHSP(hsp_)
                            print("Результат утвержден!")
                            print("Хромосома:", hsp_.hit_id)
                            print("Белок:", protein.name)
                            print("Позиция начала:", hsp_.hit_start)
                            print("Вероятный сайт связывания(для дальнейшей проверки):",
                                  genome[hsp_.hit_id][hsp_.hit_start - 6 : hsp_.hit_start + 4])
                            print("Интроны:", hsp_.hit_inter_ranges)
                            print("Позиция конца:", hsp_.hit_end)
                            print("--------------------------------------")

                        else:
                            print(buf)
                            printdbg("Последовательность между старт и стоп кодоном не кратна трем", len(r) % 3)
                    else:
                        printdbg("Стоп кодон не на месте", [protein.name, hsp_.hit_id, hsp_.hit_start, hsp_.hit_end,
                                                      genome[hsp_.hit_id][hsp_.hit_end - 3: hsp_.hit_end]])
                else:
                    printdbg("Старт кодон не на месте", [protein.name, hsp_.hit_id, hsp_.hit_start, hsp_.hit_end,
                                                      genome[hsp_.hit_id][hsp_.hit_start: hsp_.hit_start + 3]])

        return protein_alignment_result

    def CheckCodonUsageBias(genome, protein_alignment):
        print("Step 3.2 - CheckCodonUsageBias; Not working now")
        return []

    check_stack = [CheckORFAndEIS,
                   CheckCodonUsageBias]

    for i in range(len(check_stack)):
        if criterion == 0 or criterion == i + 1:
            _protein_alignment = check_stack[i](_genome, _protein_alignment)

    return _protein_alignment

def ReadProteinAlignment(filename, file_from = "exonerate"):
    if filename == "":
        Exception("Set the file name for check")

    file_from = file_from.lower()

    result = []

    # В разработке, поэтому пишем парсер спрева только для exonerate'a
    if file_from == "exonerate":
        result = GetWrappedEntries(filename, "exonerate-vulgar")
    else:
        Exception("Now available only for the exonerate")

    return result

def ReadFasta(source_fasta_file_name):
    genome_file = open(source_fasta_file_name, "r")
    genome_sorce_fasta = genome_file.read()
    genome_file.close()

    genome = {}
    pos1 = 0
    pos2 = 0
    while True:
        pos1 = genome_sorce_fasta.find(">", pos2)
        if pos1 < 0:
            break
        pos2 += 1
        pos2 = genome_sorce_fasta.find(">", pos2)
        if pos2 < 0:
            chrom_name = genome_sorce_fasta[pos1 + 1: pos1 + 6].replace("\n", "")
            genome[chrom_name] = genome_sorce_fasta[pos1 + 6 : ].replace("\n", "")
            break

        chrom_name = genome_sorce_fasta[pos1 + 1: pos1 + 6].replace("\n", "")
        genome[chrom_name] = genome_sorce_fasta[pos1 + 6 : pos2 - 1].replace("\n", "")

    return genome

####################    Основная логика    ####################
def __main__():
    # Парсим аргументы командной строки и вообще ведем себя как порядочная утилита
    parser = argparse.ArgumentParser()

    parser.add_argument('-sff', required = True,  help = 'Path to source fasta file')
    parser.add_argument('-saf', required = True, help = 'Path to source alignment file')
    parser.add_argument('-m', required = False, choices = ['0', '1', '2'], help = 'Check mode: 0 - All(default), '
                                                '1 - open read frame and exone-intron structure, 2 - codon usage bias')
    parser.add_argument('-pn', required = True, choices = ["exonerate"] , help = 'Program name what make alignment')
    parser.add_argument('-d', required = False, help = 'Destination folder for output')

    parse_result = vars(parser.parse_args())

    printdbg("Your input:", parse_result)

    source_fasta   = parse_result["sff"]
    alignment_file = parse_result["saf"]
    program_name   = parse_result["pn"]
    destination    = parse_result["d"] if parse_result["d"] else "result_of_gene-check.fasta"
    mode           = int(parse_result["m"]) if parse_result["m"] else 0

    print("Real destination for output:", destination)

    # source_fasta = "pseudomolecules.fa"
    # alignment_file = "ex_out"
    # program_name = "exonerate"
    # destination = "result_of_gene-check.fasta"
    # mode = 1

    print("Step 1 - Reading protein alignment")
    common_protein_alignment_data = ReadProteinAlignment(alignment_file, program_name)

    print("Step 2 - Reading source fasta")
    genome = ReadFasta(source_fasta)

    print("Step 3 - Gene-check")
    result = TotalCheck(genome, common_protein_alignment_data, mode)

    print()
    print("Step 4 - Save result")
    print("------------------ Результат ---------------------")
    if len(result) == 0:
        print("Ни один найденый белок не прошел заданные критерии")
    else:
        output = open(destination, "w")
        for protein in result:
            printdbg(protein.name)
            for hsp_ in protein.hsps:
                assert (hsp_.hit_end - hsp_.hit_start > 0)

                output.write(str("> " + protein.name + " " + hsp_.hit_id + " " +
                                 genome[hsp_.hit_id][hsp_.hit_start - 6 : hsp_.hit_start + 4] + " "
                                 + str(hsp_.hit_start) + " "
                                 + str(hsp_.hit_end)))

                output.write("\n")
                output.write(genome[hsp_.hit_id][hsp_.hit_start: hsp_.hit_end])
                output.write("\n")
                print(hsp_.hit_id, [hsp_.hit_start, hsp_.hit_end])


        output.close()

# Запуск
__main__()