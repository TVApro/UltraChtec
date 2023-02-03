# все функции для работы с fasta-файлами
import re
import textwrap as tw

def privet():
        print("Привет, я работаю! Я существую")
        
def fasta_obtain(text):
        from re import sub
        fasta_sorted_list = {}
        fasta_list = sorted(text.split(">"), key=len, reverse=True)
        for i in fasta_list:
                if i != ' ' and i != '':
                        end = i.find("\n")
                        fasta_name = '>'+sub('\n', '', i[:end])
                        fasta_sequence = sub('\n', '', i[end:])
                        fasta_sorted_list[fasta_name]=fasta_sequence
        return fasta_sorted_list

def full_fasta_sequence(fasta_sorted_list):
        full_fasta_sequence = str("")
        for i in fasta_sorted_list:
                try:
                        full_fasta_sequence += fasta_sorted_list[i]
                except:
                        full_fasta_sequence += i
        return full_fasta_sequence

def simple_sequence(seq):
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('a', 'A').replace('c', 'C').replace('t', 'T').replace('g', 'G').replace('n', 'N')
        return seq

def beautifull_sequence(seq):
        seq = tw.fill(seq, width=60)
        return seq
        
def complement_sequence(seq):
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('A', '%1%').replace('a', '%1%').replace('T', 'A').replace('t', 'A').replace('%1%', 'T')
        seq = seq.replace('C', '%2%').replace('c', '%2%').replace('G', 'C').replace('g', 'C').replace('%2%', 'G')
        seq = seq[::-1]
        return seq
        
def nucleotides_counter(seq):
        seq = simple_sequence(seq)
        A = seq.count('A')
        G = seq.count('G')
        T = seq.count('T')
        C = seq.count('C')
        N = seq.count('N')
        GC = (G+C)*100/(G+C+A+T)
        perA = A*100/(A+G+T+C+N)
        perG = G*100/(A+G+T+C+N)
        perT = T*100/(A+G+T+C+N)
        perC = C*100/(A+G+T+C+N)
        perN = N*100/(A+G+T+C+N)
        perA = '('+str(format(perA, '.2f'))+'%'+')'
        perG = '('+str(format(perG, '.2f'))+'%'+')'
        perT = '('+str(format(perT, '.2f'))+'%'+')'
        perC = '('+str(format(perC, '.2f'))+'%'+')'
        perN = '('+str(format(perN, '.2f'))+'%'+')'
        perGC = str(format(GC, '.2f'))+'%'
        lenth = "Длина последовательности = {len} bp \n".format(len=len(seq))
        countx = "{lenth}Количество оснований: \n A {A} bp {perA} \n G {G} bp {perG} \n T {T} bp {perT} \n C {C} bp {perC}, \n N {N} bp {perN} \n GC состав {perGC}".format(lenth=lenth, A=A, G=G, T=T, C=C, N=N, perA=perA, perG=perG, perT=perT, perC=perC, perN=perN, perGC=perGC)
        return(countx)

def N50_L50(fasta_sorted_list):
        full_sequence = full_fasta_sequence(fasta_sorted_list)
        N50 = 0
        L50 = 0
        half_seq = len(full_sequence)/2
        medium = 0
        try:
                for i in fasta_sorted_list:
                        medium += len(fasta_sorted_list[i])
                        L50 += 1
                        if medium >= half_seq:
                                N50 = len(fasta_sorted_list[i])
                                break
        # это покажется бессмысленным, но у списка контигов, получаемого дальше, нет имён,
        # поэтому первый алгоритм не может с ними работать, и требуется тот, что ниже
        except: 
                for i in fasta_sorted_list:
                        medium += len(i)
                        L50 += 1
                        if medium >= half_seq:
                                N50 = len(i)
                                break
        N50L50 = [N50, L50]
        return N50L50

def scaffold_devider(fasta_sorted_list):
        # разделяет скаффолды на контиги, считает количество и тех, и других
        contig_count = 0
        scaffold_count = 0
        contig_list = []
        for i in fasta_sorted_list:
                scaffold_count += 1
                seq = re.sub(r'NNNNNNNNN+', '$', fasta_sorted_list[i])
                contigs = seq.split('$')
                for g in contigs:
                        if g:
                                contig_count += 1
                                contig_list.append(g)
        contig_list.sort(key=lambda x: len(x), reverse=True)
        output_list = [contig_count, scaffold_count, contig_list]
        return output_list

def max_len(fasta_sorted_list):
        # вычисляет максимальную длину контига/скаффолда из данного списка
        max_seq_len = 0
        num = 0
        for i in fasta_sorted_list:
            num += 1
            if num == 1:
                try:
                    max_seq_len = len(fasta_sorted_list[i])
                except:
                    max_seq_len = len(i)
        return max_seq_len

def fasta_resolve(fasta_sorted_list):
        full_seq = str('')
        # получение количества контигов, скаффолдов и списка контигов
        output_sequence_count_list = scaffold_devider(fasta_sorted_list)
        contig_number = output_sequence_count_list[0]
        scaffold_number = output_sequence_count_list[1]
        contig_list = output_sequence_count_list[2]
        # Получение полной последовательности
        full_seq = full_fasta_sequence(fasta_sorted_list)
        # Получение максимальной длины скаффолда
        max_scaffold_len = max_len(fasta_sorted_list)
        # Получение максимальной длины контига
        max_contig_len = max_len(contig_list)
        # Получение значений N50/L50 скаффолдов
        N50_L50_scaffold = N50_L50(fasta_sorted_list)
        N50_scaffold = N50_L50_scaffold[0]
        L50_scaffold = N50_L50_scaffold[1]
        # Получение значений N50/L50 контигов
        N50_L50_contig = N50_L50(contig_list)
        N50_contig = N50_L50_contig[0]
        L50_contig = N50_L50_contig[1]
        # Вывод
        assembly_properties = ("\nПоказатели сборки:\n количество скаффолдов - {scaffold_number} шт.\n максимальная длина скаффолда {max_scaffold_len} bp \n N50 скаффолдов {N50_scaffold} bp \n L50 скаффолдов {L50_scaffold} шт.\n количество контигов - {contig_number} шт.\n максимальная длина контига {max_contig_len} bp \n N50 контигов {N50_contig} bp \n L50 контигов {L50_contig} шт.\n".format(scaffold_number=scaffold_number, \
                                max_scaffold_len=max_scaffold_len, N50_scaffold=N50_scaffold, L50_scaffold=L50_scaffold, \
                                contig_number=contig_number, max_contig_len=max_contig_len, N50_contig=N50_contig, \
                                L50_contig=L50_contig))
        return assembly_properties
