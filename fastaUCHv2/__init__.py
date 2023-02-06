# all functions for working with FASTA files
import re
import textwrap as tw
        
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
        lenth = "Sequence length = {len} bp \n".format(len=len(seq))
        countx = "{lenth}Number of bases: \n A {A} bp {perA} \n G {G} bp {perG} \n T {T} bp {perT} \n C {C} bp {perC}, \n N {N} bp {perN} \n GC composition {perGC}".format(lenth=lenth, A=A, G=G, T=T, C=C, N=N, perA=perA, perG=perG, perT=perT, perC=perC, perN=perN, perGC=perGC)
        return(countx)

def N50_L50(fasta_sorted_list):
        full_sequence = full_fasta_sequence(fasta_sorted_list)
        N50 = 0
        N90 = 0
        L50 = 0
        L90 = 0
        summ = 0
        half_seq = len(full_sequence)/2
        more_seq = (len(full_sequence)/100)*90
        try:
                for i in fasta_sorted_list:
                        summ += len(fasta_sorted_list[i])
                        L50 += 1
                        if summ >= half_seq:
                                N50 = len(fasta_sorted_list[i])
                                summ = 0
                                break
                for i in fasta_sorted_list:
                        summ += len(fasta_sorted_list[i])
                        L90 += 1
                        if summ >= more_seq:
                                N90 = len(fasta_sorted_list[i])
                                break
        # it will seem meaningless, but the list of contigs received further has no names,
        # so the first algorithm cannot work with them and the one below is required
        except: 
                for i in fasta_sorted_list:
                        summ += len(i)
                        L50 += 1
                        if summ >= half_seq:
                                N50 = len(i)
                                summ = 0
                                break
                for i in fasta_sorted_list:
                        summ += len(i)
                        L90 += 1
                        if summ >= more_seq:
                                N90 = len(i)
                                break
        N50L50 = [N50, L50, N90, L90]
        return N50L50

def scaffold_devider(fasta_sorted_list):
        # divides scaffolds into contigs, counts the number of both
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
        # calculates the maximum length of a contig/scaffold from the given list
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
        # getting the number of contigs, scaffolds and a list of contigs
        output_sequence_count_list = scaffold_devider(fasta_sorted_list)
        contig_number = output_sequence_count_list[0]
        scaffold_number = output_sequence_count_list[1]
        contig_list = output_sequence_count_list[2]
        # getting the full sequence
        full_seq = full_fasta_sequence(fasta_sorted_list)
        # getting the maximum length of the scaffold
        max_scaffold_len = max_len(fasta_sorted_list)
        # getting the maximum length of a contig
        max_contig_len = max_len(contig_list)
        # getting N50/L50 values of scaffolds
        N50_L50_scaffold = N50_L50(fasta_sorted_list)
        N50_scaffold = N50_L50_scaffold[0]
        L50_scaffold = N50_L50_scaffold[1]
        N90_scaffold = N50_L50_scaffold[2]
        L90_scaffold = N50_L50_scaffold[3]
        # getting N50/L50 values of scaffolds
        N50_L50_contig = N50_L50(contig_list)
        N50_contig = N50_L50_contig[0]
        L50_contig = N50_L50_contig[1]
        N90_contig = N50_L50_contig[2]
        L90_contig = N50_L50_contig[3]
        # conclusion
        assembly_properties = ("\nAssembly metrics:\n number of scaffolds - {scaffold_number} \n maximum scaffold length {max_scaffold_len} bp \n scaffolds N50 {N50_scaffold} bp \n scaffolds L50 {L50_scaffold} \n scaffolds N90 {N90_scaffold} bp\n scaffolds L90 {L90_scaffold} \n number of contigs - {contig_number}\n maximum contig length {max_contig_len} bp\n contigs N50 {N50_contig} bp\n contigs L50 {L50_contig} \n contigs N90 {N90_contig} bp \n contigs L90 {L90_contig} \n".format(scaffold_number=scaffold_number, \
                                max_scaffold_len=max_scaffold_len, N50_scaffold=N50_scaffold, L50_scaffold=L50_scaffold, \
                                contig_number=contig_number, max_contig_len=max_contig_len, N50_contig=N50_contig, \
                                L50_contig=L50_contig, N90_scaffold=N90_scaffold, L90_scaffold=L90_scaffold,\
                                N90_contig=N90_contig, L90_contig=L90_contig))
        return assembly_properties
