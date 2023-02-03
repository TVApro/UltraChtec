# всё для работы с последовательностями genbank
import re
import fastaUCHv1 as fa

def gb_obtain(text):
    assembly_info = {}
    LOCUS = re.compile(r'LOCUS   .*')
    locus_list = re.findall(LOCUS, text)
    for i in locus_list:
        n = locus_list.index(i)
        i = re.sub(r'LOCUS | DNA linear.*', '',re.sub(r' +|\\n', ' ', i))
        locus_list[n] = i
    text_list = re.split(LOCUS, text)
    annotations_list = []
    fasta_list = []
    for i in text_list:
        if i != '' and i != " ":
            ind_FEATURES = i.find("FEATURES ")
            ind_ORIGIN = i.find("ORIGIN")
            ind_end = i.find("//", ind_ORIGIN)
            annotations = i[ind_FEATURES:ind_ORIGIN]
            seq = fa.simple_sequence(re.sub(r'\d*', '', i[ind_ORIGIN+6:ind_end]))
            fasta_list.append(seq)
            annotations_list.append(annotations)
    # создаём связку имена контигов/последовательность фаста
    fa_list = dict(zip(locus_list, fasta_list))
    # создаём связку имена контигов/записи аннотаций
    an_list = dict(zip(locus_list, annotations_list))
    # сортируем контиги по длине
    sorted_values = sorted(fa_list.values(), key=lambda x:len(x), reverse=True)
    sort_fa_list = {}
    sort_an_list = {}
    # здесь встала такая проблема: не мог сортировать контиги и при этом оставлять им их имена, поэтому
    # решил проблему вот так, несколько неуклюже
    for i in sorted_values:
        for k in fa_list:
            if fa_list[k] == i:
                sort_fa_list[k] = i
    # синхронизировал сортированный список контигов и список аннотаций
    for i in sort_fa_list:
        for k in an_list:
            if i == k:
                sort_an_list[i] = an_list[k]
    # поскольку обычно требуются оба списка, то решил не заморачиваться и просто объединил их в один
    ogogo = [sort_fa_list, sort_an_list]
    return ogogo

def all_genes_names(text):
    # получение всех типов аннотаций, которые есть в геноме
    list1 = re.findall(r' {5}\S+  +\d+\W+\d+| {5}.+  +complement\(\d+\W+\d+\)', text)
    title = re.compile(r' {5}\S+  +')
    titles = set()
    for i in list1:
        yes = re.search(title, i)
        titles.add(yes.group())
    an_types = '\nПриведённые имена аннотаций:\n'
    for i in titles:
        d = re.sub(' +|\n', '', i)
        an_types = an_types+' '+d+'\n'
    return an_types

def coding_koeff(genome):
    fasta_list = []
    gene_seq_len = 0
    seq_len = 0
    LOCUS = re.compile(r'LOCUS\s\s+.*')
    scaffold_list = re.split(LOCUS, genome)
    for i in scaffold_list:
                    if i:
                            ind_FEATURES = i.find("FEATURES ")
                            ind_ORIGIN = i.find("ORIGIN")
                            ind_end = i.find("//", ind_ORIGIN)
                            annotations = i[ind_FEATURES:ind_ORIGIN]
                            seq = re.sub(r'\d*', '', i[ind_ORIGIN+6:ind_end]).replace(' ', '').replace('\n', '')
                            seq_len += len(seq)
                            annotations = i[ind_FEATURES:ind_ORIGIN]
                            title = re.compile(r' {5}\S+  +')
                            edges = re.compile(r'\d+\W+\d+|complement\(\d+\W+\d+?\)')
                            assembly_split = re.split(title, annotations)
                            for g in assembly_split:
                                    if g:
                                            if '/product=' in g:
                                                    gene_edges = re.search(edges, g)
                                                    gg_edges = gene_edges.group()
                                                    numbers = re.findall('\d+', gg_edges)
                                                    u = 0
                                                    for h in numbers:
                                                            u += 1
                                                            if u == 1:
                                                                    st = int(h)
                                                            if u == 2:
                                                                    fin = int(h)
                                                                    delta_l = abs(fin-st)
                                                    gene_seq_len += delta_l
    code_koef = '\n коэффициент кодирования - {0:.2f}% \n'.format((gene_seq_len/seq_len)*100)
    return code_koef

def features_counter(text):
    # просто считает имена наиболее распространённых видов аннотаций
    gene_count = text.count("   gene     ")
    CDS_count = text.count("   CDS   ")
    misc_feature_count = text.count("   misc_feature  ")
    tRNA_count = text.count("   tRNA   ")
    rRNA_count = text.count("   rRNA   ")
    ncRNA_count = text.count("   ncRNA   ")
    repeat_region_count = text.count("   repeat_region  ")
    RNA_count = text.count("   RNA   ")
    mRNA_count = text.count("     mRNA      ")
    rna_count = text.count("     rna        ")
    crispr_array = text.count("     crispr_array    ")
    crispr_repeat = text.count("     crispr_repeat   ")
    crispr_spacer = text.count("     crispr_spacer   ")
    prophage_count = text.count("     prophage       ")
    repeat_count = text.count("     repeat          ")
    assembly_gap_count = text.count("     assembly_gap    ")
    print1 = str('')
    print2 = str('')
    print3 = str('')
    print4 = str('')
    print5 = str('')
    print6 = str('')
    print7 = str('')
    print8 = str('')
    print9 = str('')
    print10 = str('')
    print11 = str('')
    print12 = str('')
    print13 = str('')
    if gene_count > 0:
        print1 = " указано генов - {gene_count} шт.\n".format(gene_count=gene_count)
    if CDS_count > 0:
        print2 = " СDS - {CDS_count} шт.\n".format(CDS_count=CDS_count)
    if RNA_count > 0:
        print3 = " РНК - {RNA_count} шт.\n".format(RNA_count=RNA_count)
    if rna_count > 0:
        print4 = " РНК - {rna_count} шт.\n".format(rna_count=rna_count)
    if tRNA_count > 0:
        print5 = " тРНК - {tRNA_count} шт.\n".format(tRNA_count=tRNA_count)
    if rRNA_count > 0:
        print6 = " рРНК - {rRNA_count} шт.\n".format(rRNA_count=rRNA_count)
    if mRNA_count > 0:
        print7 = " мРНК - {mRNA_count} шт.\n".format(mRNA_count=mRNA_count)
    if ncRNA_count > 0:
        print8 = " нкРНК - {ncRNA_count} шт.\n".format(ncRNA_count=ncRNA_count)
    if misc_feature_count > 0:
        print9 = " прочие признаки - {misc_feature_count} шт.\n".format(misc_feature_count=misc_feature_count)
    if repeat_region_count > 0 or repeat_count > 0:
        repeats = repeat_region_count + repeat_count
        print10 = " повторов - {repeats} шт.\n".format(repeats=repeats)
    if crispr_array > 0 or crispr_repeat > 0 or crispr_spacer > 0:
        print11 = " массивов CRISPR - {crispr_array}, CRISPR-повторов - {crispr_repeat}, CRISPR-спейсеров - {crispr_spacer}\n".format(crispr_array=crispr_array, crispr_repeat=crispr_repeat, crispr_spacer=crispr_spacer)
    if prophage_count > 0:
        print12 = " профагов - {prophage_count} шт.\n".format(prophage_count=prophage_count)
    if assembly_gap_count > 0:
        print13 = " пропусков в последовательности - {gap} шт.\n".format(gap=assembly_gap_count)
    features_output = print1+print2+print3+print4+print5+print6+print7+print8+print9+print10+print11+print12+print13
    return features_output


