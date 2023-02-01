#! /usr/bin/python3
# Не судите строго
from tkinter import *
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog as fd

import textwrap as tw
import re
import os
import sys
import datetime
import platform

class fasta():
    def fasta_obtain(text):
        from re import sub
        fa_list = {}
        fasta_list = sorted(text.split(">"), key=len, reverse=True)
        for i in fasta_list:
            if i != ' ' and i != '':
                end = i.find("\n")
                fasta_name = '>'+sub('\n', '', i[:end])
                fasta_sequence = sub('\n', '', i[end:])
                fa_list[fasta_name]=fasta_sequence
        return(fa_list)
    def full_seq(fa_list):
        full_seq = str("")
        for i in fa_list:
            try:
                full_seq += fa_list[i]
            except:
                full_seq += i
        return full_seq
    def simple_sequence(seq):
        import textwrap as tw
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('a', 'A').replace('c', 'C').replace('t', 'T').replace('g', 'G').replace('n', 'N')
        return(seq)
    def beautifull_sequence(seq):
        seq = tw.fill(seq, width=60)
        return seq
    def complement_sequence(seq):
        seq = seq.replace(' ', '').replace('\n', '').replace('-', '')
        seq = seq.replace('A', '%1%').replace('a', '%1%').replace('T', 'A').replace('t', 'A').replace('%1%', 'T')
        seq = seq.replace('C', '%2%').replace('c', '%2%').replace('G', 'C').replace('g', 'C').replace('%2%', 'G')
        seq = seq[::-1]
        return(seq)
    def nucleotides_counter(text):
        text = text.replace('a', 'A').replace('t', 'T').replace('c', 'C').replace('g', 'G').replace('n', 'N')
        A = text.count('A')
        G = text.count('G')
        T = text.count('T')
        C = text.count('C')
        N = text.count('N')
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
        lenth = "Длина последовательности = {len} bp \n".format(len=len(text))
        countx = "{lenth}Количество оснований: \n A {A} bp {perA} \n G {G} bp {perG} \n T {T} bp {perT} \n C {C} bp {perC}, \n N {N} bp {perN} \n GC состав {perGC}".format(lenth=lenth, A=A, G=G, T=T, C=C, N=N, perA=perA, perG=perG, perT=perT, perC=perC, perN=perN, perGC=perGC)
        return(countx)
    def N50_L50(seq_list):
        full_seq = fasta.full_seq(seq_list)
        N50 = 0
        L50 = 0
        half_seq = len(full_seq)/2
        medium = 0
        try:
            for i in seq_list:
                medium += len(seq_list[i])
                L50 += 1
                if medium >= half_seq:
                    N50 = len(seq_list[i])
                    break
        except:
            for i in seq_list:
                medium += len(i)
                L50 += 1
                if medium >= half_seq:
                    N50 = len(i)
                    break
        N50L50 = [N50, L50]
        return N50L50
    def scaffold_devider(fa_list):
        contig_number = 0
        scaffold_number = 0
        contig_list = []
        for i in fa_list:
            scaffold_number += 1
            seq = re.sub(r'NNNNNNNNN+', '$', fa_list[i])
            contigs = seq.split('$')
            for g in contigs:
                if g != '' and g != " ":
                    contig_number += 1
                    contig_list.append(g)
        contig_list.sort(key=lambda x: len(x), reverse=True)
        scaffold_list = [contig_number, scaffold_number, contig_list]
        return scaffold_list
    def max_len(fa_list):
        max_seq_len = 0
        num = 0
        for i in fa_list:
            num += 1
            if num == 1:
                try:
                    max_seq_len = len(fa_list[i])
                except:
                    max_seq_len = len(i)
        return max_seq_len
    def fasta_resolve(fa_list):
        full_seq = str('')
        full_seq = fasta.full_seq(fa_list)
        max_scaffold_len = fasta.max_len(fa_list)
        scaffold_list = fasta.scaffold_devider(fa_list)
        contig_number = scaffold_list[0]
        scaffold_number = scaffold_list[1]
        contig_list = scaffold_list[2]
        max_contig_len = fasta.max_len(contig_list)
        N50_L50_scaffold = fasta.N50_L50(fa_list)
        N50_scaffold = N50_L50_scaffold[0]
        L50_scaffold = N50_L50_scaffold[1]
        N50_L50_contig = fasta.N50_L50(contig_list)
        N50_contig = N50_L50_contig[0]
        L50_contig = N50_L50_contig[1]
        assembly_properties = ("\nПоказатели сборки:\n количество скаффолдов - {scaffold_number} шт.\n максимальная длина скаффолда {max_scaffold_len} bp \n N50 скаффолдов {N50_scaffold} bp \n L50 скаффолдов {L50_scaffold} шт.\n количество контигов - {contig_number} шт.\n максимальная длина контига {max_contig_len} bp \n N50 контигов {N50_contig} bp \n L50 контигов {L50_contig} шт.\n".format(scaffold_number=scaffold_number, \
                                max_scaffold_len=max_scaffold_len, N50_scaffold=N50_scaffold, L50_scaffold=L50_scaffold, \
                                contig_number=contig_number, max_contig_len=max_contig_len, N50_contig=N50_contig, \
                                L50_contig=L50_contig))
        return assembly_properties
class genbank():
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
                seq = fasta.simple_sequence(re.sub(r'\d*', '', i[ind_ORIGIN+6:ind_end]))
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
        an_types = '\nВиды аннотаций в геноме:\n'
        for i in titles:
            d = re.sub(' +|\n', '', i)
            an_types = an_types+' '+d+'\n'
        return an_types
    def features_counter(text):
        # просто считает имена наиболее распространённых видов аннота
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
class StartWindow(Frame):
    def path_input():
        global pathX
        platformX = platform.system()
        pathX = fd.askopenfilename(title='Выбери файл для открытия',
                                filetypes=[('All Files', '*'),
                                           ('fasta', '*.fa'),
                                           ('fasta', '*fas'),
                                           ('fasta', '*.fasta'),
                                           ('fasta', '*.fna'),
                                           ('fasta', '*faa'),
                                           ('genbank', '*.gb',),
                                           ('genbank', '*.gbk'),
                                           ('genbank', '*.gbff')])
        if pathX == '':
            return
        else:
            x = path_entry.get()
            if x != '' and x !=' ':
                path_entry.delete(0, END)
        if platformX == 'Windows':
            pathX = pathX.replace("/","\\")
        path_entry.insert(INSERT, pathX)
    def file_open():
        global path_file
        path_file = path_entry.get()
        with open(path_file, 'r') as file:
            text = file.read()
        return text
    def good_choice(text):
        gene_count = text.count("   gene     ")
        mRNA_count = text.count("     mRNA       ")
        if text[:1] == ">":
            file_format = 'FASTA'
        elif text[:5] == 'LOCUS' and gene_count > 0 and mRNA_count == 0:
            file_format = 'GB'
        elif text[:5] == 'LOCUS' and gene_count == 0:
            file_format = 'GBK'
        elif text[:5] == 'LOCUS' and gene_count > 0 and mRNA_count > 0:
            file_format = 'GBFF'
        else:
            file_format = 'bad'
        return file_format
    def input_assembly():
        try:
            global path_file
            file_format = str('')
            file_text = StartWindow.file_open()
            time_now = datetime.datetime.now()
            time = str(time_now.strftime("%d-%m-%Y %H:%M"))
            file_format = StartWindow.good_choice(file_text)
        except:
            messagebox.showinfo('ОШИБКА', 'Файл не найден')
        if file_format == 'bad':
            messagebox.showinfo('ОШИБКА', 'Формат для анализа неподдерживаемый')
        else:
            check = islogfile.get()
            if check == 1:
                pass
            else:
                StartWindow.clear()
        if file_format == 'FASTA':
            message = path_file+'\nФайл формата fasta\n'+time
        elif file_format == 'GB':
            message = path_file+'\nФайл формата genbank\n'+time
        elif file_format == 'GBK':
            message = path_file+'\nФайл формата .gbk\n'+time
        elif file_format == 'GBFF':
            message = path_file+'\nФайл формата .gbff\n'+time
        features_list.insert(INSERT, message+'\n')
        if file_format == 'FASTA':
            fa_list = fasta.fasta_obtain(file_text)
            full_assembly = str('')
            for i in fa_list:
                features_list.insert(INSERT, '\n'+i)
                nucl = fasta.nucleotides_counter(fa_list[i])
                features_list.insert(INSERT, '\n'+nucl+'\n')
            features_list.insert(INSERT, "_____________________________________________\n")
            features_list.insert(INSERT, fasta.nucleotides_counter(fasta.full_seq(fa_list))+'\n')
            features_list.insert(INSERT, fasta.fasta_resolve(fa_list))
        if file_format == 'GB' or file_format == 'GBK' or file_format == 'GBFF':
            ORGANISM = re.compile(r'ORGANISM  .*')
            organism_list = re.findall(ORGANISM, file_text)
            organism_set = set(organism_list)
            for h in organism_set:
                features_list.insert(INSERT, h+'\n')
            assemblyX = genbank.gb_obtain(file_text)
            fa_list = assemblyX.pop(0)
            assembly_annotations = assemblyX.pop()
            full_assembly = str('')
            for g in assembly_annotations:
                for i in fa_list:
                    if i == g:
                        features_list.insert(INSERT, '\n'+i)
                        nucl = fasta.nucleotides_counter(fa_list[i])
                        features_list.insert(INSERT, '\n'+nucl)
                        annot = genbank.features_counter(assembly_annotations[g])
                        features_list.insert(INSERT, '\n'+annot)
            features_list.insert(INSERT, "_____________________________________________\n")
            features_list.insert(INSERT, fasta.nucleotides_counter(fasta.full_seq(fa_list))+'\n')
            features_list.insert(INSERT, fasta.fasta_resolve(fa_list))
            features_list.insert(INSERT, genbank.features_counter(file_text))
            features_list.insert(INSERT, genbank.all_genes_names(file_text))
        features_list.insert(INSERT, '\n******************************************* \n')
    def clear():
        features_list.delete(1.0, END)
    def show_file():
        check = islogfile.get()
        if check == 1:
            pass
        else:
            StartWindow.clear()
        try:
            file_text = StartWindow.file_open()
            features_list.insert(INSERT, file_text)
        except:
            messagebox.showinfo('ОШИБКА', 'Невозможно открыть файл')
    def save_file_simple():
        check = isblocked.get()
        if check == True:
            messagebox.showinfo('ОШИБКА', 'Перезапись запрещена. Разблокируйте перезапись')
        if check == False:
            try:
                path_file = path_entry.get()
                file_info = features_list.get(1.0, END)
                with open(path_file, 'w+') as output:
                    print(file_info, file=output)
            except:
                messagebox.showinfo('ОШИБКА', 'Не удалось сохранить файл')
    def save_file_as():
        try:
            save_as = fd.asksaveasfilename(defaultextension=".txt")
            file_text = features_list.get(1.0, END)
            f = open(save_as, "w")
            f.write(file_text)
            f.close()
        except:
            messagebox.showinfo('ОШИБКА', 'Не удалось сохранить файл')
    def exit_file():
        message = "Вы уверены, что хотите закрыть программу?"
        if messagebox.askyesno(message=message):
            window.destroy()
    def spravka():
        message = 'Программа UltraЧТЕЦ 3.0 - простой текстовый редактор, созданный для быстрого просмотра файлов нуклеотидных и аминокислотных последовательностей, а также для частичного анализа геномных сборок. Права не защищены, никто не защищён, мир опасен. Берегите себя, удачного использования программы!'
        messagebox.showinfo("Справка", message)
    def help():
        path = os.getcwd()
        try:
            path_file = os.path.join('UCH_instruction.txt')
            with open(path_file, 'r') as file:
                text = file.read()
            features_list.delete(1.0, END)
            features_list.insert(INSERT, text)
        except:
            messagebox.showinfo('Упс', 'Кажется, мы потеряли инструкцию :[ Можете попытаться найти её и вернуть на место в ту же папку, где лежит программа. Её зовут UCH_instruction.txt')
class TextFindField():
    def text_find_field():
        textfindfield = tk.Toplevel()
        textfindfield.title("Введите поисковый запрос")
        global key_plus
        key_entry_plus = Entry(textfindfield, textvariable=key_plus, width=80, font = (12))
        finder_button = Button(textfindfield, text="Поиск", command=TextFindField.find_text, width=10, font = (12))
        key_entry_plus.pack(side=LEFT)
        finder_button.pack(side=LEFT)
    def reset_list():
        global search_list
        global keyword_2
        j = key_plus.get()
        if keyword_2 != j:
            keyword_2 = key_plus.get()
            search_list.clear()
            features_list.tag_remove(SEL, 1.0, "end-1c")
    def find_text():
        TextFindField.reset_list()
        global search_list
        global keyword_2
        features_list.focus_set()
        if keyword_2:
            if search_list == []:
                idx = "1.0"
            else:
                idx = search_list[-1]
            idx = features_list.search(keyword_2, idx, nocase=1, stopindex=END)
            lastidx = '%s+%dc' % (idx, len(keyword_2))

            try:
                features_list.tag_remove(SEL, 1.0,lastidx)
            except:
                pass
            try:
                features_list.tag_add(SEL, idx, lastidx)
                counter_list = []
                counter_list = str(idx).split('.')
                features_list.mark_set("insert", "%d.%d" % (float(int(counter_list[0])), float(int(counter_list[1]))))
                features_list.see(float(int(counter_list[0])))
                search_list.append(lastidx)
            except:
                messagebox.showinfo("Поиск завершён","Больше совпадений нет")
                search_list.clear()
                features_list.tag_remove(SEL, 1.0,"end-1c")
class FinderField():
    def finder_field():
        global seq_format
        global isfasta
        global key_entry
        global isscanning
        finderfield = tk.Toplevel()
        finderfield.title("Введите поисковый запрос")
        key_entry = Entry(finderfield, width=80, font = (12))
        key_entry.pack(side=TOP)
        # опции выбора
        fna = Radiobutton(finderfield, text='получить fna (для файлов анотаций)', value='fna', variable=seq_format, font = (12)).pack(anchor=SW)
        faa = Radiobutton(finderfield, text='получить faa (для файлов аннотаций)', value='faa', variable=seq_format, font = (12)).pack(anchor=SW)
        report = Radiobutton(finderfield, text='получить отчёт (для файлов аннотаций)', value='rep',variable=seq_format, font = (12)).pack(anchor=SW)
        saving_bloked = Radiobutton(finderfield, text="искать по заголовкам fasta (для fasta-файлов)", value='fas', variable=seq_format, font = (12)).pack(anchor=SW)
        scanning_button = Checkbutton(finderfield, text="исследовать текст в диалоговом окне", onvalue=1, offvalue=0, variable=isscanning, font = (12)).pack(anchor=SW)
        # небольшое описание
        global output_label
        finder_button = Button(finderfield, text="Поиск", command=FinderField.start_finder, width=10, font = (12))
        finder_button.pack(side=BOTTOM)
        output_label = Label(finderfield,text="", justify=CENTER, font=(14))
        output_label.pack(side=BOTTOM)
    def start_finder():
        scan_check = isscanning.get()
        seq_check = seq_format.get()
        keyword = key_entry.get()
        file_text = str("")
        error = False
        if scan_check == '0':
            try:
                file_text = StartWindow.file_open()
            except:
                messagebox.showinfo('ОШИБКА', 'Файл не найден')
                error = True
        if error == False:
            if scan_check == '1':
                file_text = features_list.get(1.0, END)
            assembly_format = StartWindow.good_choice(file_text)
            if assembly_format == 'FASTA':
                if seq_check == 'fas':
                    FinderField.fasta_finder(file_text, keyword)
                    output_label.configure(text = "Обнаружено совпадений: {keyword_number}".format(keyword_number=keyword_number))
                if seq_check != 'fas':
                    messagebox.showinfo('ОШИБКА', 'Выбранный файл не соответствует запрошенному формату')
            if assembly_format == 'bad':
                messagebox.showinfo('ОШИБКА', 'Этот текст не содержит аннотаций')
            if assembly_format != 'FASTA' and assembly_format != 'bad':
                if seq_check == 'fas':
                    messagebox.showinfo('ОШИБКА', "Выбранный файл не соответствует запрошенному формату")
                if seq_check != 'fas':
                    FinderField.gb_finder(file_text, keyword)
                    if seq_check == 'rep':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}".format(keyword_number=keyword_number))
                    if seq_check == 'fna':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}, записано {writing_number}".format(keyword_number=keyword_number, writing_number=writing_number))
                    if seq_check == 'faa':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}\nпоследовательности аминокислот есть у {writing_number} из них".format(keyword_number=keyword_number, writing_number=writing_number))
    def gb_finder(text, keyword):
        StartWindow.clear()
        def inner_insert(text):
            features_list.insert(INSERT, text+'\n')
        # переменные
        global keyword_number
        global writing_number
        global features_list
        seq_check = seq_format.get()
        keyword_number = 0
        writing_number = 0
        # компильные регулярки
        title = re.compile(r' {5}\S+  +')
        gene_name = re.compile(r'/product=\".*?\"', re.DOTALL)
        gene_edges = re.compile(r'\d+\W+\d+|complement\(\d+\W+\d+?\)')
        gene_id = re.compile(r'/protein_id=\".*?\"', re.DOTALL)
        gene_translation = re.compile(r'/translation=\"[\w\s]*\"')
        annot_split = re.compile(r'/.+=\\')
        # аннотации и последовательности
        assemblyX = genbank.gb_obtain(text)
        assembly_annotations = assemblyX.pop()
        assembly_sequence = assemblyX.pop(0)

        for i in assembly_annotations:
            for h in assembly_sequence:
                if h==i:
                    scaffold_list = re.split(title, assembly_annotations[i])
                    scaffold_sequence = assembly_sequence[h]
                    for j in scaffold_list:
                        if '/product="' in j:
                            g_name = re.search(gene_name, j)
                            gg_name = re.sub(r' +', ' ', (g_name.group()).replace('\n', ''))
                            if keyword in gg_name:
                                keyword_number += 1
                                g_edges = re.search(gene_edges, j)
                                gg_edges = g_edges.group()
                                if seq_check == 'rep':
                                    inner_insert(gg_name)
                                    inner_insert('/source="'+i+'"')
                                    inner_insert('/position="'+gg_edges+'"')
                                try:
                                    edges = re.findall('\d+', gg_edges)
                                    u = 0
                                    na_sequence = str('')
                                    for h in edges:
                                        u += 1
                                        if u == 1:
                                            st = int(h)
                                        if u == 2:
                                            fin = int(h)
                                            delta_l = abs(fin-st)
                                            na_sequence = scaffold_sequence[st:st+delta_l]
                                            break
                                    if 'complement' in gg_edges:
                                        na_sequence = fasta.complement_sequence(na_sequence)
                                    if seq_check == 'rep':
                                        inner_insert('/sequence="'+na_sequence+'"')
                                except:
                                    na_sequence = None
                                try:
                                    g_id = re.search(gene_id, j)
                                    gg_id = re.sub(r' +', ' ', g_id.group().replace('\n', ''))
                                    if seq_check == 'rep':
                                        inner_insert(gg_id)
                                except:
                                    gg_id = None
                                try:
                                    g_translation = re.search(gene_translation, j)
                                    gg_translation = re.sub(r'\n| *', '', g_translation.group())
                                    if seq_check == 'rep':
                                        inner_insert(gg_translation)
                                except:
                                    gg_translation = None
                                if seq_check == 'rep':
                                    features_list.insert(INSERT, '\n')
                                if seq_check == 'fna':
                                    if na_sequence != '' and na_sequence != ' ' and na_sequence != None:
                                        writing_number += 1
                                        gg_name1 = gg_name.replace('/product=\"', '>').replace('"', '')
                                        if gg_name1 == '>' or gg_name1 == '> ':
                                            gg_name1 = '>Unknown_'+str(writing_number)
                                        na_sequence1 = fasta.beautifull_sequence(na_sequence)
                                        inner_insert(gg_name1+'\n'+na_sequence1)
                                if seq_check == 'faa':
                                    if keyword_number == 1:
                                        StartWindow.clear()
                                    if gg_translation != None:
                                        writing_number += 1
                                        gg_name1 = gg_name.replace('/product=\"', '>').replace('"', '')
                                        gg_translation1 = gg_translation.replace('/translation=\"', '').replace('"', '')
                                        gg_translation = fasta.beautifull_sequence(gg_translation1)
                                        if gg_name1 == '>' or gg_name1 == '> ':
                                            gg_name1 = '>Unknown_'+str(writing_number)
                                        inner_insert(gg_name1+'\n'+gg_translation)
    def fasta_finder(text, keyword):
        fasta_assembly = fasta.fasta_obtain(text)
        global keyword_number
        keyword_number = 0
        for i in fasta_assembly:
            if keyword in i:
                keyword_number += 1
                if keyword_number == 1:
                    StartWindow.clear()
                features_list.insert(INSERT, i+'\n'+tw.fill(fasta_assembly[i], width=60)+'\n')
if __name__ == '__main__':
    # само окно
    window = Tk()
    StartWindow(window).grid()
    window.resizable(True, True)
    window.grid_columnconfigure(0, weight=1)
    window.title("UltraЧТЕЦ 3.0")

    # переменные
    seq_format = StringVar()
    isscanning = StringVar()
    islogfile = BooleanVar()
    isblocked = BooleanVar()
    isfasta = BooleanVar()
    pathX = StringVar()
    key_plus = StringVar()
    search_list = list()
    keyword_2 = str("")

    # задание значений по умолчанию
    isfasta.set(False)
    isscanning.set(0)
    isblocked.set(False)
    seq_format.set('rep')

    SIZE = 12
    def increase():
        global SIZE
        SIZE += 1
        features_list['font']=("Liberation Mono", SIZE)
        numbers['font']=("Liberation Mono", SIZE)
        path_input['font']=("Liberation Mono", SIZE)
        path_entry['font']=("Liberation Mono", SIZE)
    def decrease():
        global SIZE
        if SIZE > 1:
            SIZE -= 1
        features_list['font']=("Liberation Mono", SIZE)
        numbers['font']=("Liberation Mono", SIZE)
        path_input['font']=("Liberation Mono", SIZE)
        path_entry['font']=("Liberation Mono", SIZE)

    # создание меню
    mainmenu = Menu(window)
    window.config(menu=mainmenu)

    filemenu = Menu(mainmenu, tearoff=0)
    filemenu.add_command(label="Открыть...", command = StartWindow.path_input, font = (12))
    filemenu.add_command(label="Сохранить", command = StartWindow.save_file_simple, font = (12))
    filemenu.add_command(label="Сохранить как...", command = StartWindow.save_file_as, font = (12))
    filemenu.add_command(label="Выход", command = StartWindow.exit_file, font = (12))

    editmenu = Menu(mainmenu, tearoff=0)
    editmenu.add_command(label="Просмотр содержимого", command = StartWindow.show_file, font = (12))
    editmenu.add_command(label="Увеличить текст", command = increase, font = (12))
    editmenu.add_command(label="Уменьшить текст", command = decrease, font = (12))
    editmenu.add_command(label="Очистить поле", command = StartWindow.clear, font = (12))

    specmenu = Menu(mainmenu, tearoff=0)
    specmenu.add_command(label="Анализ сборки/генома", command = StartWindow.input_assembly, font = (12))
    specmenu.add_command(label="Простой поиск", command = TextFindField.text_find_field, font = (12))
    specmenu.add_command(label="Поиск по аннотациям", command = FinderField.finder_field, font = (12))
    specmenu.add_checkbutton(label="Запретить перезапись файла", onvalue=1, offvalue=0, variable=isblocked, font=(12))
    specmenu.add_checkbutton(label="Режим свитка", onvalue=1, offvalue=0, variable=islogfile, font=(12))

    helpmenu = Menu(mainmenu, tearoff=0)
    helpmenu.add_command(label="Помощь", command = StartWindow.help, font = (12))
    helpmenu.add_command(label="О программе", command = StartWindow.spravka, font = (12))

    mainmenu.add_cascade(label="Файл", menu=filemenu, font = (12))
    mainmenu.add_cascade(label="Правка", menu=editmenu, font = (12))
    mainmenu.add_cascade(label="Специальные возможности", menu=specmenu, font = (12))
    mainmenu.add_cascade(label="Справка", menu=helpmenu, font = (12))

    # фреймы
    path_input = LabelFrame(window, text="Путь для поиска файла", labelanchor=NE, font=(12), padx=5, pady=5)
    path_input.grid_columnconfigure(0, weight=1)
    output = Frame(window)
    output.grid_columnconfigure(0, weight=1)

    # элементы окна
    path_entry = Entry(path_input, textvariable=pathX, font = ("Liberation Mono", SIZE))

    # Магия, если честно, сам не до конца понимаю
    numbers = Text(output, width=6, bg='lightgray', state=DISABLED, relief=FLAT, font = ("Liberation Mono", SIZE))
    scroll_y = ttk.Scrollbar(output, orient=tk.VERTICAL)
    scroll_x = ttk.Scrollbar(output, orient=tk.HORIZONTAL)

    def on_yscrollcommand(*args):
        scroll_y.set(*args)  # Синхронизация скролбара с текстовым полем
        numbers.yview_moveto(args[0])  # Синхронизация поля с номерами с текстовым полем
    def on_xscrollcommand(*args):
        scroll_x.set(*args)  # Синхронизация скролбара с текстовым полем

    features_list = Text(output, xscrollcommand=on_xscrollcommand, yscrollcommand=on_yscrollcommand, wrap=NONE, font = ("Liberation Mono", SIZE))

    def scroll_command_y(*args):
        # Движение скролбара управляет отображением текста в обоих текстовых полях
        features_list.yview(*args)
        numbers.yview(*args)

    def scroll_command_x(*args):
        # Движение скролбара управляет отображением текста в обоих текстовых полях
        features_list.xview(*args)

    scroll_y.config(command=scroll_command_y)
    scroll_x.config(command=scroll_command_x)

    def insert_numbers():
        count_of_lines = features_list.get(1.0, END).count('\n') + 1
        #numbers.tag_configure('RIGHT', justify=RIGHT)
        numbers.config(state=NORMAL)
        numbers.delete(1.0, END)
        numbers.insert(1.0, '\n'.join(map(str, range(1, count_of_lines))))
        numbers.config(state=DISABLED)

    insert_numbers()

    def on_edit(event):
        # Срабатывает при изменениях в текстовом поле
        insert_numbers()
        features_list.edit_modified(0)  # Сбрасываем флаг изменения текстового поля

    features_list.bind('<<Modified>>', on_edit)

    # Нужно чтобы текстовое поле автоматически меняло размер при изменении размера окна
    window.grid_rowconfigure(1, weight=1)
    output.grid_columnconfigure(1, weight=50)
    output.grid_rowconfigure(0, weight=60)

    # Не трогать!
    path_input.grid(column=0, row=0, sticky = 'NSWE')
    output.grid(column=0, row=1, sticky = 'NSWE')

    path_entry.grid(column=0, row=0, sticky = 'NSWE')
    path_entry.focus()

    numbers.grid(column=0, row=0, sticky='NSWE')
    features_list.grid(column=1, row=0, sticky='NSWE')
    scroll_y.grid(column=2, row=0, sticky='NSWE')
    scroll_x.grid(column=1, row=2, sticky='NSWE')

    window.mainloop()
