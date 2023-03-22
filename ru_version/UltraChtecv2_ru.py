#! /usr/bin/python3

import re
import os
import sys
import datetime
import platform
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog as fd
from tkinter import font

import fastaUCHv2 as fa
import genbankUCHv2 as gb

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
            fa_list = fa.fasta_obtain(file_text)
            full_assembly = str('')
            for i in fa_list:
                features_list.insert(INSERT, '\n'+i)
                nucl = fa.nucleotides_counter(fa_list[i])
                features_list.insert(INSERT, '\n'+nucl+'\n')
            features_list.insert(INSERT, "_____________________________________________\n")
            features_list.insert(INSERT, fa.nucleotides_counter(fa.full_fasta_sequence(fa_list))+'\n')
            features_list.insert(INSERT, fa.fasta_resolve(fa_list))
        if file_format == 'GB' or file_format == 'GBK' or file_format == 'GBFF':
            ORGANISM = re.compile(r'ORGANISM  .*')
            organism_list = re.findall(ORGANISM, file_text)
            organism_set = set(organism_list)
            for h in organism_set:
                features_list.insert(INSERT, h+'\n')
            assemblyX = gb.gb_obtain(file_text)
            fa_list = assemblyX.pop(0)
            assembly_annotations = assemblyX.pop()
            full_assembly = str('')
            for g in assembly_annotations:
                for i in fa_list:
                    if i == g:
                        features_list.insert(INSERT, '\n'+i)
                        nucl = fa.nucleotides_counter(fa_list[i])
                        features_list.insert(INSERT, '\n'+nucl)
                        annot = gb.features_counter(assembly_annotations[g])
                        features_list.insert(INSERT, '\n'+annot)
            features_list.insert(INSERT, "_____________________________________________\n")
            features_list.insert(INSERT, fa.nucleotides_counter(fa.full_fasta_sequence(fa_list))+'\n')
            features_list.insert(INSERT, fa.fasta_resolve(fa_list))
            features_list.insert(INSERT, gb.features_counter(file_text))
            features_list.insert(INSERT, gb.coding_koeff(file_text))
            features_list.insert(INSERT, gb.all_genes_names(file_text))
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
        message = 'Программа UltraЧТЕЦ v.2 - простой текстовый редактор, созданный для быстрого просмотра файлов нуклеотидных и аминокислотных последовательностей, а также для частичного анализа геномных сборок'
        messagebox.showinfo("Справка", message)
    def help():
        path = os.getcwd()
        try:
            path_file = os.path.join('Manual_ru.txt')
            with open(path_file, 'r') as file:
                text = file.read()
            features_list.delete(1.0, END)
            features_list.insert(INSERT, text)
        except:
            messagebox.showinfo('Упс', 'Кажется, мы потеряли инструкцию :[ Можете попытаться найти её и вернуть на место в ту же папку, где лежит программа')
        

# окно простого поиска в тексте
class TextFindField():
    def text_find_field():
        global key_plus
        # имя поля
        textfindfield = tk.Toplevel()
        textfindfield.title("Введите поисковый запрос")
        # задание параметров окна
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

# окно усложнённого поиска по аннотациям
class FinderField():
    # задание окна
    def finder_field():
        global seq_format
        global isfasta
        global key_entry
        global isscanning
        global isscroll
        global SIZE
        global SHRIFT
        # имя окна - finderfield
        finderfield = tk.Toplevel()
        finderfield.title("Введите поисковый запрос")
        key_entry = Entry(finderfield, width=80, font = (SHRIFT, SIZE))
        key_entry.pack(side=TOP)
        # опции выбора
        fna = Radiobutton(finderfield, text='получить fna (для файлов анотаций)', value='fna', variable=seq_format, font = (SIZE)).pack(anchor=SW)
        faa = Radiobutton(finderfield, text='получить faa (для файлов аннотаций)', value='faa', variable=seq_format, font = (SIZE)).pack(anchor=SW)
        report = Radiobutton(finderfield, text='получить отчёт (для файлов аннотаций)', value='rep',variable=seq_format, font = (SIZE)).pack(anchor=SW)
        saving_bloked = Radiobutton(finderfield, text="искать по заголовкам fasta (для fasta-файлов)", value='fas', variable=seq_format, font = (SIZE)).pack(anchor=SW)
        scanning_button = Checkbutton(finderfield, text="исследовать текст в диалоговом окне", onvalue=1, offvalue=0, variable=isscanning, font = (SIZE)).pack(anchor=SW)
        scroll_mode = Checkbutton(finderfield, text="режим свитка при поиске", onvalue=1, offvalue=0, variable=isscroll, font = (SIZE)).pack(anchor=SW)

        # небольшое описание
        global output_label
        finder_button = Button(finderfield, text="Поиск", command=FinderField.start_finder, width=10, font = (SIZE))
        finder_button.pack(side=BOTTOM)
        output_label = Label(finderfield,text="", justify=CENTER, font=(SIZE))
        output_label.pack(side=BOTTOM)

    def start_finder():
        # опция сканирования поля
        scan_check = isscanning.get()
        # определение формата
        seq_check = seq_format.get()
        # очищает поле, если выключен режим свитка
        scroll_check = isscroll.get()
        if scroll_check == False:
            StartWindow.clear()
        # определение ключевого слова
        keyword = key_entry.get()
        key = re.compile(r'[A-Za-z0-9]*'+keyword+r'[A-Za-z0-9]*') ###ВНИМАНИЕ!!! Это выражение используется в последующем для сложного поиска
        islogfile
        file_text = str("")
        error = False
        # открытие файла
        if scan_check == '0':
            try:
                file_text = StartWindow.file_open()
            except:
                messagebox.showinfo('ОШИБКА', 'Файл не найден')
                error = True
        # если открывается успешно, то
        if error == False:
            if scan_check == '1':
                file_text = features_list.get(1.0, END)
            assembly_format = StartWindow.good_choice(file_text)
            if assembly_format == 'FASTA':
                if seq_check == 'fas':
                    # если последовательность fasta, и выбрана правильная опция "поиск по заголовкам fasta", то должен искать
                    FinderField.fasta_finder(file_text, key)
                    output_label.configure(text = "Обнаружено совпадений: {keyword_number}".format(keyword_number=keyword_number))
                if seq_check != 'fas':
                    messagebox.showinfo('ОШИБКА', 'Выбранный файл не соответствует запрошенному формату')
            if assembly_format == 'bad':
                messagebox.showinfo('ОШИБКА', 'Этот текст не содержит аннотаций')
            # если не ошибка и не fasta, то геном
            if assembly_format != 'FASTA' and assembly_format != 'bad':
                if seq_check == 'fas':
                    messagebox.showinfo('ОШИБКА', "Выбранный файл не соответствует запрошенному формату")
                if seq_check != 'fas':
                    FinderField.gb_finder(file_text, key, seq_check)
                    if seq_check == 'rep':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}".format(keyword_number=keyword_number))
                    if seq_check == 'fna':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}, записано {writing_number}".format(keyword_number=keyword_number, writing_number=writing_number))
                    if seq_check == 'faa':
                        output_label.configure(text = "Обнаружено совпадений: {keyword_number}\nпоследовательности аминокислот есть у {writing_number} из них".format(keyword_number=keyword_number, writing_number=writing_number))

    def gb_finder(text, key, seq_check):
        def inner_insert(text):
            features_list.insert(INSERT, text+'\n')
        # переменные
        global keyword_number
        global writing_number
        global features_list
        keyword_number = 0
        writing_number = 0
        # регулярные выражения
        title = re.compile(r' {5}\S+  +')
        gene_name = re.compile(r'/product=\".*?\"', re.DOTALL)
        gene_edges = re.compile(r'\d+\W+\d+|complement\(\d+\W+\d+?\)')
        gene_id = re.compile(r'/protein_id=\".*?\"', re.DOTALL)
        gene_translation = re.compile(r'/translation=\"[\w\s]*\"')
        annot_split = re.compile(r'/.+=\\')
        # поиск всех встречающихся выражений с keyword
        key_list = set(re.findall(key, text))
        # аннотации и последовательности
        assemblyX = gb.gb_obtain(text)
        assembly_annotations = assemblyX.pop()
        assembly_sequence = assemblyX.pop(0)
        # сопоставление списка последовательностей и списка аннотаций
        for i in assembly_annotations:
            for h in assembly_sequence:
                if h==i:
                    # получение заголовков
                    scaffold_list = re.split(title, assembly_annotations[i])
                    # получение последовательности
                    scaffold_sequence = assembly_sequence[h]
                    for j in scaffold_list:
                        if '/product="' in j:
                            g_name = re.search(gene_name, j)
                            key_presence = False
                            gg_name = re.sub(r' +', ' ', (g_name.group()).replace('\n', '')).replace(' ', '_')
                            # проверка наличия выражения с keyword в имени гена
                            for key_expression in key_list:
                                if key_expression in gg_name:
                                    key_presence = True
                                    break
                            # если есть, то
                            if key_presence == True:
                                keyword_number += 1
                                g_edges = re.search(gene_edges, j)
                                gg_edges = g_edges.group()
                                if seq_check == 'rep':
                                    inner_insert(gg_name)
                                    inner_insert('/source="'+i+'"')
                                    inner_insert('/position="'+gg_edges+'"')
                                try:
                                    # поиск границ генов в нуклеотидной последовательности
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
                                            na_sequence = scaffold_sequence[st-1:st+delta_l]
                                            break
                                    if 'complement' in gg_edges:
                                        na_sequence = fa.complement_sequence(na_sequence)
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
                                    if na_sequence:
                                        writing_number += 1
                                        gg_name1 = gg_name.replace('/product=\"', '>').replace('"', '')
                                        if gg_name1 == '>' or gg_name1 == '> ':
                                            gg_name1 = '>Unknown_'+str(writing_number)
                                        na_sequence1 = fa.beautifull_sequence(na_sequence)
                                        inner_insert(gg_name1+'\n'+na_sequence1)
                                if seq_check == 'faa':
                                    if gg_translation:
                                        writing_number += 1
                                        gg_name1 = gg_name.replace('/product=\"', '>').replace('"', '')
                                        gg_translation1 = gg_translation.replace('/translation=\"', '').replace('"', '')
                                        gg_translation = fa.beautifull_sequence(gg_translation1)
                                        if gg_name1 == '>' or gg_name1 == '> ':
                                            gg_name1 = '>Unknown_'+str(writing_number)
                                        inner_insert(gg_name1+'\n'+gg_translation)
    def fasta_finder(fasta_sorted_list, key):
        # специфический поиск по заголовкам fasta
        global keyword_number
        global islogfile
        keyword_number = 0
        for i in fasta_sorted_list:
            key_presence = re.search(key, i)
            if key_presence:
                keyword_number += 1
                features_list.insert(INSERT, i+'\n'+fa.beautifull_sequence(fasta_sorted_list[i])+'\n')

class Fonts_list():
    def choosing_fonts():
        global SHRIFT
        global SIZE
        def choose_fonts(event):
            global SHRIFT
            selected_indice = font_listbox.curselection()
            SHRIFT = font_listbox.get(selected_indice)
        def right_button():
            global SHRIFT
            features_list['font']=(SHRIFT, SIZE)
            numbers['font']=(SHRIFT, SIZE)
            path_input['font']=(SHRIFT, SIZE)
            path_entry['font']=(SHRIFT, SIZE)
            fontslist.destroy()
        def exit_fonts():
            fontslist.destroy()
        fontslist = tk.Toplevel()
        fontslist.title("Моноширинные шрифты")
        font_entry = Entry(fontslist, width=40, font=(SHRIFT, SIZE))
        font_listbox = Listbox(font_entry, width=30, listvariable=SHRIFT, font=(SHRIFT, SIZE))
        choose_button = Button(fontslist, text = "Выбрать", command = right_button)
        close_button = Button(fontslist, text = "Отмена", command = exit_fonts)
        list_fonts = list(font.families())
        for i in list_fonts:
            if "mono" in i or "Mono" in i:
                font_listbox.insert(END, i)
        font_entry.pack()
        font_listbox.pack()
        font_listbox.bind("<<ListboxSelect>>", choose_fonts)
        choose_button.pack()
        close_button.pack()
  
if __name__ == '__main__':
    # увеличить размер шрифта
    def increase():
        global SIZE
        global SHRIFT
        SIZE += 1
        features_list['font']=(SHRIFT, SIZE)
        numbers['font']=(SHRIFT, SIZE)
        path_input['font']=(SHRIFT, SIZE)
        path_entry['font']=(SHRIFT, SIZE)

    # уменьшить размер шрифта
    def decrease():
        global SIZE
        global SHRIFT
        if SIZE > 1:
            SIZE -= 1
        features_list['font']=(SHRIFT, SIZE)
        numbers['font']=(SHRIFT, SIZE)
        path_input['font']=(SHRIFT, SIZE)
        path_entry['font']=(SHRIFT, SIZE)

    # отмена действия в текстовом поле
    def undo():
        features_list.edit_undo()
        
    # возврат действия в текстовом поле
    def redo():
        features_list.edit_redo()

    # возможность выделения текста нажатием Ctrl-A
    def select_text(event):
        def selection(widget):
            widget.select_range(0, 'end')
            widget.icursor('end')
        window.after(10, selection, event.widget)

    # скопировать весь текст в поле вывода
    def all_text_copy():
        text = features_list.get('1.0', END)
        window.clipboard_clear()
        window.clipboard_append(text)
        window.update()
        
    # скопировать выделенный фрагмент
    def copy_text():
        selection = features_list.tag_ranges(SEL)
        if selection:
            window.clipboard_clear()
            window.clipboard_append(features_list.get(*selection))

    # вставить выделенный фрагмент
    def past_text():
        text = window.clipboard_get()
        if text:
            features_list.insert(INSERT, text)

    # вырезать текст
    def cut_text():
        copy_text()
        selection = features_list.tag_ranges(SEL)
        if selection:
            features_list.delete(*selection)

    # щелчок правой кнопкой мыши
    def right_click(event):
        global x, y
        x = event.x
        y = event.y
        clickmenu.post(event.x_root, event.y_root)

    # щелчок левой кнопкой мыши (убирает контекстное меню)
    def left_click(event):
        clickmenu.unpost()

    def on_yscrollcommand(*args):
        scroll_y.set(*args)  # Синхронизация скролбара с текстовым полем
        numbers.yview_moveto(args[0])  # Синхронизация поля с номерами с текстовым полем

    def on_xscrollcommand(*args):
        scroll_x.set(*args)  # Синхронизация скролбара с текстовым полем

    def scroll_command_y(*args):
        # Движение скролбара управляет отображением текста в обоих текстовых полях
        features_list.yview(*args)
        numbers.yview(*args)

    def scroll_command_x(*args):
        # Движение скролбара управляет отображением текста в обоих текстовых полях
        features_list.xview(*args)

    # вставка чисел в поле numbers
    def insert_numbers():
        count_of_lines = features_list.get(1.0, END).count('\n') + 1
        #numbers.tag_configure('RIGHT', justify=RIGHT)
        numbers.config(state=NORMAL)
        numbers.delete(1.0, END)
        numbers.insert(1.0, '\n'.join(map(str, range(1, count_of_lines))))
        numbers.config(state=DISABLED)

    # Срабатывает при изменениях в текстовом поле
    def edit_numbers(event):
        insert_numbers()
        features_list.edit_modified(0)  # Сбрасываем флаг изменения текстового поля
        
    # само окно
    window = Tk()
    StartWindow(window).grid()
    window.resizable(True, True)
    window.grid_columnconfigure(0, weight=1)
    window.title("UltraЧТЕЦ v.2")

    # переменные
    seq_format = StringVar()
    isscanning = StringVar()
    islogfile = BooleanVar()
    isblocked = BooleanVar()
    isscroll = BooleanVar()
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
    isscroll.set(False)
    SIZE = 14
    SHRIFT = "Liberation Mono"
    
    # создание меню
    mainmenu = Menu(window)
    window.config(menu=mainmenu)

    filemenu = Menu(mainmenu, tearoff=0)
    filemenu.add_command(label="Открыть...", command = StartWindow.path_input, font = (SHRIFT, SIZE))
    filemenu.add_command(label="Сохранить", command = StartWindow.save_file_simple, font = (SHRIFT, SIZE))
    filemenu.add_command(label="Сохранить как...", command = StartWindow.save_file_as, font = (SHRIFT, SIZE))
    filemenu.add_separator()
    filemenu.add_command(label="Выход", command = StartWindow.exit_file, font = (SHRIFT, SIZE))

    editmenu = Menu(mainmenu, tearoff=0)
    editmenu.add_command(label="Просмотр содержимого", command = StartWindow.show_file, font = (SHRIFT, SIZE))
    editmenu.add_command(label='Отменить действие', command = undo, font = (SHRIFT, SIZE))
    editmenu.add_command(label='Вернуть действие', command = redo, font = (SHRIFT, SIZE))
    editmenu.add_command(label='Скопировать текст', command = all_text_copy, font = (SHRIFT, SIZE))
    editmenu.add_command(label="Увеличить текст", command = increase, font = (SHRIFT, SIZE))
    editmenu.add_command(label="Уменьшить текст", command = decrease, font = (SHRIFT, SIZE))
    editmenu.add_command(label="Изменить шрифт", command=Fonts_list.choosing_fonts, font = (SHRIFT, SIZE))
    editmenu.add_command(label="Очистить поле", command = StartWindow.clear, font = (SHRIFT, SIZE))

    specmenu = Menu(mainmenu, tearoff=0)
    specmenu.add_command(label="Анализ генома", command = StartWindow.input_assembly, font = (SHRIFT, SIZE))
    specmenu.add_command(label="Простой поиск", command = TextFindField.text_find_field, font = (SHRIFT, SIZE))
    specmenu.add_command(label="Поиск по аннотациям", command = FinderField.finder_field, font = (SHRIFT, SIZE))

    opcmenu = Menu(mainmenu, tearoff=0)
    opcmenu.add_checkbutton(label="Запретить перезапись файла", onvalue=1, offvalue=0, variable=isblocked, font=(SHRIFT, SIZE))
    opcmenu.add_checkbutton(label="Режим свитка", onvalue=1, offvalue=0, variable=islogfile, font=(SHRIFT, SIZE))

    helpmenu = Menu(mainmenu, tearoff=0)
    helpmenu.add_command(label="Помощь", command = StartWindow.help, font = (SHRIFT, SIZE))
    helpmenu.add_command(label="О программе", command = StartWindow.spravka, font = (SHRIFT, SIZE))

    mainmenu.add_cascade(label="Файл", menu=filemenu, font = (SHRIFT, SIZE))
    mainmenu.add_cascade(label="Правка", menu=editmenu, font = (SHRIFT, SIZE))
    mainmenu.add_cascade(label="Возможности", menu=specmenu, font = (SHRIFT, SIZE))
    mainmenu.add_cascade(label="Режимы", menu=opcmenu, font = (SHRIFT, SIZE))
    mainmenu.add_cascade(label="Справка", menu=helpmenu, font = (SHRIFT, SIZE))

    clickmenu = Menu(mainmenu, tearoff=0)
    clickmenu.add_command(label="Копировать", command = copy_text, font = (SHRIFT, SIZE))
    clickmenu.add_command(label="Вставить", command = past_text, font = (SHRIFT, SIZE))
    clickmenu.add_command(label="Вырезать", command = cut_text, font = (SHRIFT, SIZE))

    # фреймы
    path_input = LabelFrame(window, text="Путь поиска файла", labelanchor=NE, font=(SHRIFT, SIZE), padx=5, pady=5)
    path_input.grid_columnconfigure(0, weight=1)
    output = Frame(window)
    output.grid_columnconfigure(0, weight=1)

    # элементы окна и действия с ними
    path_entry = Entry(path_input, textvariable=pathX, font = (SHRIFT, SIZE))
    path_entry.focus_set
    path_entry.bind('<Control-a>', select_text)
    path_entry.bind('<Control-z>', undo)
    path_entry.bind('<Control-Shift-z>', redo)

    # создание нумерации строк
    numbers = Text(output, width=6, bg='lightgray', state=DISABLED, relief=FLAT, font = (SHRIFT, SIZE))
    scroll_y = Scrollbar(output, orient=VERTICAL)
    scroll_x = Scrollbar(output, orient=HORIZONTAL)
    scroll_y.config(command=scroll_command_y)
    scroll_x.config(command=scroll_command_x)

    # задание тестового поля
    features_list = Text(output, xscrollcommand = on_xscrollcommand, yscrollcommand = on_yscrollcommand, wrap = NONE, undo = True, font = (SHRIFT, SIZE))
    features_list.focus_set

    insert_numbers()

    # Привязка действий в текстовом поле
    features_list.bind('<<Modified>>', edit_numbers)
    features_list.bind('<Control-z>', undo)
    features_list.bind('<Control-Shift-z>', redo)
    features_list.bind('<Button-3>', right_click)
    features_list.bind('<Button-1>', left_click)

    # Автоматическое изменение размеров окна под текст
    window.grid_rowconfigure(1, weight=1)
    output.grid_columnconfigure(1, weight=50)
    output.grid_rowconfigure(0, weight=60)

    # Расположение элементов в окне
    path_input.grid(column=0, row=0, sticky = 'NSWE')
    output.grid(column=0, row=1, sticky = 'NSWE')

    path_entry.grid(column=0, row=0, sticky = 'NSWE')
    path_entry.focus()

    numbers.grid(column=0, row=0, sticky='NSWE')
    features_list.grid(column=1, row=0, sticky='NSWE')
    scroll_y.grid(column=2, row=0, sticky='NSWE')
    scroll_x.grid(column=1, row=2, sticky='NSWE')

    window.mainloop()
