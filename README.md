# UltraChtec v2.0. Program manual
February 3, 2023.

This is a simple portable text editor for fast analysis of FASTA and GenBank genomic sequences formates. The program is written entirely in Python 3.10 and can be run from any operating system with Python3 installed. It consists of a main file (`UltraChtecv2.py`) with a size of about 679 lines and two modules:

`fastaUCHv2` is a small module for analyzing fasta sequences. The default module size is about 168 lines. 

`genbankUCHv2` is a module for parsing GenBank files. The module size is about 157 lines. Requires the previous module to work.

## USAGE
To work without installing Python3, there is a beta version for Windows in .exe format, so far only in Russian and with a limited set of functions. Along with the program, in the "samples" folder, there are two test genomic sequences as an example of incoming files. To run, you may need to install additional modules `re` and `tkinter`.

Modules can be installed through the console with the command: 
`pip3 install <module_name>` 

You can run the program, like any other Python code, using an interpreter program or from the console: 

`python3 UltraChtecv2.py` (from the folder where the file is located) 

or (for example):

`python3 C:\\Downloads\UltraChtecv2.py` (Windows) 

`python3 ~/Downloads/UltraChtecv2.py` (Linux)

## FUNCTIONS
### 1) Viewing and editing the content of genomic files. 
  To open a file, select the `Open...` item in the Menu and select it in the window that opens, or enter its address in the system manually in the `File path` field. The contents of the file will not appear in the text field until the `Edit – Viewing the content` option is selected. To save changes to a file, choose `File - Save` or `File - Save As...`. To prevent accidental saving of the changed text under the same name, select `Modes - Prevent file overwriting`. Managing actions in the text field is carried out using the menu: `Edit - Undo action` and `Edits - Redo action`. `Edit - Change font` changes the font of the workspace to any of the monospace fonts available in the tkinter.
  
  If the `Edit - Copy all text` option is selected, all text in the field is copied to the clipboard (you do not need to select it for full copying). You can reduce and enlarge the text using the corresponding options `Edit - Increase font` and `Edit - Decrease font`. All text can be removed using `Edit - Clear text field`. With any manipulations in the program, except for saving, the original file with the genome is not changed. By default, the program opens the file every time, determines its format, looks for matches in the headers. The program itself determines the file format, not focusing on the extension specified in the file name. This does not prevent the researcher from viewing and editing third-party files, but blocks the possibility of analysis.
  
### 2) Calculation of the main parameters of genomes

    • length (total, maximum for contigs, scaffolds)
    
    • number of scaffolds and contigs
    
    • A+G+T+C composition (%)
    
    • N50, N90
    
    • L50, L90
    
    • number and types of annotations (if the genome is annotated)
    
    • genome coding ratio (if the genome is annotated)
    
To carry out all of the above actions, you need to select `Features - Genome analysis`. The result will be visible in the text field, and include (in the order listed) the file path, file type, genome name, analysis time, analysis of contigs and/or scaffolds included in the genome, general analysis of the assembly.

### 3) Search for target genes
#### simple search in a file

To open the selected genome or report in the text window (`Edit – Viewing the content`) or paste it there after copying, select `Features - Simple search`, enter a word / part of a word in the text field and click `Search`.

#### search by annotations

This function allows you to find the name of a gene in the annotations by a keyword, and display its amino acid or nucleotide sequence (`Features — Annotation search`). The search is carried out both by the whole word and by its part. Looks for the necessary expression in the annotation headers (or FASTA file headers) and outputs (depending on the format and the selected target).

    input: .gb, .gbk or .gbff

    outputs: FASTA (`get fna` parameter), amino acids FASTA (`get faa`) or report* (`get report`)

Report includes the main elements of GenBank files: sequence name, contig (scaffold), individual gene number, nucleotide sequence, amino acid translation.

    input: FASTA

    output: FASTA (`search by FASTA headers`)
    
The option is relevant for FASTA files (.fna, .faa, etc.). For other types of files, an attempt to run this search will fail. 
However, if immediately after another search that resulted in a set of FASTA sequences, you select `Features — Annotation search – explore text in dialog box` at the same time as `Features — Annotation search – search FASTA headers`, then the search will be performed, allowing you to discard part of the unusable results.
     
The program cannot directly search for genes by their sequences, however, when you run `Features — Annotation search` with the `get report` option and with an empty search field, the program will return all genes in the form of reports with sequences without spaces and line breaks. By selecting `Features — Simple search` later, you can try to find a sequence that repeats the selected one.

 
### 4) Obtaining entire files with amino acid and nucleotide sequences from genomes
If you do not fill in the search field when searching by annotations (`Features — Annotation search`), then at startup the program will write out all the genes of the annotated genome in one of the required forms.

### 5) Reporting 
When the scroll mode is enabled, the program will make new records without erasing the old file. This function is available in the menu in the main field (`Modes - Scroll Mode`) and in the annotation search field (`Annotation search - scroll mode when searching`).

### 6) Naming unknown genes
If the genes do not have a name (typical for some annotating services), it is assigned automatically by the program. New names are formed according to the principle "Unknown_{number of unknown gene given by program}". This is convenient when generating FASTA files for further analysis in other programs.

## Contacts
For any questions feel free to email: 

lichoradkin43@gmail.com

Good luck to everyone!
