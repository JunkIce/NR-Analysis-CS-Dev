import threading
import time
def spin(stop:threading.Event):
    i=0
    strs=['   ','.  ','.. ','...']
    while not stop.is_set():
        print(f'\033[2;33mLoading{strs[i%len(strs)]}\033[0m',end='\r')
        time.sleep(0.15)
        i+=1
    print('          ',end='\r')
    
stop=threading.Event()
spinner=threading.Thread(target=spin, args=(stop,))
spinner.start()

import master as nrat
import json
import os
import subprocess
from typing import Iterable


def clear(args=None) -> None:
    '''
    **Quick macro that clears the terminal**
    '''

    os.system('cls' if os.name == 'nt' else 'clear')
    # 'clear' works for unix and powershell
    # just have to check for windows incase some psycho decides to keep using cmd

    return


def terminal_width() -> int:
    '''
    **Gets the terminal witdth (columns)**

    If terminal size cannot be determined, falls back to 80 (std term width)
    '''
    width:int = 80 # default size
    try:
        terminal_size = os.get_terminal_size()
        width = terminal_size.columns
    except OSError:
        pass # If py.os freaks out, just keep default
    
    return width


def printu(string:str) -> None:
    '''
    Prints while interpreting unicode escape characters
    
    they just make json all screwy, I'll probabally get rid of this at some point
    '''
    print(string.encode('utf-8').decode('unicode_escape'))
    return


config_filename='cli_stuff.json' # Probabally going to change up the config format at some point, but for some reason i thought json was a good idea
config_file=open(config_filename,'r')
config=json.load(config_file)


def title(args=None) -> None:
    '''
    Clears the terminal and prints the out the title text defined in config
    '''
    clear()
    for line in config['init_seq']:
        printu(line)
    return

def not_loaded() -> None:
    '''
    Default for if invalid arguments are given for a number of functions
    ''' 
    print('File(s) not in loaded files. Load using "load"')


'''
*** COMMANDS ***
'''

def help(args:Iterable) -> None:
    '''
    Displays Help Menu
    '''
    clear()
    # Basic Help Menu. I'll throw this in a TXT file later because this sucks lol. Good for a basic reference ig
    if not args:
        print(
f'''\033[1;32mNRAT Help Menu\033[0m

Command:                    │ Function:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿{'━'*(terminal_width()-29)}
Commands                    │
\033[1mProgram Control\033[0m             │
> h|help                    │ This Menu
> exit|q|quit               │ Exits the program
────────────────────────────┼{'─'*(terminal_width()-29)}
\033[1mStructure Management:\033[0m       │
> load:                     │
    Source Control:         │
    -pdb [pdb id]           │ Searches for and loads a structure from the PDB
    -loc/-local [path]      │ Loads a locally stored coordinate file
    -mass [dir path]        │ Loads all coordinate files from a given folder
    -list [path]            │ Loads coordinate files from newline-separated PDB 
                            │   codes. Add arguments after codes.
    Programatic Control:    │
    -n|-name [string]       │ Specifies structure name(s). Defaults to PDB 
                            │   code/File name.
> l|list                    │ Lists all loaded structure files
> sel|select [struct name]  │ Select a structure
> group [struct names]      │ Define a group of structures
────────────────────────────┼{'─'*(terminal_width()-29)}
\033[1mAnalysis:\033[0m                   │
> testfor                   │
    any                     │
    
'''
        )
        input('press enter to exit...')
        clear()
    return


def pdb_load(codes:Iterable) -> None:
    '''
    For initializing PDB entries
    '''
    for code in codes:
        try:
            nrat.init(code)
        except ValueError as e:
            print(f'Could not retrieve "{code}" from PDB.')
        if not nrat.selected:
            select([codes[0]])
    return

def local_load(paths:Iterable) -> None:
    '''
    For initializing local struct files
    '''
    for path in paths:
        try:
            nrat.init(path)
        except ValueError as e:
            print(f'Could not initialize "{path}"')
        if not nrat.selected:
            select([paths[0]])
    return

def load(args:Iterable)-> None:
    '''
    Base init command

    Passes args to pdb or local load
    '''
    load_valid={
        '-pdb':pdb_load,
        '-local':local_load
    }

    if not args or args[0] not in load_valid: # Fallback
        print('''For loading and parsing coordinate files. use 'h' for more info\n''')
    else:
        load_valid[args[0]](args[1:]) # pdb/local_load have same formatting. Honestly should probabally just be the same thing but oh well
    return


def list(args:Iterable) -> None:
    '''
    Gives list of loaded structures, or substructures if a loaded structure is passed
    '''
    if args: # If a structure is passed
        if args[0] in nrat.loaded_files:
            print('Chains in '+(target:=nrat.loaded_files[args[0]]).name + ':')
            i=1
            for chain in (chains:=target.chains):
                ligand_str=''
                if (ligands:=chains[chain].ligands):
                    ligand_names=[ligand.full_name() for ligand in ligands]

                    condensed_ligand_names=[]
                    for name in ligand_names:
                        if (count:=ligand_names.count(name)) > 1: # want to avoid redundant ligand names (ie. 'Glycerol (x3)' instead of 'Glycerol, Glycerol, Glycerol')
                            condensed_ligand_names.append(name+f' (x{count})')
                            while name in ligand_names:
                                ligand_names.pop(ligand_names.index(name))
                        else:
                            condensed_ligand_names.append(name)
                    
                    ligand_str=' with '+', '.join(condensed_ligand_names) # ie ' with GW4064, Glycerol (x3)'


                print(f'{i}) '+chains[chain].type() + ligand_str) #ie. '1) FXR with GW4064, Glycerol (x3)'
                i+=1
            print()

        # Displays all structures in a given group    
        elif args[0] in nrat.groups:
            print(f'files in group "{args[0]}":')
            i=1
            for file in nrat.groups[args[0]]:
                print(f'{i}) '+file.name)
                i+=1

    elif nrat.loaded_files: # Displays all loaded structures
        print('Loaded Structures:')
        i=1
        for file in nrat.loaded_files:
            print(('*'if nrat.selected.name==nrat.loaded_files[file].name else ' ')+f'{i}) "{file}"')
            i+=1
    else: # Fallback if nothing loaded
        print('No Loaded Structures. Use "load" to initialize a structure.')
    
    return


def group(args:Iterable) -> None:
    '''
    For making a group of structures to perform bulk actions
    '''
    if len(args)<2:return # can't initialize a group with no structures

    name,*files=args
    print('Defining group: '+name)
    structs=[]
    for file in files:
        if file in nrat.loaded_files:
            structs.append(nrat.loaded_files[file])
    nrat.group(name,structs)
    return


def select(args:Iterable) -> None:
    '''
    Selects an item (First is default)
    '''
    if len(args)<1:return # can't select nothing

    target=args[0]
    if target in nrat.loaded_files:
        nrat.selected=nrat.loaded_files[target]
        print(f'Selected: "{nrat.selected.name}"')
    else:
        not_loaded()
    return


def display(args) -> None:
    '''
    **WIP** Loads a structure file in ChimeraX (if on PATH)
    '''
    if len(args)<1:return
    target=args[0]
    if target in nrat.loaded_files:
        try:
            subprocess.Popen(['chimerax',nrat.loaded_files[target].filePath])
            print(f'Loaded "{target}" in ChimeraX. (This might take a second)')
        except:
            print('Could not initialize ChimeraX. Is it on PATH?')
    else:
        not_loaded()
    return


def analyze(args):
    '''
    For getting calculated info about a structure (distances, etc.)
    '''
    if len(args)<2:return
    targets,*funcs=args[0]
    valid_funcs={
        'pi'
    }

    if not all(func in valid_funcs for func in funcs):
        print('Invalid Args')
        return
    
    if targets not in nrat.groups:
        targets=[targets]

    if any(target not in nrat.loaded_files for target in targets):
        not_loaded()
        return
    
    for target in targets:
        file=nrat.loaded_files[target]
        for chain in file.chains:
            if (NR:=file.chains[chain]).family() != 'NR':
                continue

            NR.align()
            #NR.getInteractions()


    

# Debug to make sure the terminal is still alive
def test(args):
    print(args)

#def reload(args):
#    pass


valid_commands={
    # Help Commands
    'h':help,
    'help':help,
    'clear':clear,
    'c':clear,
    'version':title,
    'v':title,

    'load':load,
    'ld':load,
    'l':list,
    'list':list,
    'group':group,
    'g':group,
    'select':select,
    'sel':select,
    'disp':display,

    'analyze':analyze,
    'a':analyze,
}
'''List of all valid public commands. Some are aliases'''

dev_commands={
    # Dev Tools
    'test':test,
    #'reload':reload,
    #'r':reload
}
'''List of Dev Commands'''

stop.set()
spinner.join()

def main():
    '''
    Primary Run Loop
    '''
    run:bool = True
    dev_mode:bool = False
    
    # Run title for the begining of the program
    title()
    
    while run:
        input_text=input('> ')
        if not input_text:
            continue
        command,*args=input_text.split()

        if command in ('exit','quit','q'):
            break
        elif command=='DEV_MODE':
            dev_mode=True
        elif dev_mode and command in dev_commands:
            dev_commands[command](args)
        elif command in valid_commands:
            valid_commands[command](args)

        else:
            print('\033[1;32mCommand Not Recognized. Enter \'h\' for a list of valid commands.\033[0m')
    clear()
    print('\033[1;32mExiting...\033[0m')
    print('thanks for using my tool!                  '.rjust(terminal_width()))
    print('i hope it helps you figure something out :)'.rjust(terminal_width()))
    print('-cs                                        '.rjust(terminal_width()))

if __name__=='__main__':
    main()