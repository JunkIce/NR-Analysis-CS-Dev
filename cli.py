import master as nrat
import json
import os
import importlib
import sys

def clear(args=None):
    os.system('cls' if os.name == 'nt' else 'clear')

def terminal_width():
    try:
        terminal_size = os.get_terminal_size()
        return terminal_size.columns
    except OSError:
        return 80  # fallback width

def printu(string:str):
    print(string.encode('utf-8').decode('unicode_escape'))


config_filename='cli_stuff.json'
config_file=open(config_filename,'r')

config=json.load(config_file)
clear()
for line in config['init_seq']:
    printu(line)


def help(args):
    clear()
    # Basic Help Menu. I'll honestly throw this in a TXT file later because this sucks lol. Good for a basic reference ig
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

def pdb_load(codes):
    for code in codes:
        nrat.init(code)

def local_load(paths):
    for path in paths:
        nrat.init(path)

def load(args):

    load_valid={
        '-pdb':pdb_load,
        '-local':local_load
    }

    if not args or args[0] not in load_valid:
        print('''For loading and parsing coordinate files. use 'h' for more info\n''')
        return
    else:
        load_valid[args[0]](args[1:])

def list(args):
    if args:
        if args[0] in nrat.loaded_files:
            print('Chains in '+(target:=nrat.loaded_files[args[0]]).name + ':')
            i=1
            for chain in (chains:=target.chains):
                ligand_str=''
                if chains[chain].ligands:
                    ligand_list={}
                    for ligand in chains[chain].ligands:
                        if ligand not in ligand_list:
                            ligand_list[ligand]=1
                            continue
                        else:
                            ligand_list[ligand]+=1
                    
                    ligand_str=' with '
                    for ligand in ligand_list:
                        ligand_str+=ligand.full_name() + (f'(x{ligand_list[ligand]})' if ligand_list[ligand]>1 else '') + ', '

                print(f'{i}) '+chains[chain].type() + ligand_str)
                i+=1
            print()
    elif nrat.loaded_files:
        print('Loaded Structures:\n>','\n> '.join(nrat.loaded_files.keys()),'\n')
    else:
        print('No Loaded Structures. Use "load" to initialize a structure.')


def test(args):
    print(args)

def reload(args):
    pass

valid_commands={
    # Help Commands
    'h':help,
    'help':help,
    'clear':clear,

    'load':load,
    'list':list,

    # Dev Tools
    'test':test,
    #'reload':reload,
    #'r':reload
}


def main():
    run=True
    dev_mode=False

    while run:
        input_text=input('> ')
        if not input_text:
            continue
        command,*args=input_text.split()

        if command in ('exit','quit','q'):
            break
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