#!python3
from labtools.defaults import *


from exp import c


def main():
    tasks = {
        'c': c
    }

    preview = run_task_list(tasks)
    
    if preview:
        print_notes()
    else:
        write_notes('results/notes.txt')
    

if __name__ == '__main__':
    main()
