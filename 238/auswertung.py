#!python3

import labtools

import experiment1
import experiment2

def main():
    tasks = {
        'b': experiment1.b,
        'c': experiment2.c,
    }

    if labtools.task_list.run_task_list(tasks):
        labtools.defaults.print_notes()
    else:
        labtools.defaults.write_notes('results/notes.txt')

    labtools.defaults.merge_all()


if __name__ == '__main__':
    main()
