#!python3

import labtools as tools
from sys import argv

import exp_1
import exp_2


def main():

	task_list = {
		'c': exp_1.c,
		'h': exp_2.h,
	}

	preview = tools.task_list.run_task_list(task_list)

	if not preview:
		tools.notes.write_notes('results/notes.txt')
	else:
		tools.notes.print_notes()


if __name__ == '__main__':
	main()