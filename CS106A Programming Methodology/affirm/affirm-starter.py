#!/usr/bin/env python3

"""
Stanford CS106A Command Line + Printing Example/Exercise
Nick Parlante

affirm.py          - working version
affirm-exercise.py - code in main() TBD

This program takes 3 command line arg forms:

1.
python3 affirm.py -affirm *name*

2.
python3 affirm.py -hello *name*

3.
python3 affirm.py -n *number* *name*
"""

import sys
import random

AFFIRMATIONS = [
    'Looking good',
    'All hail',
    'Horray for',
    'Today is the day for',
    'I have a good feeling about',
    'A big round of applause for',
    'Everything is coming up',
]


def print_hello(name):
    """
    Given name string, print 'Hello' with that name.
    """
    print('Hello', name)


def print_affirm(name):
    """
    Given name string, print a random affirmation for that name.
    """
    affirmation = random.choice(AFFIRMATIONS)
    print(affirmation, name)


def print_n_copies(n, name):
    """
    Given int n and name string, print n copies of that name.
    """
    for i in range(n):
        # Print each copy of the name with space instead of \n after it.
        print(name, end=' ')
    # Print a single \n after the whole thing.
    print()


def main():
    """
    This program takes 3 command line arg forms:

    1.
    python3 affirm.py -affirm *name*

    2.
    python3 affirm.py -hello *name*

    3.
    python3 affirm.py -n *number* *name*
    """
    # standard - load "args" list with cmd-line-args
    args = sys.argv[1:]

    # args is a list of the command line argument strings that follow
    # the program.py file.
    # So if the command is: python3 affirm.py aaa bbb ccc
    # The args list  will be:
    #   args == ['aaa', 'bbb', 'ccc']
    # The args are all strings, whatever was typed in the terminal.
    # print(args)

    # 1. Check for the -affirm arg pattern:
    #   python3 affirm.py -affirm Bart
    #   e.g. args[0] is '-affirm' and args[1] is 'Bart'
    pass

    # 2 - Check for -hello args pattern:
    #   python3 affirm.py -hello Lisa
    pass

    # 3 - Check for -n args pattern:
    #   python3 affirm.py -n 100 Maggie
    # Note: the command line arg is a string, convert to int
    pass


# Python boilerplate.
if __name__ == '__main__':
    main()
