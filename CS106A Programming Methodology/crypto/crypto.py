#!/usr/bin/env python3

"""
Stanford CS106A Crypto Project
"""

import sys

# provided ALPHABET constant - list of the regular alphabet
# in lowercase. Refer to this simply as ALPHABET in your code.
# This list should not be modified.
ALPHABET = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u',
            'v', 'w', 'x', 'y', 'z']


def compute_slug(key):
    """
    Given a key string, compute and return the len-26 slug list for it.
    >>> compute_slug('z')
    ['z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
    >>> compute_slug('Bananas!')
    ['b', 'a', 'n', 's', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'o', 'p', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    >>> compute_slug('Life, Liberty, and')
    ['l', 'i', 'f', 'e', 'b', 'r', 't', 'y', 'a', 'n', 'd', 'c', 'g', 'h', 'j', 'k', 'm', 'o', 'p', 'q', 's', 'u', 'v', 'w', 'x', 'z']
    >>> compute_slug('Zounds!')
    ['z', 'o', 'u', 'n', 'd', 's', 'a', 'b', 'c', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'p', 'q', 'r', 't', 'v', 'w', 'x', 'y']
    """
    code = []
    remainder = list.copy(ALPHABET)
    key = key.lower()
    # Pull out the letters used in the key and then add
    # the remaining alphabet letters to the end of the list
    for i in range(len(key)):
        if key[i].isalpha() and key[i] not in code:
            code.append(key[i])
        if key[i] in remainder:
            to_delete = remainder.index(key[i])
            del remainder[to_delete]
    return code + remainder


def encrypt_char(source, slug, char):
    """
    Given source and slug lists,
    if char is in source return
    its encrypted form. Otherwise
    return the char unchanged.
    >>> # Using 'z' slug for testing.
    >>> # Can set a var within a Doctest like this.
    >>> slug = ['z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
    >>> encrypt_char(ALPHABET, slug, 'A')
    'Z'
    >>> encrypt_char(ALPHABET, slug, 'n')
    'm'
    >>> encrypt_char(ALPHABET, slug, 'd')
    'c'
    >>> encrypt_char(ALPHABET, slug, '.')
    '.'
    >>> encrypt_char(ALPHABET, slug, '\\n')
    '\\n'
    """
    output = ''
    lower = char.lower()
    # Find the position of the letter in the source and then
    # match it to the corresponding position in the slug then
    # return the corresponding slug letter
    for i in range(len(lower)):
        if lower[i] in source:
            position = int(source.index(lower[i]))
            if position < len(slug):
                output = slug[position]
        if char[i].isupper():
            output = output.upper()
        if lower[i] not in source:
            output = char
    return output


def encrypt_str(source, slug, s):
    """
    Given source and slug lists and string s,
    return a version of s where every char
    has been encrypted by source/slug.
    >>> slug = ['z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
    >>> encrypt_str(ALPHABET, slug, 'And like a thunderbolt he falls.\\n')
    'Zmc khjd z sgtmcdqanks gd ezkkr.\\n'
    """
    # Apply encrypt_char to each character in the given string
    output_str = ''
    for i in range(len(s)):
        output_str += encrypt_char(source, slug, s[i])
    return output_str


def decrypt_str(source, slug, s):
    """
    Given source and slug lists, and encrypted string s,
    return the decrypted form of s.
    >>> slug = ['z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
    >>> decrypt_str(ALPHABET, slug, 'Zmc khjd z sgtmcdqanks gd ezkkr.\\n')
    'And like a thunderbolt he falls.\\n'
    """
    # Apply encrypt_char to each character in the given string
    # BUT NOW: slug is the source and source is the slug
    result = ''
    for i in range(len(s)):
        result = encrypt_str(slug, source, s)
    return result


def encrypt_file(filename, key):
    """
    Given filename and key, compute and
    print the encrypted form of its lines.
    """
    # Open a file and apply encrypt_str to each line of the file
    with open(filename) as f:
        for line in f:
            # look at line in loop
            print(encrypt_str(ALPHABET, compute_slug(key), line), end='')


def decrypt_file(filename, key):
    """
    Given filename and key, compute and
    print the decrypted form of its lines.
    """
    # Open a file and apply decrypt_str to each line of the file
    with open(filename) as f:
        for line in f:
            # look at line in loop
            print(decrypt_str(ALPHABET, compute_slug(key), line), end='')


def main():
    args = sys.argv[1:]
    # 2 command line argument patterns:
    # -encrypt key filename
    # -decrypt key filename
    # Call encrypt_file() or decrypt_file() based on the args.
    if len(args) == 3 and args[0] == '-encrypt':
        encrypt_file(args[2], args[1])
    if len(args) == 3 and args[0] == '-decrypt':
        decrypt_file(args[2], args[1])


# Python boilerplate.
if __name__ == '__main__':
    main()
