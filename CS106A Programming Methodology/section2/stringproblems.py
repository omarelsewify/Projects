


def remove_all_occurrences(s, to_remove):
    """
    Implements remove_all_occurrences as specified in the handout

    >>> remove_all_occurrences('This is a test', 't')
    'This is a es'
    >>> remove_all_occurrences('Summer is here!', 'e')
    'Summr is hr!'
    >>> remove_all_occurrences('----O----', '-')
    'O'

    :param s: the string to remove characters from
    :param to_remove: the string to be removed (which has a length of 1)
    :return: the filtered string
    """
    pass


def is_palindrome(s):
    """
    Implements is_palindrome as specified in the handout

    >>> is_palindrome('brahm')
    False
    >>> is_palindrome('maps spam')
    True
    >>> is_palindrome('Maps spam')
    False

    :param s: a string
    :return: whether s is a palindrome
    """
    pass



def exclaim(msg, end, n):
    """
    Implements exclaim as specified in the handout
    """
    pass



def main():
    print("\nChecking solutions to problems...\n")

    print("\n")
    print("Calling remove_all_occurrences('This is a test', 't'): returns " + str(remove_all_occurrences('This is a test', 't')))
    print("Calling remove_all_occurrences('Summer is here!', 'e'): returns " + str(remove_all_occurrences('Summer is here!', 'e')))
    print("Calling remove_all_occurrences('----O----', '-'): returns " + str(remove_all_occurrences('----O----', '-')))

    print("\n")
    print("Calling is_palindrome('brahm'): returns " + str(is_palindrome('brahm')))
    print("Calling is_palindrome('maps spam'): returns " + str(is_palindrome('maps spam')))
    print("Calling is_palindrome('Maps spam'): returns " + str(is_palindrome('Maps spam')))

    print("\n")
    print("Calling exclaim('Did you see the new Star Wars movie', '?!', 6)")
    exclaim('Did you see the new Star Wars movie', '?!', 6)
    print("Calling exclaim('106A is awesome', '!', 5)")
    exclaim('106A is awesome', '!', 5)



if __name__ == "__main__":
    main()