
def in_range(n, low, high):
    """
    An example doctest is below. Add more to test your function.
    >>> in_range(3, 0, 100)
    True
    """
    if n <= high and n >= low:
        return True
    return False


def is_even(n):
    """
    Returns True if n is even and False otherwise.
    >>> is_even(69)
    False
    """
    if n % 2 == 0:
        return True
    return False


def only_one_even(num1, num2):
    """
    Returns True if exactly one input is even and False otherwise.
    >>> only_one_even(3,3)
    False
    """
    if is_even(num1) and not is_even(num2):
        return True
    if is_even(num2) and not is_even(num1):
        return True
    return False


def is_prime(n):
    """
    Returns True if n is prime and False otherwise.
    """
    pass


def main():
    print("\nChecking solutions to problems...\n")

    print("Calling in_range(4, 1, 5): returns " + str(in_range(4, 1, 5)))
    print("Calling is_even(9): returns " + str(is_even(9)))
    print("Calling only_one_even(2, 16): returns " + str(only_one_even(2, 16)))
    print("Calling is_prime(11): returns " + str(is_prime(11)))



if __name__ == "__main__":
    main()