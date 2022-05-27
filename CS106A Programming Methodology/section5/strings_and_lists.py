def censor_chat(message, pattern):
    """
    >>> msg = "My name is Karel, and I live in Wilbur Hall!"
    >>> ptn = "           xxxxx                xxxxxxxxxxx "
    >>> censor_chat(msg, ptn)
    "My name is xxxxx, and I live in xxxxxxxxxxx!"
    """
    result = ""
    for i in range(len(pattern)):
        if pattern[i] == "x":
            result += "x"
        if pattern[i] != "x":
            result += message[i]
    return results


def encrypt_double(src1, slug1, src2, slug2, char):
    """
    >>> # Make len-4 lists for testing
    >>> src1  = ['a', 'b', 'c', 'd']
    >>> slug1 = ['s', 't', 'u', 'v']
    >>> src2  = ['e', 'f', 'g', 'h']
    >>> slug2 = ['w', 'x', 'y', 'z']
    >>> encrypt_double(src1, slug1, src2, slug2, 'A')
    's'
    >>> encrypt_double(src1, slug1, src2, slug2, 'e')
    'z'
    >>> encrypt_double(src1, slug1, src2, slug2, 'f')
    'y'
    >>> encrypt_double(src1, slug1, src2, slug2, '$')
    '$'
    """
    low = char.lower()
    if low in src1:
        pos = src1.index(low)
        return slug1[pos]
    if low in src2:
        pos = src2.index(low)
        return slug2[len(slug2) - pos - 1]
    else:
        return char
