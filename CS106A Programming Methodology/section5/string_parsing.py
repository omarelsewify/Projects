def parse_liked(s, liked):
    """
    >>> parse_liked('xxabad..A..BB', ['a', 'b', 'c'])
    ['aba', 'BB']
    >>> parse_liked('^^AB.BBD.CD.', ['a', 'b', 'c'])
    ['AB', 'BB']
    """
    search = 0
    result = []
    while search < len(s):
        low = s.lower()
        word = ''
        while search < len(s) and low[search] in liked:
            word += s[search]
            search += 1
        if len(word) > 1:
            result.append(word)
        search += 1
    return result


def parse_backwards_hashtags(s):
    """
    >>> parse_backwards_hashtags('Check out how hip# I am!!#. Tell the 106A students# good luck!# on their quiz#')
    ['hip#', 'am!!#', 'students#', 'luck!#', 'quiz#']
    """
    tags = []
    search = 0
    while search < len(s):
        end = s.find('#', search)
        if end == -1:
            break
        start = end - 1
        while start >= 0 and (s[start].isalpha() or s[start] == '!'):
            start -= 1
        tag = s[start+1:end+1]
        if len(tag) >= 1:
            tags.append(tag)
        search = end + 1
    return tags
