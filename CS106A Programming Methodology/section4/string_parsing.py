
def exclaim_words(s):
    """
    This function takes in a string and returns a list of the exclamatory words in that string.
    Consider an exclamatory word the "word" substring of one or more alphabetic chars
    which are immediately to the left of the '!' in s.

    >>> exclaim_words('x123hello!32happy!day')
    ['hello', 'happy']
    """
    search = 0
    result = []
    while search < len(s):
        # Find '[', at search index
        ex = s.find('!', search)
        # No '[' -> exit loop
        if ex == -1:
            break
        result.append(s[ex-5:ex])
        # Update search for next iteration
        search = ex + 1
    return result




def parse_unique_categories(s):
    """
    This function takes in a line from a Pylib, return a list of all unique categories listed in
    that line (omitting the delimiting brackets).

    >>> parse_unique_categories("The first place I am going is [place]. I cannot wait to [verb] when I get to [place]!")
    ['place', 'verb']
    """
    pass

def parse_words(s):
    """
    This function takes in a string representing an email in the format firstname.lastname@domainname.com, and returns
    a list of the email broken into three parts: the first name, the last name, and the domain name.

    >>> parse_words('theodor.geisel@seuss.com')
    ['theodor', 'geisell', 'seuss']

    >>> parse_words('maya.angelou@gmail.com')
    ['maya', 'angelou', 'gmail']
    """
    pass

def parse_out_hashtags(s):
    """
    This function takes in a string representing a single tweet and returns a list of all hashtags in the tweet.

    >>> parse_out_hashtags('I am going to #wear #my #mask')
    ['wear', 'my', 'mask']

    >>> parse_out_hashtags("I miss Stanford #backtocampus #ResX?")
    ['backtocampus', 'ResX']
    """
    pass

def parse_phone_number(s):
    """
    This function uses a while loop accompanied with str.find() to parse a phone number into a list of strings
    where the first string is the area code of the given phone number and the second string is the rest of
    the phone number with the dash characters removed. This function will work for all phone numbers,
    not just 10 digit ones.

    >>> parse_phone_number('212-225-9876')
    ['212', '2259876']
    """
    pass

def main():

    # Add code here to test the functions
    pass

if __name__ == "__main__":
    main()