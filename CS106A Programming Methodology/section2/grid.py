#!/usr/bin/env python3

"""
Stanford CS106A Grid class
Nick Parlante
Provides simple 2d storage
"""


def grid_demo():
    """
    Demonstrate use of the Grid class.
    """
    # Create width 4 by height 2 grid, filled with None initially.
    grid = Grid(4, 2)

    # loop over contents in usual way,
    # setting every location to 6
    for y in range(grid.height):
        for x in range(grid.width):
            grid.set(x, y, 6)


    # access 0,0
    val = grid.get(0, 0)

    # verify that 3,1 is in bounds
    if grid.in_bounds(3, 1):
        # set 3,1
        grid.set(3, 1, 11)

    print(grid)
    # print uses nested-list format
    # showing row-0, then row-1
    # [[6, 6, 6, 6], [6, 6, 6, 11]]

    # Grid.build() also uses nested-list format, allowing
    # you to construct a grid on the fly.
    grid2 = Grid.build([[6, 6, 6, 6], [6, 6, 6, 11]])


class Grid:
    """
    Grid with y,x internal storage
    Has .width .height size properties
    """
    def __init__(self, width, height):
        """
        Create grid width by height.
        Initially all locations hold None.
        """
        # Pretty agro use of comprehensions!
        self.array = [[None for x in range(width)] for y in range(height)]
        self.width = width
        self.height = height

    @staticmethod
    def build(lst):
        """
        Utility.
        Construct Grid using a nested-lst literal
        e.g. this makes a 3 by 2 grid:
        Grid.build([[1, 2, 3], [4, 5 6]])
        >>> Grid.build([[1, 2, 3], [4, 5, 6]])
        [[1, 2, 3], [4, 5, 6]]
        """
        check_list_malformed(lst)
        height = len(lst)
        width = len(lst[0])
        grid = Grid(width, height)
        grid.array = lst  # slight waste, but keeps ctor params simple
        return grid

    def get(self, x, y):
        """Gets the stored value at x,y. None by default."""
        return self.array[y][x]

    def set(self, x, y, val):
        """Sets a new value into the grid at x,y."""
        self.array[y][x] = val

    def in_bounds(self, x, y):
        """Returns True if the x,y is in bounds of the grid. False otherwise."""
        return x >= 0 and x < self.width and y >= 0 and y < self.height

    def __str__(self):
        return repr(self.array)

    # In particular Doctest seems to use this, so crucial to make
    # Grid work in Doctests.
    def __repr__(self):
        return repr(self.array)


def check_list_malformed(lst):
    """
    Given a list that represents a 2-d nesting, checks that it has the
    right type and the sublists are all the same len.
    Raises exception for malformations.
    """
    if not lst or type(lst) != list:
        raise Exception('Expecting list but got:' + str(lst))

    if len(lst) >= 2:
        size = len(lst[0])
        for sub in lst:
            if len(sub) != size:
                raise Exception("Sub-lists are not all the same length:" + str(lst))




def main():
    grid_demo()


if __name__ == '__main__':
    main()


