from grid import Grid
import random


def jumpdown(grid, x, n):
    """
    >>> grid = Grid.build([['s', 's'], ['r', None], [None, 'r']])
    >>> jumpdown(grid, 0, 2)
    [[None, 's'], ['r', None], ['s', 'r']]
    >>> grid = Grid.build([['s', 's'], ['r', None], [None, 'r']])
    >>> jumpdown(grid, 1, 2)
    [['s', None], ['r', None], [None, 'r']]
    >>> grid = Grid.build([['s', None], ['s', None], [None, None]])
    >>> jumpdown(grid, 0, 1)
    [[None, None], ['s', None], ['s', None]]
    """
    for y in reversed(range(grid.height)):
        if grid.get(x, y) == 's':
            grid.set(x, y, None)
            dest_y = y + n
            if grid.in_bounds(x, dest_y):
                if grid.get(x, dest_y) is None:
                    grid.set(x, dest_y, 's')
    return grid


def is_scared(grid, x, y):
    """
    >>> grid = Grid.build([[None, None], ['y', 'b'], [None, 'p']])
    >>> is_scared(grid, 0, 1)
    True
    >>> grid = Grid.build([[None, None], ['y', 'b'], [None, 'p']])
    >>> is_scared(grid, 1, 2)
    True
    >>> grid = Grid.build([[None, None], ['y', None], [None, 'p']])
    >>> is_scared(grid, 0, 1)
    False
    """
    if grid.get(x, y) == 'y' or grid.get(x, y) == 'p':
        if grid.in_bounds(x, y - 1) and grid.get(x, y - 1) == 'b':
            return True
        if grid.in_bounds(x + 1, y) and grid.get(x + 1, y) == 'b':
            return True
    return False


def run_away(grid):
    """
    >>> grid = Grid.build([[None, None, None], [None, 'y', 'b'], [None, None, 'p']])
    >>> run_away(grid)
    [[None, None, None], ['y', None, 'b'], [None, 'y', None]]
    """
    for x in range(grid.width):
        for y in range(grid.height):
            if is_scared(grid, x, y):
                grid.set(x, y, 'y')
                if grid.in_bounds(x - 1, y) and grid.get(x - 1, y) is None:
                    grid.set(x, y, None)
                    grid.set(x - 1, y, 'y')
    return grid

def add_mine_counts(grid):
    """
    This function adds mine counts to the minesweeper grid. Spaces with
    mines ('m') are left alone. If the space is empty, counts the number of
    mines in the 8 surrounding spaces and places that number in the space.

    >>> grid=Grid.build([[None, None, None], ['m', None, None], [None, None, 'm']])
    >>> add_mine_counts(grid)
    [[1, 1, 0], ['m', 2, 1], [1, 2, 'm']]
    """
    for x in range(grid.width):
        for y in range(grid.height):
            space_val = grid.get(x, y)
            if space_val != 'm':
                num_mines = get_num_mines(grid, x, y)
                grid.set(x, y, num_mines)
    return grid


# initialize the grid with random values
def initialize_grid(grid, choices):
    random.seed(1)
    for i in range(grid.width):
        for j in range(grid.height):
            grid.set(i, j, random.choice(choices))
    return grid


def get_num_mines(grid, x, y):
    """
    This function has been implemented for you. It takes in a grid filled
    with mines and a target square. It returns the number of mines in the
    8 spaces surrounding the target square.
    """
    num_mines = 0
    for change_x in range(-1, 2):  # need to evaluate changes of -1, 0, and 1
        for change_y in range(-1, 2):
            # don't assess the target space
            if change_x != 0 or change_y != 0:
                curr_x = x + change_x
                curr_y = y + change_y
                # don't count out of bounds or non-mine spaces
                if grid.in_bounds(curr_x, curr_y) and grid.get(curr_x, curr_y) == 'm':
                    num_mines += 1
    return num_mines
