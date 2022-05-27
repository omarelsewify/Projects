from grid import Grid

def bubble_up(grid, y):
    """
    Implement this function as described in the handout. Add some more doctests.
    >>> grid = Grid.build([[None, 'r', 'b'], ['b', 'b', 'b']])
    >>> bubble_up(grid, 1)
    [['b', 'r', 'b'], [None, 'b', 'b']]
    """
    for i in range(grid.width):
        if grid.get(i, y) == 'b':
            if grid.in_bounds(i, y-1) and grid.get(i, y-1) == None:
                grid.set(i, y, None)
                grid.set(i, y-1, 'b')
    return grid

def main():
    print("\nChecking solutions to problems...\n")
    grid = Grid.build([[None, 'r', 'b'], ['b', 'b', 'b']])
    print("grid = Grid.build([[None, 'r', 'b'], ['b', 'b', 'b']])")
    bubble_up(grid, 1)
    print("Calling bubble_up(grid, 1) returns " + str(grid))

if __name__ == "__main__":
    main()