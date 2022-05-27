from grid import Grid
from drawcanvas import DrawCanvas

ARM_LENGTH = 100 # given constant for draw snowman

CELL_SIZE = 100 # given constant for visualize grid

def fancy_lines(canvas, x, y, width, height, n):
    """
    This function takes in a left, top, width, height and number of lines n.
    The function should leave a 10-pixel high margin across the top and bottom
    of the figure without any drawing. Given int n which is 2 or more, draw n lines,
    as follows: The first line should begin at the left edge, 10 pixels from the top,
    and extend to the right edge, 10 pixels from the bottom. The last line should start
    at the left edge, 10 pixels from the bottom, and extend to the right edge 10 pixels
    from the top, with the other lines distributed proportionately.
    """
    for i in range(n):
        y_add = i * ((height - 1 - 20)//(n-1))
        canvas.draw_line(x, y + 10 + y_add, x + width - 1, y + height -1 - 10 - y_add)



def diagonal_lines(canvas, x, y, width, height, n):
    """
    This function takes in a left, top, width, height
    and number of lines n. Given int n which is 2 or more,
    draw n lines, each starting on the left edge and ending
    at the lower right corner of the figure. The first line
    should start at the upper left corner, the last line should
    start at the lower left corner, with the other lines distributed evenly.
    """
    for i in range(n):
        y_add = i * ((height - 1)//(n - 1))
        canvas.draw_line(x, y + y_add, x + width - 1, y + height - 1)

def draw_bullseye(canvas, width, height, radii, colors):
    """
    This function draws a bullseye centered on the given canvas. The radius and color
    of each ring are given by the passed in lists radii and colors, respectively.
    If there are more rings than colors, then colors are repeated.
    """
    mid_x = (width - 1) // 2
    mid_y = (height - 1) // 2
    for i in range(len(radii)):
        curr_color = colors[i % len(colors)]
        canvas.fill_oval(mid_x - radii[i], mid_y - radii[i], radii[i]*2, radii[i]*2, curr_color)


def draw_snowman(canvas, width, height, radii):
    """
    This function draws a snowman centered along the bottom of the given canvas.
    The snowman consists of circles, one on top of another, with the smallest
    on top and the largest on the bottom; the radius for each is given by the
    passed in radii list (smallest radius is first, largest is last). Two arms
    are drawn protruding from the middle circle of the snowman.
    """
    mid_x = (width - 1)//2
    new_bottom = height
    for i in reversed(range(len(radii))):
        x = mid_x - radii[i]
        y = new_bottom - (2*radii[i])
        canvas.draw_oval(x, y, 2*radii[i], 2*radii[i])
        if i % 3 == 1:
            canvas.draw_line(x, y + radii[i], x - ARM_LENGTH, y + radii[i])
            canvas.draw_line(x + 2*radii[i], y + radii[i], x + 2*radii[i] + ARM_LENGTH, y + radii[i])
        new_bottom = new_bottom - 2*radii[i]


def visualize_grid(grid):
    """
    Working with data structures like Grids can get a little tricky when
    you can't visualize what you're working with!
​
    Let's code the visualize_grid(grid) function, that takes a Grid as input
    and draws the grid to the canvas.
​
    Each box in the grid should have height and width CELL_SIZE. We also want
    to print the contents of each grid cell as a string on the canvas!
​
    You can assume that the grid's contents are all strings.
    >>> grid = Grid.build([['1', '2', '3'], ['4', '5', '6'], ['7', '8', '9']])
    >>> visualize_grid(grid)
    """
    canvas = DrawCanvas(grid.width*CELL_SIZE, grid.height*CELL_SIZE)
    for row in range(grid.width):
        for col in range(grid.height):
            canvas.draw_rect(CELL_SIZE*col, CELL_SIZE*row, CELL_SIZE, CELL_SIZE)
            canvas.draw_string(CELL_SIZE*col, CELL_SIZE*row, grid.get(col, row))

def main():
    width = 1000
    height = 600
    canvas = DrawCanvas(width, height)
    grid = Grid.build([['1', '2', '3'], ['4', '5', '6'], ['7', '8', '9']])
    visualize_grid(grid)
    # draw_snowman(canvas, 500, 300, [25, 50, 100])
    # draw_bullseye(canvas, 500, 300, [200, 125, 100, 50, 20, 5], ['red', 'green', 'blue', 'yellow', 'orange'])
    canvas.mainloop()

if __name__ == "__main__":
    main()