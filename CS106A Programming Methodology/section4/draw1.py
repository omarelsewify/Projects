import tkinter
from drawcanvas import DrawCanvas

# given constants
SQUARE_SIDE_LENGTH_RATIO = 0.27
NUM_STRIPES = 9

def draw_stanford_flag(canvas, left_x, top_y, width, height):
    """
    This function draws the outline of the new Stanford flag at the given location
    with the given width and height. The parameters are:
        - canvas: the canvas on which to draw the flag
        - left_x: the x coordinate for the upper left corner of the flag
        - top_y: the y coordinate for the upper left corner of the flag
        - width: the total width of the flag
        - height: the total height of the flag
    """
    pass

def draw_flags(canvas, width, height):
    """
    This function draws three identical flags on the given canvas, one in the bottom
    left corner, one in the top center, and one in the bottom right corner. The parameters
    are:
        - canvas: the canvas on which to draw the flags
        - width: the total width of the canvas
        - height: the total height of the canvas
    """
    pass

def main():

    # These test the draw_flags function
    width = 900
    height = 475
    canvas = DrawCanvas(width, height)
    draw_flags(canvas, width, height)
    tkinter.mainloop()


if __name__ == "__main__":
    main()