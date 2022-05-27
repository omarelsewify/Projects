#!/usr/bin/env python3

"""
Stanford CS106A Python DrawCanvas
Provides an on screen canvas with basic drawing functions.
Nick Parlante

Minor cleanup Jan-16-2020
Updated Windows-10 fix Oct-9-2019
Draft version - Oct-3-2019
See test_main() below for sample client code.
"""

import tkinter


"""
DrawCanvas - create one of these and draw to it like:

    canvas = DrawCanvas(200, 100, title='Draw Test')
    canvas.draw_rect(0, 0, 50, 60, color='red')
    ...
    
Call DrawCanvas.mainloop() at the bottom of main()
"""


class DrawCanvas(object):
    def __init__(self, width, height, fast_draw=True, title=None):
        """
        Creates a new on-screen drawing canvas.
        """
        self.canvas = make_canvas(width, height, title)
        self.auto_update = not fast_draw
        # Could add: background color

    def draw_line(self, x1, y1, x2, y2, color='black'):
        """
        Draws a black line between points x1,y1 and x2,y2
        Optional color='red' parameter can specify a color.
        """
        self.canvas.create_line(x1, y1, x2, y2, fill=color)
        if self.auto_update:
            self.canvas.update()

    def draw_rect(self, x, y, width, height, color='black'):
        """
        Draws a 1 pixel rectangle frame with its upper left at x,y
        and covering width, height pixels.
        Takes optional color='red' parameter.
        """
        if type(color) == tuple:
            color = DrawCanvas.color_name(color)
        self.canvas.create_rectangle(x, y, x + width - 1, y + height - 1, outline=color)
        if self.auto_update:
            self.canvas.update()

    def fill_rect(self, x, y, width, height, color='black'):
        """
        Draws a solid black rectangle with its upper left at x,y
        and covering width, height pixels.
        Takes optional color='red' parameter.
        """
        if type(color) == tuple:
            color = DrawCanvas.color_name(color)
        self.canvas.create_rectangle(x, y, x + width - 1, y + height - 1, outline=color, fill=color)
        # tricky: spec both fill and outline to get simple filled rect
        if self.auto_update:
            self.canvas.update()

    def draw_oval(self, x, y, width, height, color='black'):
        """
        Draws a 1 pixel oval frame with its upper left bounding rect at x,y
        and covering width, height pixels.
        Takes optional color='red' parameter.
        """
        if type(color) == tuple:
            color = DrawCanvas.color_name(color)
        self.canvas.create_oval(x, y, x + width - 1, y + height - 1, outline=color)
        if self.auto_update:
            self.canvas.update()

    def fill_oval(self, x, y, width, height, color='black'):
        """
        Draws a solid black oval with its upper left bounding rect at x,y
        and covering width, height pixels.
        Takes optional color='red' parameter.
        """
        if type(color) == tuple:
            color = DrawCanvas.color_name(color)
        self.canvas.create_oval(x, y, x + width - 1, y + height - 1, outline=color, fill=color)
        if self.auto_update:
            self.canvas.update()

    def draw_string(self, x, y, text, color='black'):
        """
        Draws a black text string with its upper left at x,y
        Takes optional color='red' parameter.
        """
        if type(color) == tuple:
            color = DrawCanvas.color_name(color)
        self.canvas.create_text(x, y, text=text, anchor=tkinter.NW, fill=color, font=('Courier', 24))
        if self.auto_update:
            self.canvas.update()

    def erase(self):
        """
        Erases all the canvas contents
        """
        self.canvas.delete('all')
        if self.auto_update:
            self.canvas.update()

    def set_fast_draw(self, fast_draw):
        """
        Sets the fast_draw boolean
        """
        self.auto_update = not fast_draw

    def update(self):
        """
        Update the onscreen pixels to reflect all the drawing.
        Normally drawing code does not need to do this.
        """
        self.canvas.update()

    @staticmethod
    def color_name(rgb):
        """
        Internal Utility. Given rgb tuple,
        form the '#ff2233' form used by TK. Generates readable
        exceptions when values not in 0..255
        We will coerce to int silently
        >>> DrawCanvas.color_name((255, 1, 0))
        '#ff0100'
        """
        if len(rgb) != 3:
            raise Exception('RGB error, expected 3-items but got:' + str(rgb))

        # Could coerce to int, then do checks

        if rgb[0] < 0 or rgb[1] < 0 or rgb[2] < 0:
            raise Exception('RGB error, negative value:' + str(rgb))

        if rgb[0] > 255 or rgb[1] > 255 or rgb[2] > 255:
            raise Exception('RGB error, value over 255:' + str(rgb))

        return '#{:02x}{:02x}{:02x}'.format(int(rgb[0]), int(rgb[1]), int(rgb[2]))

    @staticmethod
    def mainloop():
        """
        Calls the tkinter.mainloop(), typically last line of main().
        This version checks that there is a window on screen first,
        doing nothing if there is no window.
        """
        if tkinter._default_root:
            tkinter._default_root.update()
            tkinter.mainloop()

    """A few of the TK color constant names"""
    COLORS = ['red', 'orange', 'yellow', 'green', 'blue', 'lightblue', 'purple',
              'darkred', 'darkgreen', 'darkblue', 'pink', 'black', 'gray']


def make_canvas(width, height, title=None):
    """
    Creates and returns a TK drawing canvas
    of the given int size.
    This code can be used within a TK application
    to make a window suitable for TK drawing.
    Optional title parameter setting the title of the window.
    """
    top = tkinter.Tk()

    # The magic numbers below are a total hack, as TK incorrectly clips
    # the drawing canvas. Tested on Mac OS X and Windows 10
    # with the drawcanvas main() and the quilt solution.
    # Had to change these numbers around until it works on both
    # Windows 10 and Mac OS X.
    top.minsize(width=width + 10, height=height + 10)
    if title:
        top.title(title)

    canvas = tkinter.Canvas(top, width=width + 2, height=height + 2)
    canvas.pack()
    canvas.xview_scroll(8, 'units')  # hack so (0, 0) works correctly
    canvas.yview_scroll(8, 'units')  # otherwise it's clipped off

    return canvas


def test_canvas(width, height):
    """
    Creates and draws on DrawCanvas as a test.
    """
    canvas = DrawCanvas(width, height, title='Draw Test')

    canvas.draw_rect(0, 0, width, height, color='red')
    canvas.fill_oval(0, 0, width, height, color=(100, 100, 200))  # rgb tuple form

    n = 30
    for i in range(n):
        x = (i / (n - 1)) * (width - 1)
        canvas.draw_line(0, 0, x, height - 1, color='blue')
    canvas.draw_string(10, 10, 'Behold my pixels ye mighty and despair!')


def main():
    test_canvas(800, 400)
    DrawCanvas.mainloop()


if __name__ == '__main__':
    main()
