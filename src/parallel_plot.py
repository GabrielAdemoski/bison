import tkinter as tk
import numpy as np
from PIL import Image
import os


def circle(canvas, x, y, rad, **kwargs):
    """
    Creates a tkinter oval object in the shape of a circle with its center at x,y and with
    a radius of rad.
    :param canvas: the canvas object on which the circle will be drawn
    :param x: x coordinate for the center of the circle
    :param y: y coordinate for the center of the circle
    :param rad: the radius of the resultant circle
    :param kwargs: tkinter create_oval arguments
    :return:
    """
    return canvas.create_oval(x - rad, y - rad, x + rad, y + rad, kwargs)


def from_to(matrix):
    """
    Creates tuples containing two nodes and the edge value between them. Only creates tuples for nodes with a nonzero
    edge. The first node is the row index and the second node is the column index.
    :param matrix: A 2D numpy array corresponding to an adjacency matrix
    :return: A list of tuples containing two nodes and the edge value between them
    """
    from_to = []
    r = 0
    if np.any(matrix):
        # matrix = matrix**3
        matrix /= np.max(abs(matrix[np.nonzero(matrix)]))
        for row in matrix:
            c = 0
            for val in row:
                if val != 0:
                    from_to.append((r, c, val))
                c += 1
            r += 1
    return from_to

def parallel_plot(matrix, taxa, disp=True, path="figs/", fname='parallel_image', title='X'):
    os.makedirs(path, exist_ok=True)

    width = 1000
    height = width / 4 * 3

    x_border = 100
    y_border = 200
    y_top = 50
    nsamples = matrix.shape[0]
    n_nodes = matrix.shape[1]
    offset = 5

    # start graphing lines from the top y_diff pixels down
    y_spacing = (height - 2 * y_border) / (nsamples)
    x_spacing = (width - 2 * x_border) / n_nodes
    htemp = y_border

    ys = [y_border - y_top + y_spacing * i for i in range(nsamples+1)]
    xs = [(x_border + (x_spacing / 2)) + x_spacing * i for i in range(n_nodes)]

    txt = ["2019", "2020", "2021"]
    OTUs = [i+1 for i in range(n_nodes)]
    # blue is negative and orange is positive
    line_color = ["blue", "orange"]

    # Create the main window
    window = tk.Tk()
    window.title("Example Years Canvas")
    # Create a canvas widget
    canvas = tk.Canvas(window, width=width, height=height, bg="white")
    canvas.pack()
    # Canvas title, centered with a downwards offset
    canvas.create_text(width / 2, y_border / 2, text=title, font=('Helvetica 25 bold'))

    for i in range(nsamples+1):
        # Draw a line representing the sample
        canvas.create_line((x_border, ys[i]), (width - x_border, ys[i]), fill='black', width=2)
        # Draw text with sample name next to sample line., the first argument of create_text()
        # is the x coord of the middle of the text
        canvas.create_text(x_border / 2, ys[i], text=txt[i], fill="black", font=('Helvetica 25 bold'))

        # Draw nodes
        for n in range(n_nodes):
            # Create vertical hash marks
            canvas.create_line((xs[n], ys[i]+offset), (xs[n], ys[i]-offset), fill='black', width=2)
            # Create a circle
            # circle(canvas, xs[n], ys[i], 5, fill='black', width=1)

            # Add tax names or OTU's at bottom of the canvas
            if i == nsamples:
                # Rotate text downwards with a downwards offset down
                canvas.create_text(xs[n], ys[i]+offset, text = taxa[n], anchor="w", angle=-90, font=('Helvetica 15 bold'))
                # No rotation, text left to right with a downwards offset
                # canvas.create_text(xs[n], ys[i] + offset, text=OTUs[n])
        # Draw edges between nodes
        for i in range(nsamples):
            for m in from_to(matrix[i]):
                weight = m[-1]
                line = canvas.create_line((xs[m[0]], ys[i]), (xs[m[1]], ys[i + 1]), fill=line_color[int(weight > 0)],
                                          width = abs(weight))
                canvas.lower(line)

        htemp += y_spacing

    canvas.update()
    path += fname
    canvas.postscript(file=path + '.eps', colormode='color')
    # use PIL to convert to PNG
    img = Image.open(path + '.eps')
    img.load(scale=5)
    img.save(path + '.png', 'png', quality=95)

    # Start the main event loop
    if disp:
        window.mainloop()


