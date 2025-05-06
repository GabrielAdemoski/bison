import tkinter as tk
import numpy as np
from PIL import Image
import os


def circle(canvas, x, y, rad, **kwargs):
    """
    Creates a tkinter oval object in the shape of a circle with its center at x,y and with a radius of rad.
    :param canvas: the canvas object on which the circle will be drawn
    :param x: x coordinate for the center of the circle
    :param y: y coordinate for the center of the circle
    :param rad: the radius of the resultant circle
    :param kwargs: tkinter create_oval arguments
    :return: A tkinter oval object
    """
    return canvas.create_oval(x - rad, y - rad, x + rad, y + rad, kwargs)


def from_to(matrix):
    """
    Creates tuples containing two nodes and the edge value between them. Only creates tuples for nodes with a nonzero
    edge. All edges have been normalized to between -1 and 1. The first node is the row index and the second node is the
    column index.
    :param matrix: A 2D numpy array corresponding to an adjacency matrix
    :return: A list of tuples containing two nodes and the edge value between them
    """
    from_to = []
    r = 0
    if np.any(matrix):
        # normalize to between -1 and 1
        matrix /= np.max(abs(matrix[np.nonzero(matrix)]))
        for row in matrix:
            c = 0
            for val in row:
                if val != 0:
                    from_to.append((r, c, val))
                c += 1
            r += 1
    return from_to


def quadrant(theta):
    """
    Returns the quadrant number based on the angle provided.
    :param theta: An angle provided in radians
    :return: q, the quadrant number
    """
    q = 0
    theta = theta % (2 * np.pi)
    if theta >= 0 and theta < np.pi / 2:
        q = 1
    elif theta >= np.pi / 2 and theta < np.pi:
        q = 2
    elif theta >= np.pi and theta < (3 * np.pi / 2):
        q = 3
    else:
        q = 4
    return q


def anchor(q):
    """
    Returns the anchoring position of the text box based on the quadrant it is in
    :param q: The quadrant the text box is in
    :return: The anchor position of the text box
    """
    match q:
        case 0:
            a = "center"
        case 1:
            # Bottom center
            a = "s"
        case 2:
            # Left center
            a = "w"
        case 3:
            # Top center
            a = "n"
        case 4:
            # Right center
            a = "e"
    return a


def circular_plot(matrix, taxa, disp=True, path="figs/", fname='circular_image', title='X'):
    os.makedirs(path, exist_ok=True)

    width = 800
    height = width
    border = 150
    wc = width / 2
    hc = height / 2
    inner_rad = 50
    max_rad = min(wc, hc)
    outer_rad = max_rad - border
    line_length = outer_rad - inner_rad
    nsamples = matrix.shape[0]
    nlines = nsamples + 1
    n_nodes = matrix.shape[1]
    node_offset = 0

    txt = ["Fall", "Winter", "Spring", "Summer"]
    OTUs = [i for i in range(n_nodes)]
    # blue is negative and orange is positive
    line_color = ["blue", "orange"]
    tick_offset = 5

    # In degrees
    deg_offset = 0
    # In radians because all trig functions use radians
    rad_offset = np.radians(deg_offset)
    thetas = np.linspace(0, 2 * np.pi, nlines)
    lable_offset = 15

    dist = np.linspace(inner_rad + node_offset, outer_rad - node_offset, n_nodes)

    # Create the main window
    window = tk.Tk()
    window.title(fname+" Canvas")
    # Create a canvas widget
    canvas = tk.Canvas(window, width=width, height=height, bg="white")
    canvas.pack()
    draw_otu = False
    # Canvas title, centered with a downwards offset
    canvas.create_text(width / 2, border / 3, text=title, font=('Helvetica 25 bold'))

    for theta in thetas[:-1]:
        # calculate the start and end vertices of the base lines
        ix = -inner_rad * np.sin(-theta - rad_offset) + wc
        iy = -inner_rad * np.cos(-theta - rad_offset) + hc
        ox = -outer_rad * np.sin(-theta - rad_offset) + wc
        oy = -outer_rad * np.cos(-theta - rad_offset) + hc

        canvas.create_line((ix, iy), (ox, oy), fill='black', width=3)
        # Arrays of x and y coordinates for the current line
        xs = -dist * np.sin(-theta - rad_offset) + wc
        ys = -dist * np.cos(-theta - rad_offset) + hc

        idx = list(thetas).index(theta)
        # Arrays of x and y coordinates for the next line
        xs1 = -dist * np.sin(-(thetas[(idx + 1) % nsamples]) - rad_offset) + wc
        ys1 = -dist * np.cos(-(thetas[(idx + 1) % nsamples]) - rad_offset) + hc

        # Draw category text starting at the top and moving clockwise
        q = quadrant(theta)
        a = anchor(q)
        textx = -(outer_rad + (lable_offset)) * np.sin(-theta - rad_offset) + wc
        texty = -(outer_rad + (lable_offset)) * np.cos(-theta - rad_offset) + wc

        canvas.create_text(textx, texty, text=txt[idx % len(txt)], fill="black", font=('Helvetica 25 bold'), anchor=a)

        # Draw lines representing graph edges
        for m in from_to(matrix[idx]):
            weight = m[-1]
            line = canvas.create_line((xs[m[0]], ys[m[0]]), (xs1[m[1]], ys1[m[1]]), fill=line_color[int(weight > 0)],
                                      width=3*abs(weight))#width=2 * weight**2)#width=(abs(weight) * 4) + 1)
            canvas.lower(line)

        for i in range(n_nodes):
            # Draw nodes on their respective lines starting from the center and moving outwards
            # circle(canvas, xs[i], ys[i], 2, fill="black")
            # canvas.create_line((xs[i], ys[i] + tick_offset), (xs[i], ys[i] - tick_offset), fill='black', width=2)

            if draw_otu == True:
                # Draw node numbering text starting from the center and moving outwards
                # No rotation, text left to right with a downwards offset
                if (theta % np.radians(180)) >= np.radians(0) and (theta % np.radians(225)) < np.radians(45):
                    # canvas.create_line((xs[i], ys[i] + tick_offset), (xs[i], ys[i] - tick_offset), fill='black',
                    #                    width=2)
                    canvas.create_text(xs[i] + tick_offset, ys[i], text=taxa[i], font=("Helvetica 8 bold"))
                elif (theta % np.radians(225)) >= np.radians(45) and (theta % np.radians(270)) < np.radians(90):
                    # canvas.create_line((xs[i] + tick_offset, ys[i]), (xs[i] - tick_offset, ys[i]), fill='black',
                    #                    width=2)
                    canvas.create_text(xs[i], ys[i] + tick_offset, text=taxa[i], font=("Helvetica 8 bold"))
                elif (theta % np.radians(270)) >= np.radians(90) and (theta % np.radians(315)) < np.radians(135):
                    # canvas.create_line((xs[i] + tick_offset, ys[i]), (xs[i] - tick_offset, ys[i]), fill='black',
                    #                    width=2)
                    canvas.create_text(xs[i], ys[i] + tick_offset, text=taxa[i], font=("Helvetica 8 bold"))
                else:
                    # canvas.create_line((xs[i], ys[i] + tick_offset), (xs[i], ys[i] - tick_offset), fill='black',
                    #                    width=2)
                    canvas.create_text(xs[i] + tick_offset, ys[i], text=taxa[i], font=("Helvetica 8 bold"))

    # Coords for an arrow showing the direction of cycle flow
    arrowint = 25
    arrowscale = .5
    arrowx = -arrowscale * inner_rad * np.sin(np.linspace(0, 2 * np.pi, arrowint) - rad_offset) + wc
    arrowy = -arrowscale * inner_rad * np.cos(np.linspace(0, 2 * np.pi, arrowint) - rad_offset) + hc

    # Draw the central arrow
    points = [(arrowx[-n], arrowy[-n]) for n in range(int(arrowint * .75))]
    canvas.create_line(points, arrow='last', smooth=1, width=1)
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

