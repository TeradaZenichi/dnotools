# Author: Lucas Zenichi Terada
# Institution: University of Campinas

from PIL import Image, ImageTk
from gui import tools
from data import config
import tkinter as tk
import tkinter.filedialog as fd
import pandas as pd
import sys


class startpage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        render = ImageTk.PhotoImage(Image.open('gui/logo.png'))
        logo = tk.Label(self, image=render, background='white')
        logo.image = render
        label    = tk.Label(self,text='Python Energy Distribution \n Solutions',**tools.title())
        buttonpf = tk.Button(self, **tools.bgreen('PowerFlow'),command=lambda:controller.show_frame(powerfpage))
        buttonrt = tk.Button(self, **tools.bgreen('Restoration'),command=lambda: controller.show_frame(recpage))
        buttonrc = tk.Button(self, **tools.bgreen('Reconfiguration'),command=lambda: controller.show_frame(restpage))
        cancelbt = tk.Button(self, **tools.bred('Cancel'),command=lambda:sys.exit(0))

        labelpfw = tk.Label(self, text="Solve Power Flow Problems", **tools.label())
        labelrtn = tk.Label(self, text="Solve Restoration Problems", **tools.label())
        labelrce = tk.Label(self, text="Solve Reconfiguration Problems", **tools.label())

        logo.grid(row=1, column=1, padx=20, pady=20, sticky='n')
        label.grid(row=1, column=2, padx=20, pady=20, sticky='n')
        buttonpf.grid(row=2, column=1, padx=20, pady=20, sticky='n')
        labelpfw.grid(row=2, column=2, padx=20, pady=20, sticky='w')
        buttonrt.grid(row=3, column=1, padx=20, pady=20, sticky='n')
        labelrtn.grid(row=3, column=2, padx=20, pady=20, sticky='w')
        buttonrc.grid(row=4, column=1, padx=20, pady=20, sticky='n')
        labelrce.grid(row=4, column=2, padx=20, pady=20, sticky='w')
        cancelbt.grid(row=5, column=1, padx=20, pady=20, sticky='n')

class powerfpage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        buses  = tk.Button(self, **tools.blblue('Bus Data'),command=lambda:self.load_file())
        show   = tk.Button(self, **tools.blblue('Bus Data'),command=lambda:self.show())
        backbt = tk.Button(self, **tools.blblue('Back'),command=lambda:controller.show_frame(startpage))

        buses.grid(row=4, column=1, padx=20, pady=20, sticky='n')
        show.grid(row=5, column=1, padx=20, pady=20, sticky='n')
        backbt.grid(row=6, column=1, padx=20, pady=20, sticky='n')

    def load_file(self):
        fname = fd.askopenfilename()
        # pedssystems.bus = pd.read_csv(fname)
        # print(pedssystems.bus)
        config.bus = pd.read_csv(fname)
        print(config.bus)
        return

    def show(self):
        print(config.bus)
        return


class restpage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        button1 = tk.Button(self, text="back to home",
                            command=lambda: controller.show_frame(startpage))
        button1.pack()

class recpage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        button1 = tk.Button(self, text="back to home",
                            command=lambda: controller.show_frame(startpage))
        button1.pack()
