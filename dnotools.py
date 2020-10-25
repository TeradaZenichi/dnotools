# Institution: University of Campinas
# Author: Lucas Zenichi Terada

from gui.interface import *
from data import config

class application(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.geometry("800x700")
        container = tk.Frame(self)
        container.config(bg='white')
        container.pack(side="top", fill="both", expand="True")
        self.frames = {}
        for F in (startpage, powerfpage,restpage,recpage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.config(bg='white')
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(startpage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


def main():
    app = application()
    app.mainloop()
    print("Vou imprimir:\n")
    print(config.bus)

if __name__ == "__main__":
    main()
    print(config.bus)
