import pathlib, sys, os
SCRIPT_DIR = str(pathlib.Path().absolute())
sys.path.insert(0,SCRIPT_DIR)
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from functools import partial

def tkcolor(rgb, alpha=1):
    return "#" + ''.join(f"{int(255 - alpha * (255 - c)):02x}" for c in rgb)

class DataTk():
    def raise_popups(self):
        for w in self.root.winfo_children():
            if w.winfo_class() ==  "Toplevel":
                w.tkraise()
    def b_a_w(self):
        if self.app.BaW:
            self.app.BaW = False
            self.button_baw.config(text="black/white")
        else:
            self.app.BaW = True
            self.button_baw.config(text="color")
        self.plot_all()
    def __init__(self, root, app, initialdir=pathlib.Path().absolute()):
        # initialdir = pathlib.Path.home().joinpath('data_dir')
        self.initialdir = initialdir
        root.wm_title(app.__class__.__name__)
        self.root = root
        self.app = app
        self.plotters = {p.__name__:p for p in app.def_plotters()}
        nb = ttk.Notebook(self.root)
        self.nb = nb
        color = tkcolor((160, 160, 160), 0.3)
        self.nbframes, self.canvases, self.figs = {}, {}, {}
        self.dir_selected = ''
        frame_label = tk.Frame(root, bg=color)
        self.text_dir_selected = tk.StringVar()
        label = ttk.Label(frame_label, textvariable=self.text_dir_selected, background='yellow')
        label.pack(side="left", expand=True, fill=tk.X)
        frame_label.pack(fill=tk.X, expand="no", padx=5, pady=5)
        frame_button = tk.Frame(root, bg=color)
        button = ttk.Button(frame_button, text="open (dir)", command=self.open)
        button.pack(side="left", expand=True, fill=tk.X)
        self.app.BaW = False
        button = ttk.Button(frame_button, text="black/white", command=self.b_a_w)
        button.pack(side="left", expand=True, fill=tk.X)
        self.button_baw = button
        frame_button.pack(fill=tk.X, expand="no", padx=5, pady=5)
        for k,plotter in self.plotters.items():
            self.nbframes[k] = tk.Frame(nb, bg=color)
            self.nbframes[k].pack(fill="both", expand="yes", padx=5, pady=5)
            nb.add(self.nbframes[k], text=k)
            # plot notebook
            self.figs[k] = Figure(figsize=(5, 4), dpi=100)
            self.canvases[k] = FigureCanvasTkAgg(self.figs[k], master=self.nbframes[k])  # A tk.DrawingArea.
            canvas = self.canvases[k]
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        nb.pack(fill="both", expand="yes", padx=5, pady=5)
        # menubar
        menubar = tk.Menu(self.root)
        menu_file = tk.Menu(menubar, tearoff=0)
        menu_file.add_command(label="Quit", command=self.quit, accelerator="Command-q")
        menu_file.add_command(label="Save single", command=self.save_single)
        menu_file.add_command(label="Save all", command=self.save, accelerator="Command-s")
        menu_file.add_command(label="Raise PupUps", command=self.raise_popups, accelerator="Command-R")
        menubar.add_cascade(label="File", menu=menu_file)
        menu_plot = tk.Menu(menubar, tearoff=0)
        self.plotters_checked = {}
        for k,plotter in self.plotters.items():
            print(f"{k=} {plotter=}")
            item = tk.IntVar()
            item.set(1)
            self.plotters_checked[k]=item
            menu_plot.add_checkbutton(label=f"{k}", variable=item)
        menubar.add_cascade(label="Plot", menu=menu_plot)
        self.root.config(menu=menubar)

    def quit(self):
        # self.root.destroy()
        self.root.quit()

    def save_single(self):
        if not hasattr(self, 'inderem'): return
        index = self.nb.index(self.nb.select())
        k = list(self.plotters.keys())[index]
        print(f"{index=} {k=}")
        filename = filedialog.asksaveasfilename(initialdir=self.initialdir, initialfile='_'.join(self.inderem))
        filename += k + ".png"
        self.canvases[k].figure.savefig(pathlib.Path(filename))
    def save(self):
        if not hasattr(self, 'inderem'): return
        # dirname = filedialog.asksaveasfilename(initialdir=self.initialdir, initialfile='_'.join(self.inderem))
        dirname = filedialog.asksaveasfilename(initialdir=self.initialdir)
        dirname += '_'.join(self.inderem)
        # print(f"{dirname=}")
        dir = pathlib.Path(dirname)
        dir.mkdir()
        # print(f"{dir=}")
        for k,plotter in self.plotters.items():
            filename = k + ".png"
            self.canvases[k].figure.savefig(dir/filename)
    def open(self, single=False):
        # print(f"{self.initialdir=}")
        dirname = filedialog.askdirectory(initialdir=self.initialdir)
        self.dir_selected = pathlib.Path(dirname)
        self.text_dir_selected.set(str(self.dir_selected))
        self.app.init_from_directory(self.dir_selected)
        self.plot_all()
    def plot_all(self):
        for k,plotter in self.plotters.items():
            if not self.plotters_checked[k].get(): continue
            fig, canvas = self.canvases[k].figure, self.canvases[k]
            fig.clear()
            plotter(fig)
            canvas.draw()
        self.root.focus_force()


#=================================================================#
if __name__ == '__main__':
    class SmallTest:
        def def_plotters(self):
            return [self.plot1, self.plot2]
        def init_from_directory(self, dir):
            print(f"{dir=}")
        def plot1(self, fig):
            print(f"hallo plot1")
        def plot2(self, fig):
            print(f"hallo plot2")
    root = tk.Tk()
    dtk = DataTk(root=root, app=SmallTest())
    dtk.open()
    tk.mainloop()