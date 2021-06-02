# VISUALISER OF CONTAMINATION
# v0.0.1, 21.8.2020

import tkinter as tk
from tkinter import messagebox as msgb
from tkinter import filedialog as fdialog
from itertools import product
import pickle

from src.ppopt import *

# ----------------------WINDOW (AS A CLASS DEFINITION)---------------------
class Vis(tk.Frame):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.choose_and_load()


        # display tutorial messagebox
        title='Welcome to contamination visualiser tool'
        message='Welcome!\nHere is how to use the contamination visulaiser tool:\n\n'
        message+='Left-click a well - for all the tips that served it, the parts that were present in the previously-visited wells will be displayed in bold.'
        message+='\n\nMiddle-click - undo left-click'
        message+='\n\nRight-click - display well contents in a message box'
        msgb.showinfo(title=title, message=message)

        self.master.title("Construct contamination")
        self.pack(fill=tk.BOTH, expand=1)

        self.canvas = tk.Canvas(self)
        self.mode = 'initial'
        self.coord_w = {(i, j): -1 for i, j in product(range(1, 9), range(1, 13))}
        self.rows = 'ABCDEFGH'
        self.rowcoords = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8}
        self.defx_tl = 30
        self.defy_tl = 30
        self.xside = 115
        self.yside = 10 * len(self.w[0]) + 20
        self.interval = 10

        self.basicpaint()

        self.mode = 'general'

    def basicpaint(self):
        self.canvas.delete("all")

        if (self.mode != 'clicked'):
            self.bold = []
            empty1 = []
            for i in range(0, len(self.w[0])):
                empty1.append(False)
            for i in range(0, len(self.w)):
                self.bold.append(empty1.copy())

        for i in range(1, 9):
            for j in range(1, 13):
                y_tl = self.defy_tl + (i - 1) * (self.yside + self.interval)
                x_tl = self.defx_tl + (j - 1) * (self.xside + self.interval)
                self.canvas.create_rectangle(x_tl, y_tl, x_tl + self.xside, y_tl + self.yside,
                                             outline="#000")

        # label self.rows and columns
        for i in range(1, 9):
            y_tl = self.defy_tl + self.yside / 2 + (i - 1) * (self.yside + self.interval)
            self.canvas.create_text(5, y_tl, anchor=tk.W, font="Verdana",
                                    text=self.rows[i - 1])
        for j in range(1, 13):
            x_tl = self.defx_tl + self.xside / 3 + (j - 1) * (self.xside + self.interval)
            self.canvas.create_text(x_tl, 10, anchor=tk.W, font="Verdana",
                                    text=str(j))
        self.canvas.pack(fill=tk.BOTH, expand=1)

        # display constructs
        for i in range(0, len(self.w)):
            con = self.dic['constructs'][i]
            if (self.assembly == 'Start-Stop'):
                precoord = con['con_liqloc'][0].display_name
                x_coord_str = ''
                for k in range(1, len(precoord)):
                    if (precoord[k] != ' '):
                        x_coord_str += precoord[k]
                    else:
                        break
                x_coord = int(x_coord_str)
                y_coord = self.rowcoords[precoord[0]]
            elif (self.assembly == 'BASIC'):
                precoord = con['con_liqloc']
                x_coord = int(precoord[1:])
                y_coord = self.rowcoords[precoord[0]]

            y_tl = self.defy_tl + (y_coord - 1) * (self.yside + self.interval)
            x_tl = self.defx_tl + (x_coord - 1) * (self.xside + self.interval)
            self.canvas.create_text(x_tl + 2, y_tl + 10, anchor=tk.W, font=("Verdana", 6),
                                    text=con['con_name'])
            for j in range(0, len(self.w[i])):
                part_name = self.dic['parts'][self.w[i][j]]['part_name']
                if (self.bold[i][j] == True):
                    self.canvas.create_text(x_tl + 20, y_tl + 20 + j * 10, anchor=tk.W, font=('Verdana', 6, 'bold'),
                                            text=part_name, )
                else:
                    self.canvas.create_text(x_tl + 20, y_tl + 20 + j * 10, anchor=tk.W, font=('Verdana', 6),
                                            text=part_name, )

            if (self.mode == 'initial'):
                self.coord_w[(x_coord, y_coord)] = i
            elif (self.mode == 'clicked'):
                if (self.select == i):
                    self.canvas.create_rectangle(x_tl, y_tl, x_tl + self.xside, y_tl + self.yside, outline="#f00")
        # if(self.mode == 'initial')
        # self.tutorial()

    def click(self, event):
        n = self.clickedwell(event.x, event.y)
        if (n != -1):
            self.mode = 'clicked'
            self.bold_cont(n)
            self.select = n
            self.basicpaint()

    def bold_cont(self, n):
        self.bold = []
        empty1 = []
        for i in range(0, len(self.w[0])):
            empty1.append(False)
        for i in range(0, len(self.w)):
            self.bold.append(empty1.copy())

        added = np.zeros((len(self.w), len(self.w[0])))  # tells which parts were added to which well

        # for the first operation in fin
        added[self.fin[0].well][self.fin[0].part[0]] = 1
        for i in range(1, len(self.fin)):
            one_cost = cost_func_with_w(self.fin[0:i], self.fin[i], self.w, added, self.caps)
            added[self.fin[i].well][self.fin[i].part[0]] = 1

            if (self.fin[i].well == n and self.fin[i].changed != True):
                backroll = 1
                while (True):
                    previ = i - backroll
                    for j in range(0, len(self.w[self.fin[previ].well])):
                        if (added[self.fin[previ].well][j] == 1):
                            self.bold[self.fin[previ].well][j] = True
                    if (self.fin[previ].changed == 1):
                        break
                    else:
                        backroll += 1


    def clickedwell(self, x, y):
        x -= 25
        y -= 25
        x_coord = int(np.ceil(x / (self.xside + self.interval)))
        y_coord = int(np.ceil(y / (self.yside + self.interval)))
        return self.coord_w[(x_coord, y_coord)]

    def unclick(self, event):
        self.mode = 'general'
        self.basicpaint()

    def examine(self, event):
        n = self.clickedwell(event.x, event.y)

        title = 'Construct ' + self.dic['constructs'][n]['con_name']

        # for Start-Stop Assembly, we exactly know the part types
        if (self.assembly == 'Start-Stop'):
            part_type = ['Promoter: ', 'RBS: ', 'CDS: ', 'Terminator: ', 'Backbone: ']
        # otherwise, create generic insert names
        else:
            part_type = ['DNA Part 1: ']
            for i in range(1, len(self.w[n])):
                part_type += ['DNA Part ' + str(i + 1) + ': ']
        message = ''
        if (self.mode != 'clicked'):
            for i in range(0, len(self.w[n])):
                message += part_type[i] + self.dic['parts'][self.w[n][i]]['part_name']
                message += '                                    \n'
        else:
            if (n != self.select):
                present = 'Parts present before tip proceeded to selected well:\n'
                notpresent = 'Parts not present before tip proceeded to selected well:\n'
                for i in range(0, len(self.w[n])):
                    if (self.bold[n][i] == True):
                        present += part_type[i] + self.dic['parts'][self.w[n][i]]['part_name']
                        present += '\n'
                    else:
                        notpresent += part_type[i] + self.dic['parts'][self.w[n][i]]['part_name']
                        notpresent += '\n'
                message = present + '\n' + notpresent
            else:
                message = 'This is the well whose predecessors are displayed'

        msgb.showinfo(title=title, message=message)

    def tutorial(self):
        title = 'Welcome to contamination viewer'
        message = 'Left-click a well with construct to check which wells were visited befroe it\n'
        message += 'and which parts were present in these wells then\n\n'
        message += 'Middle-click anywhere to return to general view\n\n'
        message += 'Right-click a well to inspect it in detail'
        msgb.showinfo(title=title, message=message)

    def choose_and_load(self):
        # select the recording
        output = fdialog.askopenfile(
            title='Select the recording', filetypes=(("Pickle files", "*.p"),
                                         ("all files", "*.*")))

        # load the recording file
        rec = pickle.load(open(output.name, 'rb'))
        self.assembly = rec['assembly']
        self.w = rec['w']
        self.fin = rec['fin']
        self.dic = rec['dic']
        self.caps = rec['caps']

        # BASIC assembly does not have inherent names of parts or constructs stored,
        # as parts are PREPARED IN SITU on a magnetic bead well plate.
        # Thus we refer to these parts by their location on magbead plate only.
        if (rec['assembly'] == 'BASIC'):
            for part in self.dic['parts'].values():
                part['part_name'] = 'Magbead well ' + part['part_liqloc']
            for construct in self.dic['constructs'].values():
                construct['con_name'] = 'Well ' + construct['con_liqloc']


# ---------------------RECORDER OF INFORMATION---------------------
def rec(assem, w, fin, dic,caps):
    recording={'assembly': assem, 'w': w, 'fin': fin, 'dic': dic, 'caps': caps}

    if (assem=='Start-Stop'):
        name = 'ppopt_recording_StartStop.p'
    elif(assem=='BASIC'):
        name = 'ppopt_recording_BASIC_'+ str(len(w[0])) + 'parts.p'

    pickle.dump(recording, open(name, 'wb'))
    return


# ---------------------- CALL THE VISUALISER --------------------
def main():
    root = tk.Tk()
    vis = Vis()
    root.geometry("1600x900")

    root.bind("<Button-1>", vis.click)
    root.bind("<Button-2>", vis.unclick)
    root.bind("<Button-3>", vis.examine)

    root.mainloop()


if __name__ == "__main__":
    main()
