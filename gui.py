import json
import tkinter as tk
from tkinter import Tk, ttk

# dictionary that stores all values to be exported in the json file
data = {
    'features': [],
    'project_name': '',
    'filenames': {},
    'clusterpools': {}
}

# function that exports data to the json file
def export():
    data['project_name'] = projectNameEntry.get()

    with open('data.json', 'w') as outfile:
        json.dump(data, outfile)


# gui to edit data
root = tk.Tk()
root.title('Cluster Pool Comparison Tool')

# style = ttk.Style(root)
# print(style.theme_names())
# style.theme_use('clam')

tk.Label(root, text='Project Name').grid(row=0, column=0)
projectNameEntry = tk.Entry(root)
projectNameEntry.grid(row=0, column=1)

exportButton = ttk.Button(root, text='Export JSON', width=15, command=export)
exportButton.grid(row=1, column=1, padx=8, pady=8)

# style.configure('Modified.TButton', background='#ffffff')
# exportButton.config(style = "Modified.TButton")

root.mainloop()
