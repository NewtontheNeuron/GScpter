import json
import os, sys
import tkinter as tk
from tkinter import Tk, ttk

os.chdir(os.path.dirname(sys.argv[0]))

# dictionary that stores all values to be exported in the json file
data = {
    'features': [],
    'project_name': '',
    'clusterpool_names': [],
    'subgroup_names': [],
    'clusterpools': {}
}

# function that exports data to the json file
def export():
    data['project_name'] = projectNameEntry.get()
    for item in featureTreeView.get_children():
        data['features'].append(featureTreeView.item(item)['text'])

    with open('data.json', 'w') as outfile:
        json.dump(data, outfile)


root = tk.Tk()
root.title('Cluster Pool Comparison Tool')

#frame that holds all elements
rootFrame = ttk.Frame(root)
rootFrame.pack(fill="both", expand=True)

#initialize azure theme
root.tk.call("source", "guitheme/azure.tcl")
root.tk.call("set_theme", "light")

##### project name
projectNameEntry = ttk.Entry(rootFrame, width=30)
projectNameEntry.insert(0, 'Project Name')
projectNameEntry.grid(row=0, column=0, columnspan=2, padx=8, pady=8)
#####

##### features
featureFrame = ttk.LabelFrame(rootFrame, text="Features", padding=(20, 10))
featureFrame.grid(row=1, column=0, padx=(20, 10), pady=10)

featureEntry = ttk.Entry(featureFrame, width=20)
featureEntry.insert(0, 'New Feature')
featureEntry.grid(row=0, column=0)

def addFeature():
    feature = featureEntry.get()
    if len(feature) == 0:
        return
    for item in featureTreeView.get_children():
        if feature == featureTreeView.item(item)['text']:
            return
    featureTreeView.insert('', 'end', text=featureEntry.get())

addFeatureButton = ttk.Button(featureFrame, text='Add', width=6, command=addFeature, style='Accent.TButton')
addFeatureButton.grid(row=0, column=1)

def removeFeature():
    for item in featureTreeView.selection():
        featureTreeView.delete(item)

removeFeatureButton = ttk.Button(featureFrame, text='Remove', width=10, command=removeFeature)
removeFeatureButton.grid(row=2, column=1)

featureTreeView = ttk.Treeview(
    featureFrame,
    selectmode="browse",
    height=6)
featureTreeView.grid(row=1, column=0, columnspan=2, pady=10)
featureTreeView.column("#0", anchor="w", width=200)
featureTreeView.heading("#0", text="Features Added", anchor="center")
#####

exportButton = ttk.Button(rootFrame, text='Export JSON', width=10, command=export, style='Accent.TButton')
exportButton.grid(row=2, column=0, pady=10)

root.mainloop()
