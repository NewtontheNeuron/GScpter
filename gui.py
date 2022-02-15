import json
import os, sys
import tkinter as tk
from tkinter import ttk
from functools import partial

class App(ttk.Frame):
    def __init__(self, parent, data):
        ttk.Frame.__init__(self)
        self.parent = parent

        self.setupWidgets()
    
    def setupWidgets(self):
        self.projectNameEntry = ttk.Entry(self, width=30)
        self.projectNameEntry.insert(0, 'Project Name')
        self.projectNameEntry.grid(row=0, column=0, columnspan=2, padx=8, pady=8)

        self.featureFrame = ttk.LabelFrame(self, text="Features", padding=(20, 10))
        self.featureFrame.grid(row=1, column=0, padx=(20, 10), pady=10)

        self.featureEntry = ttk.Entry(self.featureFrame, width=20)
        self.featureEntry.insert(0, "Gene")
        self.featureEntry.grid(row=0, column=0)

        self.addFeatureButton = ttk.Button(self.featureFrame, text='Add', width=6, command=self.addFeature, style='Accent.TButton')
        self.addFeatureButton.grid(row=0, column=1)

        self.removeFeatureButton = ttk.Button(self.featureFrame, text='Remove', width=10, command=self.removeFeature)
        self.removeFeatureButton.grid(row=2, column=1)

        self.featureTreeView = ttk.Treeview(self.featureFrame, selectmode="extended", height=7)
        self.featureTreeView.grid(row=1, column=0, columnspan=2, pady=10)
        self.featureTreeView.column("#0", anchor="w", width=200)
        self.featureTreeView.heading("#0", text="Features Added", anchor="center")

        self.clusterpoolFrame = ttk.LabelFrame(self, text="Clusterpools", padding=(20, 10))
        self.clusterpoolFrame.grid(row=1, column=1, padx=(20, 10), pady=10)

        self.clusterpoolEntry = ttk.Entry(self.clusterpoolFrame, width=20)
        self.clusterpoolEntry.insert(0, "Clusterpool")
        self.clusterpoolEntry.grid(row=0, column=0)

        self.addClusterpoolButton = ttk.Button(self.clusterpoolFrame, text='Add', width=6, command=self.addClusterpool, style='Accent.TButton')
        self.addClusterpoolButton.grid(row=0, column=1)

        self.subgroupEntry = ttk.Entry(self.clusterpoolFrame, width=20)
        self.subgroupEntry.insert(0, "Subgroup")
        self.subgroupEntry.grid(row=1, column=0)

        self.addSubgroupButton = ttk.Button(self.clusterpoolFrame, text='Add', width=6, command=self.addSubgroup, style='Accent.TButton')
        self.addSubgroupButton.grid(row=1, column=1)

        self.clusterpoolTreeView = ttk.Treeview(self.clusterpoolFrame, selectmode="browse", height=6)
        self.clusterpoolTreeView.grid(row=2, column=0, columnspan=2, pady=10)
        self.clusterpoolTreeView.column("#0", anchor="w", width=200)
        self.clusterpoolTreeView.heading("#0", text="Clusterpool Heirarchy", anchor="center")
        self.clusterpoolTreeView.bind('<<TreeviewSelect>>', self.itemSelected)
        
        self.removeGroupButton = ttk.Button(self.clusterpoolFrame, text='Remove', width=10, command=self.removeGroup)
        self.removeGroupButton.grid(row=3, column=0)

        self.editPoolsButton = ttk.Button(self.clusterpoolFrame, text='Edit Pools', width=10, command=self.openPoolPopupWindow, state='disabled')
        self.editPoolsButton.grid(row=3, column=1)

        self.exportButton = ttk.Button(self, text='Export JSON', width=10, command=partial(self.export, data), style='Accent.TButton')
        self.exportButton.grid(row=2, column=0, pady=10, columnspan=2)

    def addFeature(self):
        feature = self.featureEntry.get()
        if len(feature) == 0:
            return
        if not feature.startswith('rna_'):
            feature = 'rna_' + feature
        for iid in self.featureTreeView.get_children():
            if feature == self.featureTreeView.item(iid)['text']:
                return
        self.featureTreeView.insert('', 'end', text=feature)

    def removeFeature(self):
        for item in self.featureTreeView.selection():
            self.featureTreeView.delete(item)

    def addClusterpool(self):
        clusterpool = self.clusterpoolEntry.get()
        if len(clusterpool) == 0:
            return
        for name in data['clusterpool_names']:
            if name == clusterpool: return
        data['clusterpool_names'].append(clusterpool)
        data['clusterpools'][clusterpool] = {}
        id = self.clusterpoolTreeView.insert('', 'end', text=clusterpool)
        for subgroup in data['subgroup_names']:
            self.clusterpoolTreeView.insert(id, 'end', text=subgroup)
            data['clusterpools'][clusterpool][subgroup] = []

    def addSubgroup(self):
        subgroup = self.subgroupEntry.get()
        if len(subgroup) == 0:
            return
        for name in data['subgroup_names']:
            if name == subgroup: return
        data['subgroup_names'].append(subgroup)
        for iid in self.clusterpoolTreeView.get_children():
            self.clusterpoolTreeView.insert(iid, 'end', text=subgroup)
        for clusterpool in data['clusterpools'].keys():
            data['clusterpools'][clusterpool][subgroup] = []
            print(data)

    def removeGroup(self):
        selectedItem = self.clusterpoolTreeView.focus()
        selectedValue = self.clusterpoolTreeView.item(selectedItem)['text']

        if(self.clusterpoolTreeView.parent(self.clusterpoolTreeView.selection()) == ''): #selection is a clusterpool
            data['clusterpool_names'].remove(selectedValue)
            data['clusterpools'].pop(selectedValue)
        else:                                                                            #selection is a subgroup
            data['subgroup_names'].remove(selectedValue)
            for clusterpool in data['clusterpools']:
                data['clusterpools'][clusterpool].pop(selectedValue)

        self.clusterpoolTreeView.delete(selectedItem)
        print(data)
    
    def itemSelected(self, event):
        if(self.clusterpoolTreeView.parent(self.clusterpoolTreeView.selection()) == ''):
            self.editPoolsButton['state'] = 'disabled'
        else:
            self.editPoolsButton['state'] = 'enabled'
        print(self.clusterpoolTreeView.selection())

    def openPoolPopupWindow(self):
        selectedItem = self.clusterpoolTreeView.focus()
        parentItem = self.clusterpoolTreeView.parent(selectedItem)

        clusterpool = self.clusterpoolTreeView.item(parentItem)['text']
        subgroup = self.clusterpoolTreeView.item(selectedItem)['text']
        print("editing", clusterpool, subgroup)

        PoolPopupWindow(self, data['clusterpools'][clusterpool][subgroup])

    def updateData(self):
        data['project_name'] = self.projectNameEntry.get()
        for iid in self.featureTreeView.get_children():
            data['features'].append(self.featureTreeView.item(iid)['text'])

    def export(self, data):
        self.updateData()
        print("--exporting data:")
        print(data)
        with open('Data/JSON/data.json', 'w') as outfile:
            json.dump(data, outfile)

class PoolPopupWindow(tk.Toplevel):
    def __init__(self, parent, poolList):
        tk.Toplevel.__init__(self)
        self.parent = parent
        self.title('Pool Editor')

        self.excCheckboxList = []
        self.excCheckboxVariables = {}
        self.inhCheckboxList = []
        self.inhCheckboxVariables = {}

        for i in range(38):
            name = "Excit-"+str(i+1)
            self.excCheckboxVariables[name] = tk.IntVar(self, 0, name)
            self.excCheckboxList.append(tk.Checkbutton(
                self, 
                text=name,
                variable=self.excCheckboxVariables[name]
            ))
            self.excCheckboxList[i].grid(
                row = i if i<19 else i-19,
                column = 0 if i<19 else 1)
            if name in poolList:
                self.excCheckboxList[i].select()
            else:
                self.excCheckboxList[i].deselect()

        for i in range(27):
            name = "Inhib-"+str(i+1)
            self.inhCheckboxVariables[name] = tk.IntVar(self, 0, name)
            self.inhCheckboxList.append(tk.Checkbutton(
                self, 
                text=name,
                variable=self.inhCheckboxVariables[name]
            ))
            self.inhCheckboxList[i].grid(
                row = i if i<19 else i-19,
                column = 2 if i<19 else 3
            )
            if name in poolList:
                self.inhCheckboxList[i].select()
            else:
                self.inhCheckboxList[i].deselect()

        copyPoolsButton = ttk.Button(self, text="Copy Pools from...", command=partial(self.openCopyPoolsWindow) )
        copyPoolsButton.grid(row=19, column=2)

        closeButton = ttk.Button(self, text="Close", command=partial(self.savePoolListAndClose, poolList) )
        closeButton.grid(row=19, column=3)

    def openCopyPoolsWindow(self):
        CopyPoolsWindow(self)

    def savePoolListAndClose(self, poolList):
        poolList.clear()
        for key in self.excCheckboxVariables:
            # print(self.excCheckboxVariables[key])
            if self.excCheckboxVariables[key].get() == 1:
                poolList.append(key)
        for key in self.inhCheckboxVariables:
            # print(self.inhCheckboxVariables[key])
            if self.inhCheckboxVariables[key].get() == 1:
                poolList.append(key)
        print(poolList)
        self.destroy()

class CopyPoolsWindow(tk.Toplevel):
    def __init__(self, parent):
        tk.Toplevel.__init__(self)
        self.parent = parent
        self.title('Copy Pools')

        self.copyPoolsLabel = ttk.Label(self, text="Copy pools from:")
        self.copyPoolsLabel.grid(row=0,column=0)

        self.copyPoolsTreeView = ttk.Treeview(self, selectmode="browse", height=6)
        self.copyPoolsTreeView.grid(row=1, column=0, columnspan=2, pady=10)
        self.copyPoolsTreeView.column("#0", anchor="w", width=200)
        self.copyPoolsTreeView.heading("#0", text="Clusterpool Options", anchor="center")

        for clusterpool in data['clusterpool_names']:
            for subgroup in data['subgroup_names']:
                name = clusterpool + ' // ' + subgroup
                self.copyPoolsTreeView.insert('', 'end', text=name)
        
        self.selectButton = ttk.Button(self, text="Select", command=self.copyPools)
        self.selectButton.grid(row=2, column=1)

    def copyPools(self):
        selectedItem = self.copyPoolsTreeView.focus()
        selectedValue = self.copyPoolsTreeView.item(selectedItem)['text']
        selectedValue = selectedValue.split(' // ')
        clusterpool = selectedValue[0]
        subgroup = selectedValue[1]

        poolList = data['clusterpools'][clusterpool][subgroup]

        for i in range(38):
            name = "Excit-"+str(i+1)
            if name in poolList:
                self.parent.excCheckboxList[i].select()
            else:
                self.parent.excCheckboxList[i].deselect()

        for i in range(27):
            name = "Inhib-"+str(i+1)
            if name in poolList:
                self.parent.inhCheckboxList[i].select()
            else:
                self.parent.inhCheckboxList[i].deselect()

        self.destroy()





if __name__ == "__main__":
    data = {
        'features': [],
        'project_name': '',
        'clusterpool_names': [],
        'subgroup_names': [],
        'clusterpools': {}
    }

    os.chdir(os.path.dirname(sys.argv[0]))
    
    root = tk.Tk()
    root.title("Cluster Pool Comparison Tool")

    #initialize azure theme
    root.tk.call("source", "guitheme/azure.tcl")
    root.tk.call("set_theme", "light")

    app = App(root, data)
    app.pack(fill="both", expand=True)

    root.mainloop()