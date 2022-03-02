import json
import os, sys
import tkinter as tk                # python GUI package
from tkinter import ttk             # allows for the use of themed widgets (in this app they use Azure, stored in guitheme)
from functools import partial


'''
App class
- the main GUI window
'''
class App(ttk.Frame):
    '''
    init
    - function called when an instance of the class is created
    '''
    def __init__(self, parent, data):
        ttk.Frame.__init__(self)
        self.parent = parent

        self.setupWidgets()
    
    '''
    setupWidgets
    - function that creates and places the widgets within the app window
    '''
    def setupWidgets(self):
        self.projectNameEntry = ttk.Entry(self, width=30)
        self.projectNameEntry.insert(0, 'Project Name')
        self.projectNameEntry.grid(row=0, column=0, columnspan=2, padx=8, pady=8)

        '''
        featureFrame
        - contains all widgets related to feature selection (i.e. adding genes to the analysis)
        '''
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

        '''
        clusterpoolFrame
        - contains all widgets related to clusterpool selection (i.e. defining pools and subpools)
        '''
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
        
        self.removeGroupButton = ttk.Button(self.clusterpoolFrame, text='Remove', width=10, command=self.removePool)
        self.removeGroupButton.grid(row=3, column=0)

        self.editPoolsButton = ttk.Button(self.clusterpoolFrame, text='Edit Pools', width=10, command=self.openPoolPopupWindow, state='disabled')
        self.editPoolsButton.grid(row=3, column=1)

        self.exportButton = ttk.Button(self, text='Export JSON', width=10, command=partial(self.export, data), style='Accent.TButton')
        self.exportButton.grid(row=2, column=0, pady=10, columnspan=2)

    '''
    addFeature
    - adds the feature in featureEntry to the featureTreeView (the list of features)
    '''
    def addFeature(self):
        feature = self.featureEntry.get()
        if len(feature) == 0:                                           # featureEntry is empty -- return
            return
        if not feature.startswith('rna_'):                              # add 'rna_' prefix if it does not already have one
            feature = 'rna_' + feature
        for iid in self.featureTreeView.get_children():                 # search featureTreeView for a gene with the same name
            if feature == self.featureTreeView.item(iid)['text']:       # if it is a duplicate, don't add the feature
                return
        self.featureTreeView.insert('', 'end', text=feature)            # otherwise add the feature

    '''
    removeFeature
    - removes the selected features from the featureTreeView (the list of features)
    '''
    def removeFeature(self):
        for item in self.featureTreeView.selection():                   # get the selected items as a list
            self.featureTreeView.delete(item)                           # delete each item

    '''
    addClusterpool
    - adds a clusterpool with the name contained in clusterpoolEntry
    - clusterpools are at the top of the heirarchy, and subgroups are consistent across all clusterpools
    '''
    def addClusterpool(self):
        clusterpool = self.clusterpoolEntry.get()
        if len(clusterpool) == 0:                                               # clusterpoolEntry is empty -- return
            return
        for name in data['clusterpool_names']:                                  # if a clusterpool with this name already exists, return
            if name == clusterpool: return
        data['clusterpool_names'].append(clusterpool)                           # add the clusterpool to data
        data['clusterpools'][clusterpool] = {}
        iid = self.clusterpoolTreeView.insert('', 'end', text=clusterpool)      # add the clusterpool to clusterpoolTreeView so that it displays to the user
        for subgroup in data['subgroup_names']:                                 # add all subgroups to the clusterpool in both data and the Treeview
            self.clusterpoolTreeView.insert(iid, 'end', text=subgroup)
            data['clusterpools'][clusterpool][subgroup] = []

    '''
    addSubgroup
    - adds a subgroup with the name contained in subgroupEntry
    - subgroups are smaller pools that divide the contents of a clusterpool, and are consistent across all clusterpools
    '''
    def addSubgroup(self):
        subgroup = self.subgroupEntry.get()
        if len(subgroup) == 0:                                                  # subgroupEntry is empty -- return
            return
        for name in data['subgroup_names']:                                     # if a subgroup with this name already exists, return
            if name == subgroup: return
        data['subgroup_names'].append(subgroup)                                 # add subgroup name to data
        for iid in self.clusterpoolTreeView.get_children():                     # add the subgroup to every clusterpool in the Treeview
            self.clusterpoolTreeView.insert(iid, 'end', text=subgroup)
        for clusterpool in data['clusterpools'].keys():                         # add the subgroup to every clusterpool in data
            data['clusterpools'][clusterpool][subgroup] = []
            print(data)

    '''
    removePool
    - removes the selected clusterpool or subgroup
    - if the selection is a clusterpool, the subgroups it contains are not deleted because they are consistent across groups
    - if the selection is a subgroup, this subgroup is removed from all clusterpools
    '''
    def removePool(self):
        selectedItem = self.clusterpoolTreeView.focus()                                     # get selected item and its value (i.e. the pool name)
        selectedValue = self.clusterpoolTreeView.item(selectedItem)['text']

        if(self.clusterpoolTreeView.parent(self.clusterpoolTreeView.selection()) == ''):    # if the pool is at the top of the heirarchy, it is a clusterpool
            data['clusterpool_names'].remove(selectedValue)                                     # remove the clusterpool from data and the Treeview
            data['clusterpools'].pop(selectedValue)
            self.clusterpoolTreeView.delete(selectedItem)
        else:                                                                               # otherwise the pool is a subgroup
            data['subgroup_names'].remove(selectedValue)                                        # remove the subgroup from every clusterpool in data
            for clusterpool in data['clusterpools']:
                data['clusterpools'][clusterpool].pop(selectedValue)
            for clusterpool in self.clusterpoolTreeView.get_children():                         # remove the subgroup from every clusterpool in the Treeview
                for subgroup in self.clusterpoolTreeView.get_children(clusterpool):
                    if(self.clusterpoolTreeView.item(subgroup)['text'] == selectedValue):
                        self.clusterpoolTreeView.delete(subgroup)

        print(data)
    
    '''
    itemSelected
    - event handling function -- called every time the user changes the selection in clusterpoolTreeView
    - enables and disables the editPoolsButton depending on whether the user is selecting a clusterpool or subgroup
    - the editPoolsButton should only be enabled if the user is selecting a subgroup, because clusterpools cannot directly contain clusters
    '''
    def itemSelected(self, event):
        if(self.clusterpoolTreeView.parent(self.clusterpoolTreeView.selection()) == ''):    # if the user is selecting a clusterpool
            self.editPoolsButton['state'] = 'disabled'
        else:                                                                               # if the user is selecting a subgroup
            self.editPoolsButton['state'] = 'enabled'
        print(self.clusterpoolTreeView.selection())

    '''
    openPoolPopupWindow
    - opens a popup window that allows the user to edit the clusters included in the selected subgroup
    - this function is called when the user clicks the editPoolsButton
    '''
    def openPoolPopupWindow(self):
        selectedItem = self.clusterpoolTreeView.focus()                             # get selected subgroup as an item
        parentItem = self.clusterpoolTreeView.parent(selectedItem)                  # get the clusterpool that the selected subgroup is in as an item

        clusterpool = self.clusterpoolTreeView.item(parentItem)['text']
        subgroup = self.clusterpoolTreeView.item(selectedItem)['text']
        print("editing", clusterpool, subgroup)

        PoolPopupWindow(self, data['clusterpools'][clusterpool][subgroup])          # open popup window that handles pool editing

    '''
    updateData
    - updates 'data' dictionary before export
    - clusterpool and subgroup info is updated automatically over the course of the program, but currently the project name and features are not and
      have to be updated before we export the data as a JSON
    '''
    def updateData(self):
        data['project_name'] = self.projectNameEntry.get()
        for iid in self.featureTreeView.get_children():                         # get_children returns the item id (iid) of every feature in the Treeview
            data['features'].append(self.featureTreeView.item(iid)['text'])

    '''
    export
    - exports the user-defined data as a JSON to be read by R scripts on the back end
    - this file is named using the project name defined by the user
    '''
    def export(self, data):
        self.updateData()
        print('--exporting data to Data/JSON/' + data['project_name'] + '.json')
        print(data)
        with open('Data/JSON/' + data['project_name'] + '.json', 'w') as outfile:
            json.dump(data, outfile)
        #run the scripts with project name :)
        runScript(data['project_name'])
        


'''
PoolPopupWindow class
- a popup window that allows the user to edit the clusters included within a specific pool
- created when the user clicks the editPoolsButton in the main window
'''
class PoolPopupWindow(tk.Toplevel):
    '''
    init
    - function called when an instance of the class is created
    - creates and places the widgets within the popup window
    '''
    def __init__(self, parent, poolList):
        tk.Toplevel.__init__(self)
        self.parent = parent
        self.title('Pool Editor')

        self.excCheckboxList = []                                           # list that stores the excitatory Checkbutton widgets
        self.excCheckboxVariables = {}                                      # dictionary that pairs a cluster name (e.g. Excit-1) with an IntVar, a ...
                                                                            # ... special integer object (value 0 or 1) used by Tkinter's Checkbutton class
        self.inhCheckboxList = []                                           # same idea but for inhibitory Checkbuttons
        self.inhCheckboxVariables = {}

        for i in range(38):                                                 # create excitatory checkboxes Excit-1 to Excit-38
            name = "Excit-"+str(i+1)
            self.excCheckboxVariables[name] = tk.IntVar(self, 0, name)      # add the checkbox to excCheckboxVariables
            self.excCheckboxList.append(tk.Checkbutton(                     # add the checkbox to excCheckboxList
                self, 
                text=name,
                variable=self.excCheckboxVariables[name]
            ))
            self.excCheckboxList[i].grid(                                   # place the checkbox in the window
                row = i if i<19 else i-19,                                  # these if statements allow for screen wrap (columns of 19)
                column = 0 if i<19 else 1)
            if name in poolList:                                            # if the cluster is already in the pool, the checkbox should be selected,
                self.excCheckboxList[i].select()                            # otherwise it should be deselected
            else:
                self.excCheckboxList[i].deselect()

        for i in range(27):                                                 # same as above but for inhibitory checkboxes
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

    '''
    openCopyPoolsWindow
    - opens a popup window that allows the user to copy the clusters from another pool
    - this function is called when the user clicks the copyPoolsButton in the PoolPopupWindow
    '''
    def openCopyPoolsWindow(self):
        CopyPoolsWindow(self)

    '''
    savePoolListAndClose
    - saves the user's selected clusters to the 'data' dictionary via the poolList
    - closes the PoolPopupWindow
    '''
    def savePoolListAndClose(self, poolList):
        poolList.clear()                                            # empty the poolList
        for key in self.excCheckboxVariables:                       # add any selected clusters to the poolList
            # print(self.excCheckboxVariables[key])
            if self.excCheckboxVariables[key].get() == 1:
                poolList.append(key)
        for key in self.inhCheckboxVariables:
            # print(self.inhCheckboxVariables[key])
            if self.inhCheckboxVariables[key].get() == 1:
                poolList.append(key)
        print(poolList)
        self.destroy()                                              # close the popup window


'''
CopyPoolsWindow class
- a popup window that allows the user to select a pool to copy clusters from
- created when the user clicks the copyPoolsButton in the PoolPopupWindow
'''
class CopyPoolsWindow(tk.Toplevel):
    '''
    init
    - function called when an instance of the class is created
    - creates and places the widgets within the popup window
    '''
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

        for clusterpool in data['clusterpool_names']:                           # populates the copyPoolsTreeView with all combinations of clusterpool and subgroup
            for subgroup in data['subgroup_names']:                             # the user can then select one of these options to copy the clusters from
                name = clusterpool + ' // ' + subgroup
                self.copyPoolsTreeView.insert('', 'end', text=name)
        
        self.selectButton = ttk.Button(self, text="Select", command=self.copyPools)
        self.selectButton.grid(row=2, column=1)

    '''
    copyPools
    - copies the clusters of the selected subgroup
    - closes the CopyPoolsWindow
    '''
    def copyPools(self):
        selectedItem = self.copyPoolsTreeView.focus()                           # get the user's selected item
        selectedValue = self.copyPoolsTreeView.item(selectedItem)['text']       # get the text of the selected item
        selectedValue = selectedValue.split(' // ')                             # we split this text by the '//' delimeter to get the clusterpool and subgroup
        clusterpool = selectedValue[0]
        subgroup = selectedValue[1]

        poolList = data['clusterpools'][clusterpool][subgroup]                  # the list of clusters in the subgroup the user chose

        for i in range(38):                                                     # for every cluster, if the poolList contains that cluster we will select the...
            name = "Excit-"+str(i+1)                                            # ... respective box in the PoolPopupWindow -- otherwise, we will deselect it
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

        self.destroy()                                                          # close the CopyPoolsWindow

def runScript(project_name):
    #run all R scripts
    #command line call

    #pass project name (JSON data file name) to R
    return

'''
main
- defines a dictionary 'data' which is used to store the user input and can be exported as a JSON
- initiates an instance of the App class
'''
if __name__ == "__main__":
    data = {
        'features': [],
        'project_name': '',
        'clusterpool_names': [],
        'subgroup_names': [],
        'clusterpools': {}
    }

    os.chdir(os.path.dirname(sys.argv[0]))                  # ??? something important Justin added lol
    
    root = tk.Tk()                                          # the root of all GUI windows and widgets in tkinter
    root.title("Cluster Pool Comparison Tool")

    root.tk.call("source", "guitheme/azure.tcl")            # applies the custom GUI theme Azure (only works for ttk widgets, not tk widgets)
    root.tk.call("set_theme", "light")

    app = App(root, data)
    app.pack(fill="both", expand=True)

    root.mainloop()                                         # begins loop that allows the GUI to update