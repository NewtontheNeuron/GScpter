<<<<<<< HEAD
import json

data = {}

data['features'] = []
data['project_name'] = ''
data['filenames'] = {}
data['clusterpools'] = {}

data['features'] = ["rna_Grin1", "rna_Grin2a", "rna_Grin2b", "rna_Grin2c", "rna_Grin2d", "rna_Grin3a", "rna_Grin3b"]
data['project_name'] = "SDH_vs_DDH"
data['filenames'] = {'pooled':'x', 'unpooled':'z', 'barplot':'h'}
data['clusterpools'] = {
    'excitatory': {
        'SDH': ["Excit-01", "Excit-02","Excit-03","Excit-08","Excit-09",
                "Excit-10","Excit-12","Excit-14","Excit-15","Excit-16",
                "Excit-18","Excit-04","Excit-05","Excit-13","Excit-19"],
        'DDH': ["Excit-5", "Excit-6","Excit-20","Excit-21","Excit-22",
                "Excit-23","Excit-24","Excit-25","Excit-26","Excit-27",
                "Excit-29","Excit-30","Excit-31","Excit-32","Excit-34",
                "Excit-35","Excit-36"]
    },
    'inhibitory': {
        'SDH': ["Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05",
                "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10", "Inhib-11",
                "Inhib-12", "Inhib-13"],
        'DDH': ["Inhib-3", "Inhib-6", "Inhib-8", "Inhib-12", "Inhib-14",
                "Inhib-15", "Inhib-16", "Inhib-18", "Inhib-19", "Inhib-20",
                "Inhib-21"]
    }
}

with open('example.json', 'w') as outfile:
    json.dump(data, outfile)
