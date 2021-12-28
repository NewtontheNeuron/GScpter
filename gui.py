import json

data = {}

data['features'] = []
data['project_name'] = ''
data['filenames'] = {}
data['clusterpools'] = {}

with open('test.json', 'w') as outfile:
    json.dump(data, outfile)