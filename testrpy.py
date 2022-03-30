# Set the project name
project_name = "data"
# R function load the data
r.load_bin("../../")
# Get some metadata from the RDSfile
meta = r.get_meta_data()
# Access metadata from python
print(meta)
# R function send the data to the queue
r.add_to_queue(project_name)
