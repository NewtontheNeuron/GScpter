# New library
library(yaml)

# Stores the absolute path to the example.yaml 
yaml_file <- paste(getwd(), '/example.yaml', sep='')

# Loads the yaml but there is an issue because of the part
# where there are multiple Exciatatory keys in the same 
# clusterpools section.
yaml_obj <- read_yaml(yaml_file, fileEncoding = "UTF-8")

print(yaml_obj)
