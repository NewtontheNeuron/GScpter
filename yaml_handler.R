# New library
library(yaml)

yaml_file <- paste(getwd(), '/example.yaml', sep='')

yaml_obj <- read_yaml(yaml_file, fileEncoding = "UTF-8")

print(yaml_obj)