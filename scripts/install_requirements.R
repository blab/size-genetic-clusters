args <- commandArgs(trailingOnly = T)
requirement_file_path <- as.character(args[1])

# Load list of required libraries
name_packages_to_download <- as.vector(read.table(requirement_file_path, header = F))$V1

# Check packages that need to be installed
vec_packages_already_installed <- installed.packages()[,1]
vec_packages_to_install <- name_packages_to_download[! name_packages_to_download %in% vec_packages_already_installed]

# Install required packages
install.packages(vec_packages_to_install)
