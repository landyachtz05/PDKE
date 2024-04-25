print("test")

# Access the arguments passed from the server script

# Get the command-line arguments passed from the main script
args <- commandArgs(trailingOnly = TRUE)
print(paste("args:", length(args), "\n"))

for (a in args) {
  print(paste("arg_v2: ", a))
}
